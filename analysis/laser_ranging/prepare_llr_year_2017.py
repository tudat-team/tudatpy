"""Freeze the calendar-2017 source subset and CRD inventory without computing residuals."""

from pathlib import Path
from collections import Counter
from decimal import Decimal
from datetime import date, timedelta
import argparse
import csv
import json
import hashlib

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument("--data-root", type=Path, required=True)
parser.add_argument("--apollo-release", type=Path, required=True)
parser.add_argument(
    "--manifest", type=Path, default=Path(__file__).with_name("llr_year_2017_dataset.json")
)

args = parser.parse_args()

ROOT = args.data_root

SOURCE = args.apollo_release

manifest = json.loads(args.manifest.read_text())

# Verify every unchanged provider product before deriving source metadata.
for product in manifest["grasse"]["files"]:
    if hashlib.sha256((ROOT / product["file"]).read_bytes()).hexdigest() != product["sha256"]:
        raise ValueError(f"Provider input changed: {product['file']}")

if hashlib.sha256(SOURCE.read_bytes()).hexdigest() != manifest["apollo"]["source_release_sha256"]:
    raise ValueError("APOLLO provider release changed")
selected = [
    (i, line)
    for i, line in enumerate(SOURCE.read_bytes().splitlines(keepends=True), 1)
    if line.startswith(b"512017")
]

subset = ROOT / "apollo_group_e_20170101_20171231.txt"

subset_bytes = b"".join(line for _, line in selected)

if (
    hashlib.sha256(subset_bytes).hexdigest() != manifest["apollo"]["sha256"]
    or [index for index, _ in selected] != manifest["apollo"]["source_line_numbers"]
):
    raise ValueError("APOLLO subset differs from the frozen source rows")
subset.write_bytes(subset_bytes)
(ROOT / "apollo_subset_sources.json").write_text(
    json.dumps(
        {
            "file": subset.name,
            "sha256": hashlib.sha256(subset.read_bytes()).hexdigest(),
            "count": len(selected),
            "source_line_numbers": [i for i, _ in selected],
        },
        indent=2,
    )
    + "\n"
)
reflectors = {"apollo11": 0, "luna17": 1, "apollo14": 2, "apollo15": 3, "luna21": 4}


def utc_label(day, seconds):
    """Keep the printed decimal seconds in a UTC calendar label."""
    h = int(seconds // 3600)

    m = int((seconds - h * 3600) // 60)

    s = seconds - h * 3600 - m * 60

    return f"{day.isoformat()}T{h:02}:{m:02}:{s:015.12f}"


# Decode printed fixed source records, preserving their original line identities.
rows = []

sessions = []

sid = 0

for path in sorted(ROOT / product["file"] for product in manifest["grasse"]["files"]):
    session = None

    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        f = line.split()

        if not f:
            continue
        tag = f[0].lower()

        if tag == "h2":
            station = f[1:]

        elif tag == "h3":
            target = f[1]

        elif tag == "h4":
            sid += 1

            session = dict(
                id=sid,
                day=date(*map(int, f[2:5])),
                start=Decimal(f[5]) * 3600 + Decimal(f[6]) * 60 + Decimal(f[7]),
                h4=f[14:22],
                normal=[],
                weather=[],
                calibrations=0,
                configs={},
            )

        elif tag == "c0":
            session["configs"][f[3]] = f[2]

        elif tag == "20":
            session["weather"].append(f[1:6])
        elif tag in ["40", "41"]:
            session["calibrations"] += 1

        elif tag == "11":
            session["normal"].append((line_number, f))
        elif tag == "h8":
            day = session["day"]

            previous = None

            normal_times = []

            for line_number, np in session["normal"]:
                seconds = Decimal(np[1])

                if previous is not None and seconds + 10000 < previous:
                    day += timedelta(days=1)

                normal_times.append((line_number, np, day, seconds))
                previous = seconds

            weather = []

            met_day = session["day"]

            previous = None

            for met in session["weather"]:
                seconds = Decimal(met[0])

                if previous is not None and seconds + 10000 < previous:
                    met_day += timedelta(days=1)

                weather.append((met_day, seconds, met))
                previous = seconds

            for line_number, np, day, seconds in normal_times:
                absolute = (day - date(2000, 1, 1)).days * 86400 + seconds

                met_times = [(d - date(2000, 1, 1)).days * 86400 + s for d, s, _ in weather]

                nearest = (
                    min(abs(float(absolute - epoch)) for epoch in met_times) if met_times else None
                )

                outside = (
                    None
                    if not met_times
                    else max(
                        float(min(met_times) - absolute), float(absolute + 5 - max(met_times)), 0.0
                    )
                )

                excluded = (
                    path.name == manifest["excluded_normal_point"]["source_file"]
                    and line_number == manifest["excluded_normal_point"]["source_line"]
                )

                if (
                    excluded
                    and utc_label(day, seconds) != manifest["excluded_normal_point"]["utc_launch"]
                ):
                    raise ValueError("Excluded source identity changed")
                row = dict(
                    source_file=path.name,
                    source_line=line_number,
                    session_id=sid,
                    station_cdp=station[1],
                    system=station[2],
                    occupancy=station[3],
                    utc_reference_code=station[4],
                    reflector_code=reflectors[target],
                    target=target,
                    utc_launch=utc_label(day, seconds),
                    seconds_of_day_printed=np[1],
                    two_way_time_of_flight_s=np[2],
                    epoch_event=np[4],
                    wavelength_nm=session["configs"][np[3]],
                    nearest_weather_s=nearest,
                    calibration_records_in_session=session["calibrations"],
                    bin_rms_ps_metadata_only=np[7],
                    weather_endpoint_hold_required_s=outside,
                    retained=not excluded,
                )

                rows.append(row)
            sessions.append(
                dict(
                    file=path.name,
                    session_id=sid,
                    start_date=session["day"].isoformat(),
                    start_s=str(session["start"]),
                    station=station,
                    target=target,
                    h4=session["h4"],
                    normal_points=len(normal_times),
                    weather_records=session["weather"],
                    calibration_records=session["calibrations"],
                    configurations=session["configs"],
                    normal_rollovers=sum(day != session["day"] for _, _, day, _ in normal_times),
                    first_np_before_h4_start_by12h=bool(
                        normal_times and normal_times[0][3] + 43200 < session["start"]
                    ),
                )
            )
            session = None

all_rows = rows

# Select the predeclared Grasse station; other stations remain in the session inventory.
rows = [row for row in all_rows if row["station_cdp"] == "7845"]

if (
    len(rows) != manifest["grasse"]["source_count"]
    or sum(row["retained"] for row in rows) != manifest["grasse"]["retained_count"]
):
    raise ValueError("Grasse inventory differs from the frozen count")

rows.sort(key=lambda r: (r["utc_launch"], r["reflector_code"], r["source_file"], r["source_line"]))
with (ROOT / "grasse_exact_inventory.csv").open("w", newline="") as stream:
    w = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")

    w.writeheader()
    w.writerows(rows)
(ROOT / "crd_sessions.json").write_text(json.dumps(sessions, indent=2) + "\n")
summary = dict(
    apollo_count=len(selected),
    crd_raw_count=len(all_rows),
    grasse_count=sum(r["station_cdp"] == "7845" for r in rows),
    retained_count=sum(r["retained"] and r["station_cdp"] == "7845" for r in rows),
    station_counts=dict(Counter(r["station_cdp"] for r in all_rows)),
    sessions=len(sessions),
    h4_flags={str(k): v for k, v in Counter(tuple(s["h4"]) for s in sessions).items()},
    epoch_events=dict(Counter(r["epoch_event"] for r in rows)),
    wavelengths=dict(Counter(r["wavelength_nm"] for r in rows)),
    weather_origins=dict(Counter(m[4] for s in sessions for m in s["weather_records"])),
    normal_rollover_count=sum(s["normal_rollovers"] for s in sessions),
    first_np_date_ambiguities=[
        s["session_id"] for s in sessions if s["first_np_before_h4_start_by12h"]
    ],
    max_nearest_weather_s=max(
        r["nearest_weather_s"] for r in rows if r["nearest_weather_s"] is not None
    ),
    weather_hold_over_1000=[
        r
        for r in rows
        if r["weather_endpoint_hold_required_s"] is None
        or r["weather_endpoint_hold_required_s"] > 1000
    ],
    first_grasse_utc=rows[0]["utc_launch"],
    last_grasse_utc=rows[-1]["utc_launch"],
)

(ROOT / "coverage-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
