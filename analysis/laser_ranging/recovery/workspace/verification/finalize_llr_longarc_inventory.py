"""Retain every fixed-interval row and describe observable metadata, without native model calls."""

from pathlib import Path
from collections import Counter, defaultdict
from decimal import Decimal
from datetime import date, timedelta
import csv, json, hashlib

WORK = Path("/home/dominic/Tudat/laser-ranging-work")
OUT = WORK / "research/longarc-20170304"
reflectors = {"apollo11": 0, "luna17": 1, "apollo14": 2, "apollo15": 3, "luna21": 4}


def label(day, seconds):
    """Preserve decimal transmit seconds in an unambiguous UTC calendar label."""
    hour = int(seconds // 3600)
    minute = int((seconds - hour * 3600) // 60)
    second = seconds - hour * 3600 - minute * 60
    return f"{day.isoformat()}T{hour:02}:{minute:02}:{second:015.12f}"


rows = []
sessionid = 0
epoch_counts = Counter()
flags = Counter()
wavelengths = Counter()
meteo_gaps = []
for path in sorted((OUT / "crd").glob("*.npt")):
    session = None
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        fields = line.split()
        if not fields:
            continue
        code = fields[0].lower()
        if code == "h2":
            station = fields[1:]
        elif code == "h3":
            target = fields[1]
        elif code == "h4":
            sessionid += 1
            session = {
                "id": sessionid,
                "day": date(*map(int, fields[2:5])),
                "start_s": int(fields[5]) * 3600 + int(fields[6]) * 60 + Decimal(fields[7]),
                "h4": fields[14:22],
                "meteo": [],
                "normal": [],
                "calibrations": 0,
                "wavelength": None,
            }
        elif code == "c0":
            session["wavelength"] = fields[2]
        elif code == "20":
            session["meteo"].append([float(v) for v in fields[1:6]])
        elif code in ["40", "41"]:
            session["calibrations"] += 1
        elif code == "11":
            session["normal"].append((lineno, fields))
        elif code == "h8":
            assert session is not None and session["meteo"] and session["calibrations"]
            flags[tuple(session["h4"])] += 1
            for lineno, normal in session["normal"]:
                seconds = Decimal(normal[1])
                day = session["day"]
                if seconds < session["start_s"] - 43200:
                    day += timedelta(days=1)
                gaps = [abs(float(seconds) - m[0]) for m in session["meteo"]]
                nearest = min(gaps)
                meteo_gaps.append(nearest)
                epoch_counts[normal[4]] += 1
                wavelengths[session["wavelength"]] += 1
                rows.append(
                    {
                        "source_file": path.name,
                        "source_line": lineno,
                        "session_id": sessionid,
                        "station_cdp": station[1],
                        "system": station[2],
                        "occupancy": station[3],
                        "utc_reference_code": station[4],
                        "reflector_code": reflectors[target],
                        "target": target,
                        "utc_launch": label(day, seconds),
                        "seconds_of_day_printed": normal[1],
                        "two_way_time_of_flight_s": normal[2],
                        "epoch_event": normal[4],
                        "wavelength_nm": session["wavelength"],
                        "nearest_weather_s": nearest,
                        "calibration_records_in_session": session["calibrations"],
                        "bin_rms_ps_metadata_only": normal[7],
                    }
                )
            session = None
assert len(rows) == 237
rows.sort(key=lambda r: (r["utc_launch"], r["reflector_code"], r["source_file"], r["source_line"]))
with (OUT / "grasse_exact_inventory.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
a = Counter()
with (OUT / "apollo_night_reflector_counts.csv").open() as f:
    for r in csv.DictReader(f):
        a[(r["utc_night"], int(r["reflector_code"]))] = int(r["normal_points"])
g = Counter((r["utc_launch"][:10], r["reflector_code"]) for r in rows)
with (OUT / "cross_station_night_reflector_overlap.csv").open("w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["utc_night", "reflector_code", "apollo_2022_group_e", "grasse_crdv1"])
    w.writerows(
        (day, reflector, a[(day, reflector)], g[(day, reflector)])
        for day, reflector in sorted(a.keys() | g.keys())
    )
summary = json.loads((OUT / "coverage-summary.json").read_text())
summary["edc_current_crd"].update(
    nights=len({r["utc_launch"][:10] for r in rows}),
    epoch_event_counts=dict(epoch_counts),
    wavelength_point_counts=dict(wavelengths),
    header_flags={
        "data_release": 0,
        "troposphere_applied": 0,
        "center_of_mass_applied": 0,
        "receive_amplitude_applied": 0,
        "station_system_delay_applied": 1,
        "spacecraft_system_delay_applied": 0,
        "range_type": 2,
        "quality_alert": 0,
    },
    meteo_missing_sessions=0,
    calibration_missing_sessions=0,
    max_nearest_weather_seconds=max(meteo_gaps),
    min_nearest_weather_seconds=min(meteo_gaps),
)
shared = [key for key in a if key in g]
summary["overlap"] = {
    "utc_nights": sorted({day for day, reflector in shared}),
    "night_reflector_pairs": len(shared),
    "apollo_points_on_shared_night_reflector_pairs": sum(a[k] for k in shared),
    "grasse_points_on_shared_night_reflector_pairs": sum(g[k] for k in shared),
}
summary["release_caution"] = (
    "EDC files are existing CRDv1 archive products, not the frozen APOLLO Group-e2022 release. No replacement/deduplication across the two stations or provider releases."
)
summary["matera"] = (
    "No Matera in the ten inspected CRDv1 monthly products; separate V2 browser query returnedHTTP500 and remains unchecked. Do not infer general absence."
)
summary["sigma"] = (
    "APOLLO provider uncertainty retained; CRD bin RMS is accepted-return spread only, no NP sigma or sqrt(N) inference."
)
(OUT / "coverage-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps({k: v for k, v in summary.items() if k != "apollo_release"}, indent=2))
