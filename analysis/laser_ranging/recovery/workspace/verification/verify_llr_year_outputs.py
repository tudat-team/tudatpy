"""Independently audit saved annual residual CSVs; no native imports or solves."""

import argparse
import collections
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(row):
    return (
        row["station"],
        "APOLLO-original-release" if row["station"] == "APOLLO" else row["source_file"],
        int(row["source_line"]),
        row["utc_launch"],
    )


def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--annual", type=Path, required=True)
parser.add_argument("--baseline", type=Path, required=True)
parser.add_argument("--manifest", type=Path, required=True)
parser.add_argument("--data-root", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
manifest = json.loads(args.manifest.read_text())
expected = []
for source in read(args.data_root / manifest["grasse"]["inventory"]["file"]):
    exclusions = [manifest["excluded_normal_point"]] + manifest.get(
        "source_invalid_normal_points", []
    )
    matched = [
        record
        for record in exclusions
        if int(source["source_line"]) == record["source_line"]
        and source["source_file"] == record["source_file"]
    ]
    if matched:
        assert len(matched) == 1 and source["utc_launch"] == matched[0]["utc_launch"]
        continue
    source["station"] = "GRASSE"
    expected.append(identity(source))
widths = (2, 4, 2, 2, 2, 2, 9, 14, 1, 5, 3, 6, 3, 1, 6, 4, 2, 5, 1, 4, 1)
for number, line in zip(
    manifest["apollo"]["source_line_numbers"],
    (args.data_root / manifest["apollo"]["file"]).read_text().splitlines(),
):
    values = []
    offset = 0
    for width in widths:
        values.append(line[offset : offset + width])
        offset += width
    year, month, day, hour, minute = map(int, values[1:6])
    seconds = Decimal(values[6]) * Decimal("1e-7")
    utc = f"{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{seconds:010.7f}"
    expected.append(("APOLLO", "APOLLO-original-release", number, utc))
annual = read(args.annual)
baseline = read(args.baseline)
actual = collections.Counter(identity(row) for row in annual)
wanted = collections.Counter(expected)
annual_map = {identity(row): row for row in annual}
overlap = [row for row in annual if "2017-03-01" <= row["utc_launch"] < "2017-05-01"]
changes = []
for row in baseline:
    key = identity(row)
    new = annual_map.get(key)
    if new is None:
        changes.append({"identity": key, "missing": True})
    elif float(row["computed_full_two_way_range_m"]) != float(new["computed_full_two_way_range_m"]):
        changes.append(
            {
                "identity": key,
                "old_m": row["computed_full_two_way_range_m"],
                "new_m": new["computed_full_two_way_range_m"],
                "change_m": float(new["computed_full_two_way_range_m"])
                - float(row["computed_full_two_way_range_m"]),
            }
        )
weather_rows = [row for row in annual if row["weather_mode"] == "held-last-long-source-limitation"]
weather_expected = {
    (q["source_file"], q["source_line"], q["utc_launch"])
    for q in manifest["grasse"]["weather_endpoint_hold_exceptions"]
}
weather_actual = {(q["source_file"], int(q["source_line"]), q["utc_launch"]) for q in weather_rows}
duplicates = [
    row
    for row in annual
    if row["source_file"] == "luna21_201708.npt" and int(row["source_line"]) in [102, 105]
]
station_counts = collections.Counter(row["station"] for row in annual)
coverage = {
    station: {
        "first_utc": min(row["utc_launch"] for row in annual if row["station"] == station),
        "last_utc": max(row["utc_launch"] for row in annual if row["station"] == station),
        "utc_nights": len({row["utc_launch"][:10] for row in annual if row["station"] == station}),
    }
    for station in station_counts
}
report = dict(
    annual_csv_sha256=digest(args.annual),
    baseline_csv_sha256=digest(args.baseline),
    manifest_sha256=digest(args.manifest),
    row_count=len(annual),
    station_counts=dict(station_counts),
    missing_identities=list((wanted - actual).elements()),
    unexpected_identities=list((actual - wanted).elements()),
    duplicate_output_identity_count=sum(n - 1 for n in actual.values()),
    requested_exclusion_present=any(
        row["station"] == "GRASSE"
        and row["source_file"] == "apollo15_201703.npt"
        and int(row["source_line"]) == 171
        for row in annual
    ),
    repeated_provider_record_source_lines=[int(row["source_line"]) for row in duplicates],
    long_weather_rows=[
        {
            k: row[k]
            for k in [
                "source_file",
                "source_line",
                "utc_launch",
                "nearest_measured_weather_s",
                "weather_mode",
            ]
        }
        for row in weather_rows
    ],
    long_weather_exact_match=weather_actual == weather_expected and len(weather_rows) == 3,
    coverage=coverage,
    march_april_count=len(overlap),
    baseline_count=len(baseline),
    march_april_identity_exact_match=collections.Counter(identity(row) for row in overlap)
    == collections.Counter(identity(row) for row in baseline),
    changed_computed_ranges=changes,
    all_epoch_events_transmit=all(row["epoch_event"] == "2" for row in annual),
)
report["coverage_matches_frozen_sources"] = all(
    coverage[station]["first_utc"] == manifest[key]["first_utc_launch"]
    and coverage[station]["last_utc"] == manifest[key]["last_utc_launch"]
    and coverage[station]["utc_nights"] == manifest[key]["nights"]
    for station, key in [("APOLLO", "apollo"), ("GRASSE", "grasse")]
)
report["source_invalid_rows_present"] = [
    identity(row)
    for row in annual
    for record in manifest.get("source_invalid_normal_points", [])
    if row["station"] == record["station"]
    and row["source_file"] == record["source_file"]
    and int(row["source_line"]) == record["source_line"]
]
report["pass"] = (
    not report["source_invalid_rows_present"]
    and report["coverage_matches_frozen_sources"]
    and len(annual) == 1579
    and station_counts == {"APOLLO": 311, "GRASSE": 1268}
    and actual == wanted
    and not report["requested_exclusion_present"]
    and sorted(report["repeated_provider_record_source_lines"]) == [102, 105]
    and report["long_weather_exact_match"]
    and report["march_april_identity_exact_match"]
    and len(overlap) == 305
    and not changes
    and report["all_epoch_events_transmit"]
)
args.output.write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
raise SystemExit(0 if report["pass"] else 1)
