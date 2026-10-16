"""Describe raw annual LLR month/night structure without fitting or changing any residual."""

import csv
from collections import defaultdict
from datetime import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

source, output = map(Path, sys.argv[1:3])
rows = list(csv.DictReader(source.open(newline="")))
output.mkdir(parents=True, exist_ok=False)


def write_rows(name, values):
    """Retain plain descriptive CSV columns."""
    with (output / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


monthly = []
for station in ("APOLLO", "GRASSE"):
    for reflector in ("ALL", "A11", "A14", "A15", "L1", "L2"):
        for month in range(1, 13):
            key = f"2017-{month:02}"
            group = [
                row
                for row in rows
                if row["station"] == station
                and row["utc_launch"].startswith(key)
                and (reflector == "ALL" or row["reflector"] == reflector)
            ]
            values = np.array([float(row["one_way_residual_m"]) for row in group])
            monthly.append(
                dict(
                    utc_month=key,
                    station=station,
                    reflector=reflector,
                    count=len(values),
                    utc_nights=len({row["utc_launch"][:10] for row in group}),
                    raw_median_m=float(np.median(values)) if len(values) else None,
                    raw_percentile10_m=float(np.percentile(values, 10)) if len(values) else None,
                    raw_percentile90_m=float(np.percentile(values, 90)) if len(values) else None,
                    long_held_weather_count=sum(
                        row["weather_mode"] == "held-last-long-source-limitation" for row in group
                    ),
                )
            )
write_rows("monthly_raw_medians.csv", monthly)

groups = defaultdict(list)
for row in rows:
    groups[(row["station"], row["reflector"], row["utc_launch"][:10])].append(row)
nightly = []
for (station, reflector, night), group in sorted(groups.items()):
    values = np.array([float(row["one_way_residual_m"]) for row in group])
    nightly.append(
        dict(
            station=station,
            reflector=reflector,
            utc_night=night,
            count=len(values),
            first_utc=min(row["utc_launch"] for row in group),
            last_utc=max(row["utc_launch"] for row in group),
            raw_mean_m=float(values.mean()),
            raw_median_m=float(np.median(values)),
            raw_minimum_m=float(values.min()),
            raw_maximum_m=float(values.max()),
        )
    )
write_rows("nightly_raw_means_and_medians.csv", nightly)

extrema = []
for station, reflector in sorted({(row["station"], row["reflector"]) for row in nightly}):
    group = [row for row in nightly if row["station"] == station and row["reflector"] == reflector]
    low = min(group, key=lambda row: row["raw_mean_m"])
    high = max(group, key=lambda row: row["raw_mean_m"])
    extrema.append(
        dict(
            station=station,
            reflector=reflector,
            utc_nights=len(group),
            minimum_raw_nightly_mean_m=low["raw_mean_m"],
            minimum_utc_night=low["utc_night"],
            minimum_count=low["count"],
            maximum_raw_nightly_mean_m=high["raw_mean_m"],
            maximum_utc_night=high["utc_night"],
            maximum_count=high["count"],
            raw_nightly_mean_peak_to_peak_m=high["raw_mean_m"] - low["raw_mean_m"],
        )
    )
write_rows("dated_nightly_extrema.csv", extrema)

pairing = Path(__file__).with_name("compare_llr_nearby_nights.py")
subprocess.run([sys.executable, str(pairing), str(source), str(output)], check=True)
pairs = list(csv.DictReader((output / "nearby_nightly_raw_means.csv").open(newline="")))
paired = []
for reflector in ("ALL", "A11", "A14", "A15", "L1", "L2"):
    group = [row for row in pairs if reflector == "ALL" or row["reflector"] == reflector]
    values = np.array([float(row["grasse_minus_apollo_raw_mean_m"]) for row in group])
    gaps = np.array([float(row["grasse_minus_apollo_hours"]) for row in group])
    paired.append(
        dict(
            reflector=reflector,
            paired_nights=len(group),
            median_raw_grasse_minus_apollo_m=float(np.median(values)),
            percentile10_raw_difference_m=float(np.percentile(values, 10)),
            percentile90_raw_difference_m=float(np.percentile(values, 90)),
            minimum_raw_difference_m=float(values.min()),
            maximum_raw_difference_m=float(values.max()),
            minimum_absolute_gap_hours=float(abs(gaps).min()),
            maximum_absolute_gap_hours=float(abs(gaps).max()),
        )
    )
write_rows("paired_station_raw_summary.csv", paired)
summary = dict(
    source_csv_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    source_count=len(rows),
    generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    pairing_source_sha256=hashlib.sha256(pairing.read_bytes()).hexdigest(),
    monthly="Raw unweighted NP median with N/UTC-night count and P10/P90 for each calendar month; ALL pools reflector composition; empty months remain empty",
    nightly="Raw unweighted mean/median, with dated extrema and no correction or centering",
    pairing="Nearest same-reflector Grasse nightly mean epoch to each APOLLO night, absolute gap<=12h; exact signed gaps retained, no temporal alignment/interpolation",
    caution="Uneven lunar phase, wavelength, reflector and station sampling can change monthly medians. Paired observations are distinct epochs, not simultaneous station biases. No residual adjusted, fitted or removed.",
    observation_solves=0,
)
(output / "temporal_provenance.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(paired, indent=2))
