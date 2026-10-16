"""Describe raw calendar-month/night LLR structure and nearby station pairs without fitting."""

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np

from plot_llr_residuals import read_points


def write_csv(path, rows):
    """Write descriptive fields without altering the source residual CSV."""

    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    """Summarize one completed calendar-2017 case; never construct or solve an observation."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    points = read_points(args.csv)
    args.output.mkdir(parents=True, exist_ok=False)
    monthly = []

    for station in ("APOLLO", "GRASSE"):
        for reflector in ("ALL", "A11", "A14", "A15", "L1", "L2"):
            for month in range(1, 13):
                group = [
                    p
                    for p in points
                    if p["station"] == station
                    and p["epoch"].year == 2017
                    and p["epoch"].month == month
                    and (reflector == "ALL" or p["reflector"] == reflector)
                ]
                values = np.array([p["residual"] for p in group])
                monthly.append(
                    dict(
                        utc_month=f"2017-{month:02}",
                        station=station,
                        reflector=reflector,
                        count=len(group),
                        utc_nights=len({p["epoch"].date() for p in group}),
                        raw_median_m=float(np.median(values)) if len(group) else None,
                        raw_percentile10_m=float(np.percentile(values, 10)) if len(group) else None,
                        raw_percentile90_m=float(np.percentile(values, 90)) if len(group) else None,
                    )
                )
    write_csv(args.output / "monthly_raw_medians.csv", monthly)
    groups = defaultdict(list)

    for point in points:
        groups[(point["station"], point["reflector"], point["epoch"].date().isoformat())].append(
            point
        )
    nightly = [
        dict(
            station=station,
            reflector=reflector,
            utc_night=night,
            count=len(group),
            mean_epoch_s=sum(p["epoch"].timestamp() for p in group) / len(group),
            raw_mean_m=sum(p["residual"] for p in group) / len(group),
        )
        for (station, reflector, night), group in sorted(groups.items())
    ]
    extrema = []

    for station, reflector in sorted({(r["station"], r["reflector"]) for r in nightly}):
        group = [r for r in nightly if (r["station"], r["reflector"]) == (station, reflector)]
        low = min(group, key=lambda r: r["raw_mean_m"])
        high = max(group, key=lambda r: r["raw_mean_m"])
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
    write_csv(args.output / "dated_nightly_extrema.csv", extrema)
    pairs, unmatched, used = [], [], set()

    for apollo in [r for r in nightly if r["station"] == "APOLLO"]:
        candidates = [
            r for r in nightly if r["station"] == "GRASSE" and r["reflector"] == apollo["reflector"]
        ]
        grasse = min(candidates, key=lambda r: abs(r["mean_epoch_s"] - apollo["mean_epoch_s"]))
        gap = (grasse["mean_epoch_s"] - apollo["mean_epoch_s"]) / 3600
        if abs(gap) > 12:
            continue
        used.add((grasse["reflector"], grasse["utc_night"]))
        pairs.append(
            dict(
                reflector=apollo["reflector"],
                apollo_night=apollo["utc_night"],
                grasse_night=grasse["utc_night"],
                apollo_mean_utc=datetime.fromtimestamp(
                    apollo["mean_epoch_s"], timezone.utc
                ).isoformat(),
                grasse_mean_utc=datetime.fromtimestamp(
                    grasse["mean_epoch_s"], timezone.utc
                ).isoformat(),
                grasse_minus_apollo_hours=gap,
                apollo_count=apollo["count"],
                grasse_count=grasse["count"],
                apollo_raw_mean_m=apollo["raw_mean_m"],
                grasse_raw_mean_m=grasse["raw_mean_m"],
                grasse_minus_apollo_raw_mean_m=grasse["raw_mean_m"] - apollo["raw_mean_m"],
            )
        )
    write_csv(args.output / "nearby_nightly_raw_means.csv", pairs)
    matched_apollo = {(r["reflector"], r["apollo_night"]) for r in pairs}

    for row in nightly:
        matched = matched_apollo if row["station"] == "APOLLO" else used
        if (row["reflector"], row["utc_night"]) in matched:
            continue
        other = [
            r
            for r in nightly
            if r["station"] != row["station"] and r["reflector"] == row["reflector"]
        ]
        nearest = min(other, key=lambda r: abs(r["mean_epoch_s"] - row["mean_epoch_s"]))
        unmatched.append(
            dict(
                station=row["station"],
                reflector=row["reflector"],
                utc_night=row["utc_night"],
                count=row["count"],
                nearest_other_station_night=nearest["utc_night"],
                nearest_other_minus_this_hours=(nearest["mean_epoch_s"] - row["mean_epoch_s"])
                / 3600,
            )
        )
    write_csv(args.output / "unmatched_station_reflector_nights.csv", unmatched)
    paired = []

    for reflector in ("ALL", "A11", "A14", "A15", "L1", "L2"):
        group = [r for r in pairs if reflector == "ALL" or r["reflector"] == reflector]
        values = np.array([r["grasse_minus_apollo_raw_mean_m"] for r in group])
        gaps = np.array([r["grasse_minus_apollo_hours"] for r in group])
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
    write_csv(args.output / "paired_station_raw_summary.csv", paired)
    provenance = dict(
        source_csv_sha256=hashlib.sha256(args.csv.read_bytes()).hexdigest(),
        source_count=len(points),
        source_sha256={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__), Path(__file__).with_name("plot_llr_residuals.py"))
        },
        matched_pairs=len(pairs),
        adjacent_date_pairs=sum(r["apollo_night"] != r["grasse_night"] for r in pairs),
        method="Raw monthly NP median/N/P10/P90; unweighted nightly means; nearest same-reflector Grasse night within12h of each APOLLO one. No fitting, centering, temporal alignment or residual adjustment.",
        caution="ALL pools reflector composition; sampling changes can move monthly medians. Paired epochs differ and are not simultaneous station biases. Empty months remain empty.",
        observation_solves=0,
    )
    (args.output / "temporal_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
