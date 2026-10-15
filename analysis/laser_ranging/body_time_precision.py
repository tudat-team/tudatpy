"""Compare computed half-ranges from two explicit fixed-orbit SLR epoch routes."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def read_result(directory):
    """Read unchanged numerical outputs and their recorded experiment provenance."""
    with (directory / "residuals.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))

    return rows, json.loads((directory / "summary.json").read_text())


def main():
    """Require the same physical experiment before comparing its computed ranges."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    baseline, old = read_result(args.baseline)
    candidate, new = read_result(args.candidate)
    keys = lambda rows: [(r["station"], r["session"], r["normal_point"]) for r in rows]
    hashes = lambda summary: {Path(r["path"]).name: r["sha256"] for r in summary["inputs"]}

    # Route comparisons preserve data, corrections, point selection and compiled code.
    if keys(baseline) != keys(candidate) or old["excluded"] != new["excluded"]:
        raise ValueError("The compared runs have different accepted points or exclusions")
    if (
        hashes(old) != hashes(new)
        or old["software"]["kernel"]["sha256"] != new["software"]["kernel"]["sha256"]
        or old["software"]["script"]["sha256"] != new["software"]["script"]["sha256"]
    ):
        raise ValueError("The compared runs require identical inputs, analysis script and kernel")
    for name in (
        "orbit",
        "station_solution",
        "sp3_interpolation_points",
        "station_coordinates",
        "ocean_loading",
    ):
        if old[name] != new[name]:
            raise ValueError(f"The compared runs differ in {name}")

    rows = []

    for left, right in zip(baseline, candidate):
        rows.append(
            {
                "station": left["station"],
                "session": left["session"],
                "normal_point": left["normal_point"],
                "computed_half_range_candidate_minus_baseline_m": 0.5
                * (
                    float(right["computed_full_two_way_range_m"])
                    - float(left["computed_full_two_way_range_m"])
                ),
            }
        )

    values = [row["computed_half_range_candidate_minus_baseline_m"] for row in rows]
    summary = {
        "baseline_routes": old["ephemeris_routes"],
        "candidate_routes": new["ephemeris_routes"],
        "point_count": len(rows),
        "same_accepted_point_keys_and_exclusions": True,
        "maximum_absolute_half_range_change_m": max(abs(value) for value in values),
        "rms_half_range_change_m": math.sqrt(sum(value * value for value in values) / len(values)),
        "kernel_sha256": old["software"]["kernel"]["sha256"],
        "source_output_sha256": {
            f"{label}/{name}": hashlib.sha256((directory / name).read_bytes()).hexdigest()
            for label, directory in (("baseline", args.baseline), ("candidate", args.candidate))
            for name in ("residuals.csv", "summary.json")
        },
        "interpretation": "Computed half-range route difference; no residual fitting or selection by improvement.",
    }
    args.output.mkdir(parents=True, exist_ok=False)

    with (args.output / "half_range_changes.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
