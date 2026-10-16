"""Save compact portable numerical evidence from the successful frozen SLR runs."""

import csv
import hashlib
import json
import math
from pathlib import Path

work = Path("/home/dominic/Tudat/laser-ranging-work")
source_root = work / "results"
target = Path(
    "/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-data/analysis/laser_ranging/results/stable_cubic_two_weeks"
)
target.mkdir(exist_ok=False)
variants = {
    "week1_no_loading": "slr-after-cubic-current-noload-attempt1",
    "week1_historical_radial": "slr-after-cubic-historical-radial-attempt1",
    "week1_geodetic": "slr-after-cubic-current-geodetic-attempt1",
    "week2_no_loading": "slr-week2-noload-attempt1",
    "week2_geodetic": "slr-week2-geodetic-attempt1",
}
raw = {}
summaries = {}
statistics = {}
input_files = {}
output_hashes = {}
for label, directory in variants.items():
    path = source_root / directory
    with (path / "residuals.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    summary = json.loads((path / "summary.json").read_text())
    raw[label] = rows
    summaries[label] = summary
    statistics[label] = {"all": summary["all"], "stations": summary["stations"]}
    for record in summary["inputs"]:
        portable = dict(record, path=Path(record["path"]).name)
        if portable["path"] in input_files:
            assert portable == input_files[portable["path"]]
        input_files[portable["path"]] = portable
    with (target / f"{label}.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            portable = dict(row, source=Path(row["source"]).name)
            writer.writerow(portable)
    output_hashes[label] = {
        name: hashlib.sha256((path / name).read_bytes()).hexdigest()
        for name in ("residuals.csv", "summary.json")
    }


def keys(rows):
    """Identify retained observations independently of their numerical residuals."""
    return [(row["station"], row["session"], row["normal_point"]) for row in rows]


with (source_root / "slr-precise-cubic30-time-attempt1/residuals.csv").open(newline="") as stream:
    old_kernel = list(csv.DictReader(stream))
with (source_root / "slr-after-cubic-first-her-radial-attempt1/residuals.csv").open(
    newline=""
) as stream:
    first_point = list(csv.DictReader(stream))
with (source_root / "slr-after-cubic-historical-noload-attempt1/residuals.csv").open(
    newline=""
) as stream:
    historical_noload = list(csv.DictReader(stream))
assert historical_noload == raw["week1_no_loading"]
assert first_point == raw["week1_historical_radial"][:1]
assert keys(old_kernel) == keys(raw["week1_no_loading"])
for label in ("week1_historical_radial", "week1_geodetic"):
    assert keys(raw[label]) == keys(raw["week1_no_loading"])
    assert summaries[label]["excluded"] == summaries["week1_no_loading"]["excluded"]
assert keys(raw["week2_no_loading"]) == keys(raw["week2_geodetic"])
assert summaries["week2_no_loading"]["excluded"] == summaries["week2_geodetic"]["excluded"]
changes = []
for week, end in ((1, "2024-04-06"), (2, "2024-04-13")):
    baseline = raw[f"week{week}_no_loading"]
    geodetic = raw[f"week{week}_geodetic"]
    for index, (left, right) in enumerate(zip(baseline, geodetic)):
        original = float(left["computed_full_two_way_range_m"])
        loaded = float(right["computed_full_two_way_range_m"])
        row = {
            "week_ending": end,
            **{key: left[key] for key in ("station", "session", "normal_point")},
            "no_loading_computed_full_two_way_m": left["computed_full_two_way_range_m"],
            "geodetic_computed_full_two_way_m": right["computed_full_two_way_range_m"],
            "geodetic_minus_no_loading_half_range_m": 0.5 * (loaded - original),
            "historical_radial_computed_full_two_way_m": "",
            "geodetic_minus_radial_half_range_m": "",
            "stable_minus_old_cubic_no_loading_half_range_m": "",
        }
        if week == 1:
            radial = float(raw["week1_historical_radial"][index]["computed_full_two_way_range_m"])
            row["historical_radial_computed_full_two_way_m"] = raw["week1_historical_radial"][
                index
            ]["computed_full_two_way_range_m"]
            row["geodetic_minus_radial_half_range_m"] = 0.5 * (loaded - radial)
            row["stable_minus_old_cubic_no_loading_half_range_m"] = 0.5 * (
                original - float(old_kernel[index]["computed_full_two_way_range_m"])
            )
        changes.append(row)
with (target / "computed_half_range_changes.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(changes[0]))
    writer.writeheader()
    writer.writerows(changes)

metrics = {}
for key in (
    "geodetic_minus_no_loading_half_range_m",
    "geodetic_minus_radial_half_range_m",
    "stable_minus_old_cubic_no_loading_half_range_m",
):
    metrics[key] = {}
    for end in ("2024-04-06", "2024-04-13"):
        values = [
            float(row[key]) for row in changes if row["week_ending"] == end and row[key] != ""
        ]
        if values:
            metrics[key][end] = {
                "count": len(values),
                "maximum_absolute_m": max(map(abs, values)),
                "rms_m": math.sqrt(sum(value * value for value in values) / len(values)),
            }
current = summaries["week1_no_loading"]
radial = summaries["week1_historical_radial"]
geodetic = summaries["week1_geodetic"]
software = {
    "kernel_sha256": current["software"]["kernel"]["sha256"],
    "kernel_native_source_commit": "b1aa772503070347e8ac0b263c8de9278d390c17",
    "current_analysis_source_commit": "9a912e85ee38d07e235f93bffcf8578a395b91d5",
    "current_analysis_script_sha256": current["software"]["script"]["sha256"],
    "historical_analysis_source_commit": "74c207ff30d355738a8f73b65334b1bc955bd826",
    "historical_analysis_script_sha256": radial["software"]["script"]["sha256"],
    "historical_radial_helper_sha256": radial["software"]["ocean_loading_helper"]["sha256"],
    "current_geodetic_helper_sha256": geodetic["software"]["ocean_loading_helper"]["sha256"],
    "numerical_threads": 1,
    "native_build_provenance_sha256": hashlib.sha256(
        (work / "logs/stable-cubic-kernel-provenance.json").read_bytes()
    ).hexdigest(),
}
assert (
    software["kernel_sha256"] == "10f8307ba799847f11d9614d2c3f3da3db23962c46a739d73f536f9dc31ae5d5"
)
weeks = {}
for week, end, raw_count in ((1, "2024-04-06", 29), (2, "2024-04-13", 134)):
    summary = summaries[f"week{week}_no_loading"]
    weeks[end] = {
        "raw_point_count": raw_count,
        "accepted_point_count": summary["all"]["count"],
        "excluded": summary["excluded"],
        "same_accepted_keys_and_exclusions_with_loading": True,
        "ephemeris_routes": summary["ephemeris_routes"],
        "sp3_holdout": summary["sp3_holdout"],
    }
loading = []
for row in geodetic["ocean_loading_metadata"]:
    portable = dict(row, file=Path(row["file"]).name)
    loading.append(portable)
record = {
    "experiment": "Two consecutive predeclared LAGEOS1 combined-postfit weeks, fixed SLRF2020 positions/velocities and dated CoM/DHF; no propagation, estimation, fitted bias, demeaning or residual clipping.",
    "default": "No ocean loading; geodetic HARPOS is an explicit alternate realization, not exact ASI BLQ/HARDISP.",
    "portable_files": "Raw CSV numerical fields unchanged; absolute source paths replaced by basenames. Historical/current duplicate no-loading CSV is stored once.",
    "provider_input_manifest": "../../dataset.json",
    "software": software,
    "statistics_one_way_observed_minus_computed_metres": statistics,
    "weeks": weeks,
    "computed_half_range_change_metrics": metrics,
    "loading_metadata": loading,
    "geodetic_basis_provenance": "../../ocean_loading_basis.json",
    "historical_radial_failure": "../got410c_harpos_radial_failure/failure.json",
    "same_polynomial_diagnosis": "../got410c_harpos_radial_failure/cubic_same_polynomial_proof.json",
    "first_historical_failing_point_rerun": {
        "status": "PASS with original50iteration/1e-12fractional settings and exact original script/helper",
        "station": "7840",
        "session": 29,
        "normal_point": 0,
        "computed_full_two_way_range_m": float(first_point[0]["computed_full_two_way_range_m"]),
        "observed_minus_computed_m": float(first_point[0]["observed_minus_computed_m"]),
        "identical_to_first_point_in_full_historical_radial_run": True,
    },
    "date_parameterization_control": "Same-kernel historical and current no-loading CSV rows match every field; original23 keys/six exclusions retained. Cross-kernel changes are measured, not required to vanish.",
    "input_files": list(input_files.values()),
    "original_output_sha256": output_hashes,
    "convergence_unchanged": {"maximum_iterations": 50, "fractional_tolerance": 1e-12},
    "precision_limits": "Time callback preserves precise epochs; rotation and BCRS state components retain double arithmetic. Source SP3 print resolution is1mm; no end-to-end submillimetre accuracy claim.",
}
(target / "summary.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
print(
    json.dumps(
        {"saved": str(target), "statistics": statistics, "half_range_metrics": metrics}, indent=2
    )
)
