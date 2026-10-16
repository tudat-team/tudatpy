"""Package the predeclared six-case matrix without fitting or modifying raw rows."""

from pathlib import Path
import csv
import hashlib
import json
import shutil
import statistics
import subprocess

work = Path("/home/dominic/Tudat/laser-ranging-work")
source = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-data")
output = source / "analysis/laser_ranging/results/coherent_station_ephemeris_controls"
output.mkdir(exist_ok=False)
runtime = json.loads((work / "logs/llr-six-controls-provenance.json").read_text())
assert len(runtime["runs"]) == 6
report = {
    "native_source_revision": runtime["native_source_revision"],
    "kernel_binary_sha256": runtime["kernel_binary_sha256"],
    "analysis_source_revision": runtime["analysis_source_revision"],
    "numerical_threads": 1,
    "fractional_light_time_tolerance": 1e-13,
    "maximum_iterations": 50,
    "no_fitting_or_point_selection": True,
    "sources": {},
    "cases": {},
}
for name in [
    "llr_residuals.py",
    "matera_case.py",
    "llr_dataset.json",
    "matera_dataset.json",
    "llr_station_profiles.json",
    "hofmann2017_station_profiles.json",
    "de421_dataset.json",
    "icrf_spice_alias.tf",
]:
    path = source / "analysis/laser_ranging" / name
    report["sources"][name] = {
        "path": str(path.relative_to(source)),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "blob": subprocess.check_output(
            [
                "git",
                "rev-parse",
                f"{runtime['analysis_source_revision']}:analysis/laser_ranging/{name}",
            ],
            cwd=source,
            text=True,
        ).strip(),
    }
rows_by_case = {}
for record in runtime["runs"]:
    label = record["label"]
    original = Path(record["output"])
    summary = json.loads((original / "summary.json").read_text())
    shutil.copyfile(original / "normal_points.csv", output / f"{label}.csv")
    with (original / "normal_points.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    rows_by_case[label] = rows
    report["cases"][label] = {
        key: summary[key]
        for key in [
            "case",
            "ephemeris",
            "count",
            "mean_m",
            "rms_about_zero_m",
            "scatter_about_mean_m",
            "weighted_mean_m",
            "normalized_rms",
            "chi_squared_no_fitted_parameters",
            "normal_point_sha256",
            "dataset_manifest_sha256",
            "kernel_manifest_sha256",
            "ephemeris_constants_source",
            "gravitational_parameters_m3_s2",
            "reflector_coordinate_metadata",
            "kernel_sha256",
            "inertial_axes",
            "station_profile",
            "station_xyz_m",
            "station_profile_manifest_sha256",
            "hofmann_station_profile_manifest_sha256",
            "station_profile_metadata",
            "matera_adapter_sha256",
            "spice_knot_step_s",
            "ocean_loading",
        ]
    }
    report["cases"][label].update(
        {
            "mean_by_reflector_m": {
                code: statistics.mean(
                    float(r["one_way_residual_m"]) for r in rows if r["reflector"] == code
                )
                for code in sorted({r["reflector"] for r in rows})
            },
            "csv_sha256": hashlib.sha256((original / "normal_points.csv").read_bytes()).hexdigest(),
            "original_summary_sha256": hashlib.sha256(
                (original / "summary.json").read_bytes()
            ).hexdigest(),
            "event_diagnostics_sha256": hashlib.sha256(
                (original / "event_diagnostics.json").read_bytes()
            ).hexdigest(),
            "ephemeris_controls_sha256": hashlib.sha256(
                (original / "ephemeris_controls.json").read_bytes()
            ).hexdigest(),
        }
    )
    if "earth_orientation_resource_sha256" not in report:
        report["earth_orientation_resource_sha256"] = summary["earth_orientation_resource_sha256"]
    else:
        assert (
            report["earth_orientation_resource_sha256"]
            == summary["earth_orientation_resource_sha256"]
        )
changes = []
for baseline, candidate in [
    ("de440-inpop17a", "de440-hofmann-apollo"),
    ("de440-matera", "de440-hofmann-matera"),
    ("de440-inpop17a", "de421-inpop17a"),
    ("de440-matera", "de421-matera"),
]:
    left = rows_by_case[baseline]
    right = rows_by_case[candidate]
    assert [(r["index"], r["utc_launch"], r["reflector"]) for r in left] == [
        (r["index"], r["utc_launch"], r["reflector"]) for r in right
    ]
    for a, b in zip(left, right):
        changes.append(
            {
                "baseline": baseline,
                "candidate": candidate,
                "index": a["index"],
                "utc_launch": a["utc_launch"],
                "reflector": a["reflector"],
                "computed_half_range_candidate_minus_baseline_m": float(a["one_way_residual_m"])
                - float(b["one_way_residual_m"]),
            }
        )
with (output / "computed_half_range_changes.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=changes[0])
    writer.writeheader()
    writer.writerows(changes)
(output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
print(output)
