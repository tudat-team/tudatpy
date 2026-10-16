"""Package four unchanged-input LLR runs and separate their numerical controls."""

import csv
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path

WORK = Path("/home/dominic/Tudat/laser-ranging-work/results")
REPO = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-llr")
DEST = REPO / "analysis/laser_ranging/results/stable_cubic_20261016"
DEST.mkdir(exist_ok=True)
CASES = {
    "apollo_survey": "apollo-20170309-survey-stable-cubic30",
    "apollo_inpop17a": "apollo-20170309-inpop17a-stable-cubic30",
    "matera": "matera-20220821-stable-cubic30",
    "apollo_inpop17a_geodetic_got410c": "apollo-20170309-inpop17a-geodetic-got410c-stable-cubic30",
}
PRIORS = {
    "apollo_survey": (
        "apollo-20170309-matera-refactor-control",
        "../apollo_20170309/survey_cubic30/normal_points.csv",
    ),
    "apollo_inpop17a": (
        "apollo-20170309-inpop17a-gir-cubic30",
        "../apollo_20170309/inpop17a_gir_cubic30/normal_points.csv",
    ),
    "matera": ("matera-20220821-de440-cubic30-committed", "../matera_20220821/normal_points.csv"),
}
FIELDS = [
    "one_way_residual_m",
    "geometric_atomic_half_range_m",
    "shapiro_increment_m",
    "atmosphere_increment_m",
    "glass_increment_m",
]


def sha(path):
    """Identify exact external and committed artifacts."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(*args):
    """Return a local Git object identifier without changing any source."""
    return subprocess.check_output(["git", *args], cwd=REPO, text=True).strip()


def read_csv(path):
    """Read a raw CSV while retaining the original field strings."""
    with path.open() as stream:
        return list(csv.DictReader(stream))


summaries = {
    case: json.loads((WORK / name / "summary.json").read_text()) for case, name in CASES.items()
}
first = summaries["apollo_survey"]
assert all(
    s["kernel_binary_sha256"] == "10f8307ba799847f11d9614d2c3f3da3db23962c46a739d73f536f9dc31ae5d5"
    for s in summaries.values()
)
shared_keys = [
    "kernel_binary_sha256",
    "script_sha256",
    "kernel_sha256",
    "earth_orientation_resource_sha256",
    "spice_knot_step_s",
]
for key in shared_keys:
    assert all(s[key] == first[key] for s in summaries.values())
report = {
    "native_source_revision": "b1aa772503070347e8ac0b263c8de9278d390c17",
    "shared": {key: first[key] for key in shared_keys},
    "runtime": {
        "package_relative_to_build": "src/tudatpy/__init__.py",
        "kernel_relative_to_build": "src/tudatpy/kernel.so",
        "interpreter_environment": "tudatpy-dev",
        "numerical_threads": 1,
    },
    "fractional_light_time_tolerance": 1e-13,
    "maximum_iterations": 50,
    "reconstruction": "git show COMMIT:PATH or git show BLOB gives each exact executed source/manifest; verify SHA256. Raw external summaries retain their original import paths and are separately identified below.",
    "sources": {},
    "cases": {},
    "cross_kernel_controls": {},
    "comparison_convention": "All differences are fresh minus prior. No-loading controls compare unchanged physical models. Historical radial-loading outputs remain untouched and are not used as pure kernel controls.",
}
for file in [
    "llr_residuals.py",
    "matera_case.py",
    "ocean_loading.py",
    "ocean_loading_basis.json",
    "llr_dataset.json",
    "matera_dataset.json",
    "llr_station_profiles.json",
]:
    path = "analysis/laser_ranging/" + file
    report["sources"][file] = {
        "commit": git("rev-parse", "HEAD"),
        "path": path,
        "blob": git("rev-parse", "HEAD:" + path),
        "sha256": sha(REPO / path),
    }
allrows = {}
for case, name in CASES.items():
    run = WORK / name
    target = DEST / (case + ".csv")
    shutil.copyfile(run / "normal_points.csv", target)
    assert target.read_bytes() == (run / "normal_points.csv").read_bytes()
    rows = read_csv(target)
    allrows[case] = rows
    s = summaries[case]
    r = {
        key: s[key]
        for key in [
            "count",
            "mean_m",
            "rms_about_zero_m",
            "scatter_about_mean_m",
            "weighted_mean_m",
            "normalized_rms",
            "chi_squared_no_fitted_parameters",
            "normal_point_sha256",
            "dataset_manifest_sha256",
            "station_profile",
            "station_xyz_m",
        ]
    }
    r["csv_sha256"] = sha(target)
    r["original_summary_sha256"] = sha(run / "summary.json")
    r["event_diagnostics_sha256"] = sha(run / "event_diagnostics.json")
    r["ephemeris_controls_sha256"] = sha(run / "ephemeris_controls.json")
    r["mean_by_reflector_m"] = {
        ref: sum(float(row["one_way_residual_m"]) for row in rows if row["reflector"] == ref)
        / sum(row["reflector"] == ref for row in rows)
        for ref in sorted({row["reflector"] for row in rows})
    }
    if s["ocean_loading"] is not None:
        r["ocean_loading"] = {
            key: value for key, value in s["ocean_loading"].items() if key != "file"
        }
    report["cases"][case] = r

changes = []
for case, (prior_name, prior_path) in PRIORS.items():
    prior = WORK / prior_name
    old = read_csv(prior / "normal_points.csv")
    new = allrows[case]
    assert len(old) == len(new)
    assert sha(DEST / prior_path) == sha(prior / "normal_points.csv")
    old_summary = json.loads((prior / "summary.json").read_text())
    for a, b in zip(old, new):
        assert all(
            a[k] == b[k]
            for k in [
                "index",
                "utc_launch",
                "reflector",
                "photons",
                "observed_rtt_s",
                "one_way_sigma_m",
            ]
        )
        changes.append(
            {
                "case": case,
                "index": b["index"],
                "reflector": b["reflector"],
                **{key: float(b[key]) - float(a[key]) for key in FIELDS},
            }
        )
    ds = [row for row in changes if row["case"] == case]
    report["cross_kernel_controls"][case] = {
        "prior_csv": prior_path,
        "prior_csv_sha256": sha(prior / "normal_points.csv"),
        "prior_kernel_sha256": old_summary["kernel_binary_sha256"],
        "prior_summary_sha256": sha(prior / "summary.json"),
        "mean_residual_change_m": sum(row["one_way_residual_m"] for row in ds) / len(ds),
        "max_abs_residual_change_m": max(abs(row["one_way_residual_m"]) for row in ds),
    }
with (DEST / "cross_kernel_differences.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(
        stream, fieldnames=["case", "index", "reflector", *FIELDS], lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(changes)
plain = allrows["apollo_inpop17a"]
loaded = allrows["apollo_inpop17a_geodetic_got410c"]
loading_differences = [
    float(b["one_way_residual_m"]) - float(a["one_way_residual_m"]) for a, b in zip(plain, loaded)
]
report["same_kernel_loading_comparison"] = {
    "mean_residual_change_m": sum(loading_differences) / len(loading_differences),
    "minimum_residual_change_m": min(loading_differences),
    "maximum_residual_change_m": max(loading_differences),
    "mean_by_reflector_m": {
        ref: sum(d for row, d in zip(plain, loading_differences) if row["reflector"] == ref)
        / sum(row["reflector"] == ref for row in plain)
        for ref in sorted({row["reflector"] for row in plain})
    },
}
(DEST / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
assert "/home/" not in (DEST / "summary.json").read_text()
print(
    "Exact copied CSVs and unchanged observation identities verified. Cases:",
    {case: len(rows) for case, rows in allrows.items()},
)
print("Portable package", DEST)
