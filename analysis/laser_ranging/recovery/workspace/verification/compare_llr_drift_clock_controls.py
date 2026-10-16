"""Compare every raw field/geometry and physical provenance for the authorized306-point pair."""

from pathlib import Path
import csv, hashlib, json, ast, subprocess

WORK = Path("/home/dominic/Tudat/laser-ranging-work")
SOURCE = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-signature")
base = WORK / "results/llr-drift-clock-baseline-306"
candidate = WORK / "results/llr-drift-clock-candidate-306"
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
b = json.loads((base / "summary.json").read_text())
c = json.loads((candidate / "summary.json").read_text())
with (base / "residuals.csv").open() as stream:
    reader = csv.DictReader(stream)
    fields = reader.fieldnames
    rows_b = list(reader)
rows_c = list(csv.DictReader((candidate / "residuals.csv").open()))
assert len(rows_b) == len(rows_c) == b["count"] == c["count"] == 306
assert (base / "residuals.csv").read_bytes() == (candidate / "residuals.csv").read_bytes()
assert rows_b == rows_c
assert (base / "nightly_raw_statistics.csv").read_bytes() == (
    candidate / "nightly_raw_statistics.csv"
).read_bytes()
keys = [
    "manifest_sha256",
    "ephemeris_manifest_sha256",
    "kernel_binary_sha256",
    "kernel_sha256",
    "earth_orientation_resource_sha256",
    "station_profiles",
    "physics",
    "fitting",
    "geometry_precision",
    "uncertainty",
    "station_statistics",
    "station_reflector_statistics",
]
for key in keys:
    assert b[key] == c[key], key
assert len(b["groups"]) == len(c["groups"])
for left, right in zip(b["groups"], c["groups"]):
    assert {k: v for k, v in left.items() if k not in ["setup_seconds", "solver_seconds"]} == {
        k: v for k, v in right.items() if k not in ["setup_seconds", "solver_seconds"]
    }
assert b["source_sha256"]["llr_long_arc.py"] == c["source_sha256"]["llr_long_arc.py"]
assert (
    b["source_sha256"]["llr_station_profiles.json"]
    == c["source_sha256"]["llr_station_profiles.json"]
)
old = subprocess.check_output(
    ["git", "show", "3dc1574d9:analysis/laser_ranging/llr_residuals.py"], cwd=SOURCE, text=True
)
new = (SOURCE / "analysis/laser_ranging/llr_residuals.py").read_text()
old_tree = ast.parse(old)
new_tree = ast.parse(new)
changed = [
    getattr(x, "name", type(x).__name__)
    for x, y in zip(old_tree.body, new_tree.body)
    if ast.dump(x) != ast.dump(y)
]
assert len(old_tree.body) == len(new_tree.body) and changed == ["station_kinematic_displacement"]
clock = json.loads(
    (WORK / "research/llr-runtime-audit/wavelength_and_drift_clock.json").read_text()
)
result = {
    "decision": "retain: all raw results/geometry/source identities preserved, clear solver-time benefit",
    "baseline_revision": "3dc1574d9cf29cf247efb040585e2e561cb4bf21",
    "kernel_sha256": b["kernel_binary_sha256"],
    "observations": 306,
    "raw_csv_columns": len(fields),
    "raw_csv_fields_compared": 306 * len(fields),
    "raw_csv_bytes_identical": True,
    "raw_csv_sha256": sha(base / "residuals.csv"),
    "nightly_statistics_bytes_identical": True,
    "geometry_fields_compared": [
        key
        for key in fields
        if "elevation" in key or "direction" in key or "los_pa" in key or "phase" in key
    ],
    "source_identity_summary_keys_equal": keys,
    "all_group_source_weather_metadata_equal": True,
    "only_changed_function": changed[0],
    "baseline_model_source_sha256": b["source_sha256"]["llr_residuals.py"],
    "candidate_model_source_sha256": c["source_sha256"]["llr_residuals.py"],
    "timings": {
        key: {
            "baseline_s": b[key],
            "candidate_s": c[key],
            "reduction_fraction": 1 - c[key] / b[key],
        }
        for key in ["setup_seconds", "solver_seconds", "total_seconds"]
    },
    "external_wall": {
        label: json.loads(
            (WORK / "logs" / f"llr-drift-clock-{label}-306-provenance.json").read_text()
        )["external_wall_seconds"]
        for label in ["baseline", "candidate"]
    },
    "existing_independent_drift_bound": {
        "source": str(WORK / "research/llr-runtime-audit/wavelength_and_drift_clock.json"),
        "sha256": sha(WORK / "research/llr-runtime-audit/wavelength_and_drift_clock.json"),
        "clock_summary": clock["clock_summary"],
    },
    "scope": "Only linear cm/year station drift uses geocentric TT. Observable clocks, Earth simultaneity, generic converter and all native code unchanged. One sequential baseline/candidate pair on same machine and pinned kernel, all numerical thread limits1; no repeat needed for observed35percent solver reduction. Timing is a measured local result, not a general performance guarantee.",
}
(WORK / "logs/llr-drift-clock-comparison.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
