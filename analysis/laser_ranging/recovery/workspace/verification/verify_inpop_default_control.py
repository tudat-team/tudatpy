"""Verify one default DE440 night exactly against the committed optimized core."""

from pathlib import Path
import hashlib, json, subprocess, sys, types
import numpy as np

source = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-llr/analysis/laser_ranging")
sys.path.insert(0, str(source))
import llr_residuals as current
import llr_long_arc as arc

root = Path("/home/dominic/Tudat/laser-ranging-work")
text = subprocess.check_output(
    ["git", "show", "39f3c7b8c:analysis/laser_ranging/llr_residuals.py"], cwd=source, text=True
)
reference = types.ModuleType("committed_reference")
reference.__file__ = str(source / "llr_residuals.py")
exec(compile(text, "39f3c7b8c:analysis/laser_ranging/llr_residuals.py", "exec"), reference.__dict__)
manifest = json.loads((source / "llr_dataset.json").read_text())
groups = arc.read_groups(root / "research/longarc-20170304", json.loads(arc.MANIFEST.read_text()))
points, weather, profile, ledger = next(
    g for g in groups if g[3]["station"] == "APOLLO" and g[3]["group"] == "2017-03-09"
)
outputs = {}
for label, model in [("reference", reference), ("modified_default", current)]:
    model.load_ephemeris_kernels(
        root / "data/llr/kernels", manifest, manifest["analysis_frame_alias"]
    )
    bodies, position, _ = model.create_environment(
        points, station_profile=profile, meteorological_history=weather
    )
    converter = bodies.get("Earth").rotation_model.angles_calculator.time_scale_converter
    values = []
    for point in points:
        epoch = converter.convert_time_object(
            model.time.utc_scale, model.time.tdb_scale, point["utc"], position
        )
        computed, diagnostic = model.computed_range(
            point, epoch, bodies, model.complete_light_time_corrections(point)
        )
        values.append(
            dict(utc=point["utc_label"], reflector=point["reflector"], computed_range_m=computed)
        )
    outputs[label] = values
result = dict(
    reference_revision="39f3c7b8c",
    reference_source_sha256=hashlib.sha256(text.encode()).hexdigest(),
    modified_source_sha256=current.sha256(source / "llr_residuals.py"),
    kernel_sha256=current.sha256(Path(current.kernel.__file__)),
    count=len(points),
    all_rows_exact=outputs["reference"] == outputs["modified_default"],
    maximum_range_difference_m=max(
        abs(a["computed_range_m"] - b["computed_range_m"]) for a, b in zip(*outputs.values())
    ),
    rows=outputs,
)
assert result["all_rows_exact"]
(root / "results/llr-inpop17a-default-control.json").write_text(json.dumps(result, indent=2) + "\n")
print({k: v for k, v in result.items() if k != "rows"})
