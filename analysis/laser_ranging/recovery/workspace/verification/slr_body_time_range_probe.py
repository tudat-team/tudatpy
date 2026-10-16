"""Measure existing cubic Time body tables against the frozen direct-SPICE SLR baseline."""

import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from tudatpy import kernel
from tudatpy.astro import time_representation as time
from tudatpy.dynamics import environment_setup
from tudatpy.dynamics.environment_setup import ephemeris
from tudatpy.interface import spice
from tudatpy.math import interpolators

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--step", type=int, choices=(30, 60), required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
expected = Path("/home/dominic/Tudat/laser-ranging-build/src/tudatpy/kernel.so")
assert Path(kernel.__file__).resolve() == expected
assert (
    hashlib.sha256(expected.read_bytes()).hexdigest()
    == "1864ed04530210fa056e78128bdbbfe1ad536cae6a9d3ade4f23295fbeda1744"
)
source = Path(
    "/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-data/analysis/laser_ranging/fixed_orbit_slr.py"
)
spec = importlib.util.spec_from_file_location("fixed_orbit_slr", source)
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)
original_defaults = environment_setup.get_default_body_settings
original_create = environment_setup.create_system_of_bodies
captured = {}
start = float(analysis.calendar_epoch(2024, 3, 31)) - 86400
end = float(analysis.calendar_epoch(2024, 4, 7)) + 86400
spline = interpolators.cubic_spline_interpolation(
    boundary_interpolation=interpolators.throw_exception_at_boundary
)


def body_time_defaults(*positional, **keywords):
    """Replace only shared body ephemerides with the existing cubic Time factory."""
    settings = original_defaults(*positional, **keywords)
    for name in ("Earth", "Sun", "Moon"):
        settings.get(name).ephemeris_settings = ephemeris.interpolated_spice(
            start, end, args.step, "SSB", "J2000", spline, name
        )
    return settings


def capture_bodies(*positional, **keywords):
    """Keep the actual analysis bodies for independent state-route controls."""
    bodies = original_create(*positional, **keywords)
    captured["bodies"] = bodies
    return bodies


environment_setup.get_default_body_settings = body_time_defaults
environment_setup.create_system_of_bodies = capture_bodies
sys.argv = [
    str(source),
    "--data-root",
    "/home/dominic/Tudat/laser-ranging-work/data/slr",
    "--output",
    str(args.output),
]
analysis.main()

bodies = captured["bodies"]
converter = bodies.get("Earth").rotation_model.angles_calculator.time_scale_converter
baseline = Path(
    "/home/dominic/Tudat/laser-ranging-work/results/slr-baseline-combined-slrf2020-20261016-attempt3"
)
with (baseline / "residuals.csv").open(newline="") as stream:
    reference_rows = list(csv.DictReader(stream))
with (args.output / "residuals.csv").open(newline="") as stream:
    table_rows = list(csv.DictReader(stream))
keys = lambda rows: [(r["session"], r["normal_point"], r["station"]) for r in rows]
assert keys(reference_rows) == keys(table_rows)
old_summary = json.loads((baseline / "summary.json").read_text())
new_summary = json.loads((args.output / "summary.json").read_text())
assert old_summary["excluded"] == new_summary["excluded"]
range_changes = [
    0.5
    * (float(new["computed_full_two_way_range_m"]) - float(old["computed_full_two_way_range_m"]))
    for old, new in zip(reference_rows, table_rows)
]
controls = {}
for name in ("Earth", "Sun", "Moon"):
    interpolation_errors = []
    precise_epoch_displacements = []
    for row in reference_rows:
        utc = time.Time(int(row["utc_hours_since_j2000"]) * 3600) + float(
            row["utc_seconds_into_hour"]
        )
        nominal = np.asarray(
            bodies.get("Earth")
            .get_ground_station(row["domes"])
            .station_state.cartesian_position_at_reference_epoch
        ).reshape(3)
        transmit = converter.convert_time_object(time.utc_scale, time.tdb_scale, utc, nominal)
        for offset in (0.0, float(row["tof_seconds"]) / 2, float(row["tof_seconds"])):
            precise = transmit + offset
            rounded = float(precise)
            exact_state = np.asarray(
                bodies.get(name).state_in_base_frame_from_ephemeris(precise)
            ).reshape(6)
            rounded_state = np.asarray(
                bodies.get(name).state_in_base_frame_from_ephemeris(time.Time(rounded))
            ).reshape(6)
            direct_state = np.asarray(
                spice.get_body_cartesian_state_at_epoch(name, "SSB", "J2000", "NONE", rounded)
            ).reshape(6)
            interpolation_errors.append(float(np.linalg.norm(rounded_state[:3] - direct_state[:3])))
            precise_epoch_displacements.append(
                float(np.linalg.norm(exact_state[:3] - rounded_state[:3]))
            )
    controls[name] = {
        "count": len(interpolation_errors),
        "maximum_table_vs_direct_at_equal_double_epoch_m": max(interpolation_errors),
        "rms_table_vs_direct_at_equal_double_epoch_m": float(
            np.sqrt(np.mean(np.square(interpolation_errors)))
        ),
        "maximum_table_time_vs_table_rounded_epoch_position_m": max(precise_epoch_displacements),
    }

result = {
    "kind": "controlled body-only cubic Time route; identical custom double satellite callback",
    "kernel_sha256": hashlib.sha256(expected.read_bytes()).hexdigest(),
    "script_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "wrapper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "step_seconds": args.step,
    "table_start_tdb": start,
    "table_end_tdb": end,
    "same_accepted_point_keys_and_exclusions": True,
    "point_count": len(range_changes),
    "half_range_table_minus_baseline_m": range_changes,
    "maximum_absolute_half_range_change_m": max(abs(value) for value in range_changes),
    "rms_half_range_change_m": float(np.sqrt(np.mean(np.square(range_changes)))),
    "body_controls": controls,
}
(args.output / "precision-comparison.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
