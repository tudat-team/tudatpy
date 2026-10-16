"""Measure two existing normal points and native function costs without changing model code."""

from pathlib import Path
from time import perf_counter
import cProfile
import hashlib
import json
import pstats
import subprocess
import sys

source = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-data/analysis/laser_ranging")
sys.path.insert(0, str(source))
import llr_residuals as model
from tudatpy import kernel
from tudatpy.astro import time_representation as time

output = Path("/home/dominic/Tudat/laser-ranging-work/research/llr-runtime-audit")
output.mkdir(exist_ok=True)
manifest = json.loads((source / "llr_dataset.json").read_text())
started = perf_counter()
points = model.read_points(
    Path("/home/dominic/Tudat/laser-ranging-work/data/llr") / manifest["normal_points"]["file"],
    manifest["normal_points"]["sha256"],
)
model.load_ephemeris_kernels(
    Path("/home/dominic/Tudat/laser-ranging-work/data/llr/kernels"),
    manifest,
    manifest["analysis_frame_alias"],
)
load_seconds = perf_counter() - started
started = perf_counter()
bodies, nominal, _ = model.create_environment(points, station_profile="inpop17a-gir")
setup_seconds = perf_counter() - started
converter = bodies.get("Earth").rotation_model.angles_calculator.time_scale_converter
profile = cProfile.Profile()
ranges = []
for point in points[:2]:
    epoch = converter.convert_time_object(time.utc_scale, time.tdb_scale, point["utc"], nominal)
    correction = model.complete_light_time_corrections(point)
    started = perf_counter()
    profile.enable()
    computed, diagnostic = model.computed_range(point, epoch, bodies, correction)
    profile.disable()
    ranges.append(
        dict(index=point["index"], complete_range_m=computed, seconds=perf_counter() - started)
    )
profile.dump_stats(str(output / "two_points.prof"))
stats = pstats.Stats(profile)
entries = [
    dict(
        file=Path(key[0]).name,
        line=key[1],
        name=key[2],
        primitive_calls=value[0],
        calls=value[1],
        self_seconds=value[2],
        cumulative_seconds=value[3],
    )
    for key, value in stats.stats.items()
]
epoch = converter.convert_time_object(time.utc_scale, time.tdb_scale, points[0]["utc"], nominal)
advancing = [epoch + i * 0.013 for i in range(100)]
repeated = [epoch for _ in range(100)]
earth = bodies.get("Earth")
moon = bodies.get("Moon")
reflector = model.REFLECTORS[points[0]["reflector"]][0]
functions = {
    "earth_ephemeris": earth.state_in_base_frame_from_ephemeris,
    "earth_rotation": earth.rotation_model.body_fixed_to_inertial_rotation,
    "earth_rotation_derivative": earth.rotation_model.time_derivative_body_fixed_to_inertial_rotation,
    "moon_rotation": moon.rotation_model.body_fixed_to_inertial_rotation,
    "tdb_to_tt": lambda t: converter.convert_time_object(time.tdb_scale, time.tt_scale, t),
    "tdb_to_utc": lambda t: converter.convert_time_object(time.tdb_scale, time.utc_scale, t),
    "earth_intrinsic_position": lambda t: earth.get_ground_station(
        "APOLLO_intrinsic"
    ).station_state.get_cartesian_position(t, "Earth"),
    "earth_final_position": lambda t: earth.get_ground_station(
        "APOLLO"
    ).station_state.get_cartesian_position(t, "Earth"),
    "moon_intrinsic_position": lambda t: moon.get_ground_station(
        reflector + "_intrinsic"
    ).station_state.get_cartesian_position(t, "Moon"),
    "moon_final_position": lambda t: moon.get_ground_station(
        reflector
    ).station_state.get_cartesian_position(t, "Moon"),
}
measurements = {}
for name, function in functions.items():
    measurements[name] = {}
    for label, epochs in [("advancing", advancing), ("same_epoch", repeated)]:
        function(epochs[0])
        started = perf_counter()
        for t in epochs:
            function(t)
        measurements[name][label + "_microseconds_per_call"] = (
            (perf_counter() - started) * 1e6 / len(epochs)
        )
result = dict(
    scope="Two frozen complete APOLLO NPs only; unchanged current model.100 native calls per microprobe; no optimization or implementation edit.",
    source_commit=subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=source, text=True
    ).strip(),
    script_sha256=hashlib.sha256((source / "llr_residuals.py").read_bytes()).hexdigest(),
    kernel_sha256=hashlib.sha256(Path(kernel.__file__).read_bytes()).hexdigest(),
    input_and_kernel_load_seconds=load_seconds,
    environment_setup_seconds=setup_seconds,
    normal_points=ranges,
    microprobes=measurements,
    profile_top_self=sorted(entries, key=lambda e: e["self_seconds"], reverse=True)[:25],
    profile_top_cumulative=sorted(entries, key=lambda e: e["cumulative_seconds"], reverse=True)[
        :25
    ],
)
(output / "measurements.json").write_text(json.dumps(result, indent=2) + "\n")
print(
    json.dumps(
        {
            key: result[key]
            for key in [
                "source_commit",
                "kernel_sha256",
                "input_and_kernel_load_seconds",
                "environment_setup_seconds",
                "normal_points",
                "microprobes",
                "profile_top_self",
            ]
        },
        indent=2,
    )
)
