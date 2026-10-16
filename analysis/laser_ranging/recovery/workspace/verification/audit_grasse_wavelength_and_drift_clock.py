"""Evaluate source wavelength sensitivity and bound station drift-clock approximation, without range solves."""

from pathlib import Path
import ast
import hashlib
import json
import math
import sys
import numpy as np

source = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-data/analysis/laser_ranging")
sys.path.insert(0, str(source))
import llr_residuals as model
import llr_long_arc as arc
from tudatpy.astro import element_conversion, time_representation as time

research = Path("/home/dominic/Tudat/laser-ranging-work/research")
reference = research.parent / "verification/llr_optical_actual_input_audit.py"
function = next(
    node
    for node in ast.parse(reference.read_text()).body
    if isinstance(node, ast.FunctionDef) and node.name == "optical"
)
namespace = {"np": np, "math": math}
exec(compile(ast.Module(body=[function], type_ignores=[]), str(reference), "exec"), namespace)
optical = namespace["optical"]
manifest = json.loads(arc.MANIFEST.read_text())
groups = arc.read_groups(research / "longarc-20170304", manifest)
data = json.loads((source / "llr_dataset.json").read_text())
model.load_ephemeris_kernels(
    research.parent / "data/llr/kernels", data, data["analysis_frame_alias"]
)
short = model.read_points(
    research.parent / "data/llr" / data["normal_points"]["file"], data["normal_points"]["sha256"]
)
bodies, _, _ = model.create_environment(short, station_profile="inpop17a-gir")
converter = bodies.get("Earth").rotation_model.angles_calculator.time_scale_converter
optical_rows = []
clock_rows = []
for points, weather, profile, ledger in groups:
    config = model.STATION_PROFILES[profile]
    nominal = np.asarray(config["position_m"])
    velocity = np.asarray(config["velocity_m_per_julian_year"]) / config["julian_year_s"]
    for point in points:
        tdb = converter.convert_time_object(time.utc_scale, time.tdb_scale, point["utc"], nominal)
        tt_site = converter.convert_time_object(time.tdb_scale, time.tt_scale, tdb, nominal)
        tt_geocentre = converter.convert_time_object(time.tdb_scale, time.tt_scale, tdb)
        difference = float(tt_geocentre - tt_site)
        clock_rows.append(
            dict(
                station=ledger["station"],
                utc=point["utc_label"],
                dt_geocentre_minus_station_seconds=difference,
                velocity_norm_m_per_second=float(np.linalg.norm(velocity)),
                position_difference_bound_m=abs(difference) * float(np.linalg.norm(velocity)),
            )
        )
        if ledger["station"] != "GRASSE":
            continue
        times = sorted(weather)
        temperature, pressure, humidity = [
            float(np.interp(float(point["utc"]), times, [weather[t][i] for t in times]))
            for i in range(3)
        ]
        xyz = nominal + velocity * float(tt_site)
        height, latitude, longitude = element_conversion.convert_cartesian_to_geodetic_coordinates(
            xyz, 6378137.0, 1 / 298.257223563, 1e-8
        ).ravel()
        row = dict(
            utc=point["utc_label"],
            source=point["source_file"],
            wavelength_nm=point["wavelength_m"] * 1e9,
            pressure_hpa=pressure / 100,
            temperature_k=temperature,
            humidity=humidity,
        )
        for elevation in (90.0, 20.0, 3.0):
            arguments = [
                latitude,
                height,
                temperature,
                pressure,
                humidity,
                1064e-9,
                math.radians(elevation),
            ]
            first = sum(optical(*arguments))
            arguments[5] = 1064.2e-9
            second = sum(optical(*arguments))
            row[f"1064p2_minus_1064p0_{elevation:g}deg_delay_m"] = float(second - first)
        optical_rows.append(row)
output = dict(
    scope="No native observation calculation or model adjustment. All306 source points for drift clock bound; all237 Grasse weather states for independent IERS formula wavelength sensitivity. Elevations90/20/3deg are declared scenarios, not measured NP geometry.",
    model_source_sha256=model.sha256(Path(model.__file__)),
    arc_manifest_sha256=model.sha256(arc.MANIFEST),
    independent_formula_source_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
    wavelengths_nm=sorted({row["wavelength_nm"] for row in optical_rows}),
    grasse_count=len(optical_rows),
    optical_summary={},
    clock_summary={},
)
for elevation in (90.0, 20.0, 3.0):
    key = f"1064p2_minus_1064p0_{elevation:g}deg_delay_m"
    values = [row[key] for row in optical_rows]
    output["optical_summary"][f"{elevation:g}deg"] = dict(
        minimum_m=min(values), maximum_m=max(values), maximum_absolute_m=max(map(abs, values))
    )
for station in sorted({row["station"] for row in clock_rows}):
    rows = [row for row in clock_rows if row["station"] == station]
    output["clock_summary"][station] = dict(
        count=len(rows),
        maximum_abs_clock_difference_seconds=max(
            abs(row["dt_geocentre_minus_station_seconds"]) for row in rows
        ),
        maximum_position_difference_bound_m=max(row["position_difference_bound_m"] for row in rows),
        maximum_velocity_m_per_second=max(row["velocity_norm_m_per_second"] for row in rows),
    )
output["first_grasse_weather"] = optical_rows[0]
output["probe_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(research / "llr-runtime-audit/wavelength_and_drift_clock.json").write_text(
    json.dumps(output, indent=2) + "\n"
)
print(json.dumps(output, indent=2))
