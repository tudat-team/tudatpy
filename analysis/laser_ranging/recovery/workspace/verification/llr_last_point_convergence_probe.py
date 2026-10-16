"""Tighten only the final-point solve stopping condition; keep physics unchanged."""

import importlib.util
import json
from pathlib import Path

import numpy as np
from tudatpy import constants, kernel

source = Path(
    "/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-llr/analysis/laser_ranging/llr_residuals.py"
)
spec = importlib.util.spec_from_file_location("llr_diagnostic", source)
llr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(llr)
dataset = json.loads(source.with_name("llr_dataset.json").read_text())
data = Path("/home/dominic/Tudat/laser-ranging-work/data/llr")
points = llr.read_points(
    data / dataset["normal_points"]["file"], dataset["normal_points"]["sha256"]
)
llr.spice.clear_kernels()
for record in dataset["kernels"]:
    path = data / "kernels" / record["file"]
    assert llr.sha256(path) == record["sha256"]
    llr.spice.load_kernel(str(path))
bodies, position = llr.create_environment(points, 30)
point = points[-1]
converter = bodies.get("Earth").rotation_model.angles_calculator.time_scale_converter
epoch = converter.convert_time_object(
    llr.time.utc_scale, llr.time.tdb_scale, point["utc"], position
)
shapiro = llr.ltc.first_order_relativistic_light_time_correction(["Sun", "Earth", "Moon"])
optical = llr.ltc.mendes_pavlis_tropospheric_light_time_correction(point["wavelength_m"])
signature = llr.ltc.LaserTargetSignature(
    [np.zeros(3)], [1.0], internal_optical_path_length=1.4555 * 0.0275
)
target = llr.ltc.laser_target_light_time_correction("Moon", signature)
factory = llr.ltc.light_time_convergence_settings
result = {
    "point_index": point["index"],
    "utc_launch": point["utc_label"],
    "spice_knot_step_s": 30,
    "kernel_path": kernel.__file__,
    "kernel_sha256": llr.sha256(Path(kernel.__file__)),
    "source_sha256": llr.sha256(source),
    "probe_sha256": llr.sha256(Path(__file__)),
    "nh_m": 1.4555 * 0.0275,
    "runs": {},
}
for tolerance in (1e-13, 1e-14):
    llr.ltc.light_time_convergence_settings = (
        lambda failure, iterations, ignored, handling: factory(
            failure, iterations, tolerance, handling
        )
    )
    try:
        atmospheric, atmospheric_events = llr.computed_range(
            point, epoch, bodies, [shapiro, optical]
        )
        complete, complete_events = llr.computed_range(
            point, epoch, bodies, [shapiro, optical, target]
        )
        result["runs"][str(tolerance)] = {
            "atmospheric_range_m": atmospheric,
            "complete_range_m": complete,
            "glass_increment_m": 0.5 * (complete - atmospheric),
            "residual_m": 0.5 * (constants.SPEED_OF_LIGHT * point["observed_rtt_s"] - complete),
            "atmospheric_events": atmospheric_events,
            "complete_events": complete_events,
        }
    except RuntimeError as error:
        result["runs"][str(tolerance)] = {"error": str(error)}
llr.ltc.light_time_convergence_settings = factory
print(json.dumps(result, indent=2))
