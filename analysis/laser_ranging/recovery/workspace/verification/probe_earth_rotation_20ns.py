"""Locate the 20ns finite-difference discrepancy without changing accuracy settings."""

import ctypes, json, math
from pathlib import Path
import numpy as np
import tudatpy.kernel
from tudatpy.dynamics import environment_setup as setup
from tudatpy.astro import time_representation as time

s = setup.BodyListSettings("SSB", "GCRS")
s.add_empty_settings("Earth")
s.get("Earth").rotation_model_settings = setup.rotation_model.gcrs_to_itrs()
r = setup.create_system_of_bodies(s).get("Earth").rotation_model
a = time.Time(765201600)
b = a + 20e-9
angles_a, ua = r.angles_calculator.get_gcrs_to_itrs_rotation_angles(a)
angles_b, ub = r.angles_calculator.get_gcrs_to_itrs_rotation_angles(b)
f = ctypes.CDLL(tudatpy.kernel.__file__).iauEra00
f.argtypes = [ctypes.c_double, ctypes.c_double]
f.restype = ctypes.c_double


def era(epoch):
    """Call the already-loaded SOFA symbol with the same native split-JD convention."""
    hours, seconds = epoch.__getstate__()
    return f(2451545 + hours // 24, ((hours % 24) * 3600 + seconds) / 86400)


era_diff = math.remainder(era(ub) - era(ua), 2 * math.pi)
spin = 2 * math.pi * 1.00273781191135448 / 86400
station = np.array([6378137.0, 0.0, 0.0])
ra = r.body_fixed_to_inertial_rotation(a)
rb = r.body_fixed_to_inertial_rotation(b)
difference = (rb - ra) @ station
prediction = r.time_derivative_body_fixed_to_inertial_rotation(a) @ station * float(b - a)
quantized_era_prediction = np.cross(np.array([0.0, 0.0, 1.0]), ra @ station) * era_diff
result = {
    "tdb_interval_s": float(b - a),
    "ut1_interval_s": float(ub - ua),
    "era_interval_rad": era_diff,
    "continuous_era_interval_rad": spin * float(ub - ua),
    "era_interval_discrepancy_equatorial_m": (era_diff - spin * float(ub - ua)) * 6378137,
    "station_displacement_m": difference.tolist(),
    "derivative_prediction_m": prediction.tolist(),
    "max_discrepancy_m": float(max(abs(difference - prediction))),
    "orientation_angle_differences": (np.asarray(angles_b) - np.asarray(angles_a))
    .reshape(-1)
    .tolist(),
    "native_era_quantization_equivalent_time_s": era_diff / spin - float(ub - ua),
    "kernel": tudatpy.kernel.__file__,
    "ut1_split_a": ua.__getstate__(),
    "ut1_split_b": ub.__getstate__(),
}
Path("/home/dominic/Tudat/laser-ranging-work/logs/earth-rotation-20ns-probe.json").write_text(
    json.dumps(result, indent=2) + "\n"
)
print(json.dumps(result, indent=2))
