"""Compare existing native ephemeris routes at two precise, equal-double epochs."""

import json
from pathlib import Path

import numpy as np

from tudatpy.astro import time_representation as time
from tudatpy.data_input import resource_paths
from tudatpy.dynamics import environment_setup
from tudatpy.dynamics.environment_setup import ephemeris
from tudatpy.interface import spice
from tudatpy.math import interpolators

# Integer knots avoid source-SPICE epoch quantization during table construction.
reference = 765201600
left = time.Time(reference)
right = left + 20e-9
settings = interpolators.lagrange_interpolation(
    8,
    boundary_interpolation=interpolators.throw_exception_at_boundary,
    lagrange_boundary_handling=interpolators.lagrange_no_boundary_interpolation,
)
spice.clear_kernels()
spice.load_kernel(
    str(Path(resource_paths.get_spice_kernel_path()) / "inpop19a_TDB_m100_p100_spice.bsp")
)
results = {
    "left_time_split": left.__getstate__(),
    "right_time_split": right.__getstate__(),
    "dt_seconds": float(right - left),
    "equal_double_epochs": float(left) == float(right),
}

for name in ("Earth", "Moon"):
    states = {}

    for route in ("direct", "time_interpolated"):
        body_settings = environment_setup.BodyListSettings("SSB", "J2000")
        body_settings.add_empty_settings(name)
        body_settings.get(name).ephemeris_settings = (
            ephemeris.direct_spice("SSB", "J2000", name)
            if route == "direct"
            else ephemeris.interpolated_spice(
                reference - 3600, reference + 3600, 60, "SSB", "J2000", settings, name
            )
        )
        body = environment_setup.create_system_of_bodies(body_settings).get(name)
        first = np.asarray(body.state_in_base_frame_from_ephemeris(left)).reshape(6)
        second = np.asarray(body.state_in_base_frame_from_ephemeris(right)).reshape(6)
        expected = first[3:] * float(right - left)
        actual = second[:3] - first[:3]
        states[route] = {
            "first_position_m": first[:3].tolist(),
            "velocity_m_per_s": first[3:].tolist(),
            "actual_displacement_m": actual.tolist(),
            "velocity_dt_m": expected.tolist(),
            "displacement_norm_m": float(np.linalg.norm(actual)),
            "difference_from_velocity_dt_m": (actual - expected).tolist(),
            "difference_from_velocity_dt_norm_m": float(np.linalg.norm(actual - expected)),
        }

    results[name] = states

print(json.dumps(results, indent=2, allow_nan=False))
