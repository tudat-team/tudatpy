"""gcrs_to_itrs rotation model: the analytical rotation-matrix derivative must agree with a numerical derivative of
the same model's rotation matrix.

The derivative used to evaluate the Earth rotation angle at UT1 - 1.42 ms, while the rotation itself uses UT1. That
shift gives an elementwise error of ~7.6e-12 /s (about 48 um/s in station velocity). The derivative still omits the
precession-nutation and UT1 rates (by design), which leaves up to ~4e-12 /s; the limit below sits between the two.
"""

import numpy as np
import pytest

from tudatpy.dynamics import environment_setup
from tudatpy.interface import spice


@pytest.mark.parametrize("t0", [0.0, 1.0e7, 1.0e8, 4.0e8, 7.0e8])
def test_gcrs_to_itrs_rotation_derivative_matches_numerical(t0):
    spice.load_standard_kernels()
    settings = environment_setup.get_default_body_settings(["Earth"], "Earth", "J2000")
    settings.get("Earth").rotation_model_settings = environment_setup.rotation_model.gcrs_to_itrs(
        environment_setup.rotation_model.IAUConventions.iau_2006, "J2000"
    )
    rotation = environment_setup.create_system_of_bodies(settings).get("Earth").rotation_model

    def matrix(t):
        return np.asarray(rotation.body_fixed_to_inertial_rotation(t))

    errors = []
    for t in t0 + np.arange(0.0, 86400.0, 3600.0):
        # fourth-order central difference
        numerical = (
            4.0 * (matrix(t + 5.0) - matrix(t - 5.0)) / 10.0
            - (matrix(t + 10.0) - matrix(t - 10.0)) / 20.0
        ) / 3.0
        analytical = np.asarray(rotation.time_derivative_body_fixed_to_inertial_rotation(t))
        errors.append(np.abs(analytical - numerical).max())
    assert max(errors) < 6.0e-12
