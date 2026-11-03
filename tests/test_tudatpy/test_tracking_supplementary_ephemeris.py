"""Supplementary spacecraft ephemerides extrapolate silently at their boundaries."""

import numpy as np
import pytest

from tudatpy.astro import time_representation as time
from tudatpy.data_input.tracking_data import (
    TrackingSupplementaryData,
    TranslationalStateSupplementaryData,
)
from tudatpy.dynamics import environment_setup
from tudatpy.estimation import observations


@pytest.mark.parametrize("existing_ephemeris", [False, True])
@pytest.mark.parametrize("time_scale", ["TDB", "UTC"])
def test_supplementary_ephemeris_extrapolates_without_boundary_warnings(
    existing_ephemeris, time_scale, capfd
):
    settings = environment_setup.BodyListSettings("SSB", "J2000")
    settings.add_empty_settings("Earth")
    earth_state = np.array([100.0, 200.0, 300.0, 0.0, 0.0, 0.0])
    settings.get("Earth").ephemeris_settings = environment_setup.ephemeris.constant(
        earth_state, "SSB", "J2000"
    )
    start, end = 5.0e8, 5.0e8 + 10.0
    initial_state = np.array([1.0, 2.0, 3.0, 0.1, 0.2, 0.3])
    final_state = initial_state + np.arange(1.0, 7.0)
    if existing_ephemeris:
        settings.add_empty_settings("Gaia")
        settings.get("Gaia").ephemeris_settings = environment_setup.ephemeris.tabulated(
            {start + i: initial_state for i in range(8)}, "Earth", "J2000"
        )
    bodies = environment_setup.create_system_of_bodies(settings)
    supplementary = TrackingSupplementaryData("Gaia", "")
    supplementary.translational_state_supplementary_data = TranslationalStateSupplementaryData(
        state_history={start: initial_state, end: final_state},
        frame_origin="Earth",
        is_velocity_defined=True,
        time_scale=time_scale,
        frame_orientation="J2000",
    )
    observations.set_tracking_supplementary_data_in_bodies(bodies, [supplementary])
    if time_scale == "UTC":
        converter = time.default_time_scale_converter()
        start = converter.convert_time(time.utc_scale, time.tdb_scale, start)
        end = converter.convert_time(time.utc_scale, time.tdb_scale, end)
    capfd.readouterr()
    # Exercise fractions of a double ULP and genuine extrapolation on both ends.
    for epoch in [
        time.Time(start) - time.Time(1.0e-9),
        time.Time(end) + time.Time(1.0e-9),
        time.Time(start - 1.0),
        time.Time(end + 1.0),
    ]:
        fraction = float(epoch - time.Time(start)) / (end - start)
        expected = initial_state + fraction * (final_state - initial_state)
        np.testing.assert_allclose(
            bodies.get("Gaia").state_in_base_frame_from_ephemeris(epoch),
            expected + earth_state,
            rtol=0.0,
            atol=1.0e-12,
        )
    assert not capfd.readouterr().err
