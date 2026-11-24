from tudatpy.data_input.environment_data.horizons import HorizonsBatch, HorizonsQuery
from tudatpy.dynamics.environment_setup.ephemeris import (
    add_horizons_batch_ephemerides,
    jpl_horizons_from_query,
    jpl_horizons,
)
from tudatpy.dynamics import environment_setup


from tudatpy.data_input.environment_data import spice
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
import datetime
import inspect
from astroquery.jplhorizons import Horizons

spice.load_standard_kernels()

pytestmark = pytest.mark.remote_data

# 87% test coverage for tudatpy/data/horizons.py. Remainder is rare user input validation
# 100% test coverage for tudatpy/numerical_simulation/environment_setup/ephemeris/horizons_wrapper.py.


def test_compare_horizons_spice(pos_tolerance=1e3, vel_tolerance=1e-3):
    """Compares the states from an array of times between HorizonsQuery and SPICE"""
    start = 20 * 365 * 86400
    end = 23 * 365 * 86400

    npoints = 100

    times = np.linspace(start, end, npoints)

    query = HorizonsQuery(
        query_id="299",
        location="500@399",
        epoch_list=list(times),
        extended_query=True,
    )
    ref = "ECLIPJ2000"

    horizons_states = query.cartesian(frame_orientation=ref)[:, 1:]
    spice_states = []

    for time in times:
        # from spice
        state_spice = spice.get_body_cartesian_state_at_epoch("Venus", "Earth", ref, "NONE", time)
        spice_states.append(state_spice)

    diff = np.abs(horizons_states - spice_states)
    max_diff_pos = diff.max(axis=0)[0:3].max()
    max_diff_vel = diff.max(axis=0)[3:6].max()

    assert max_diff_pos < pos_tolerance
    assert max_diff_vel < vel_tolerance


targets = [
    ("433;", "500@399", "J2000", "Earth"),
    ("299", "500@SSB", "J2000", "SSB"),
    ("-64", "500@SSB", "J2000", "SSB"),
    ("C/1995 O1", "500@SSB", "J2000", "SSB"),
]

user_input = [
    # start end stepsize + epoch formats
    (
        datetime.datetime(2022, 1, 1, 0, 0, 0),
        datetime.datetime(2022, 2, 1, 0, 0, 0),
        "5d",
        None,
        False,
        does_not_raise(),
    ),
    (
        datetime.datetime(2022, 1, 1, 0, 0, 0),
        datetime.datetime(2022, 5, 1, 0, 0, 0),
        "1m",
        None,
        True,
        does_not_raise(),
    ),
    (
        694267200.000,
        datetime.datetime(2022, 5, 1, 0, 0, 0),
        "1m",
        None,
        True,
        does_not_raise(),
    ),
    (
        694267200.000,
        datetime.datetime(2022, 2, 1, 0, 0, 0),
        "5d",
        None,
        False,
        does_not_raise(),
    ),
    (
        datetime.datetime(2022, 1, 1, 0, 0, 0),
        696945600.000,
        "5d",
        None,
        False,
        does_not_raise(),
    ),
    (
        694267200.000,
        704635200.000000,
        "1m",
        None,
        True,
        does_not_raise(),
    ),
    # start end num steps
    (
        datetime.datetime(2022, 1, 1, 0, 0, 0),
        datetime.datetime(2022, 2, 1, 0, 0, 0),
        "20",
        None,
        False,
        does_not_raise(),
    ),
    (
        datetime.datetime(2022, 1, 1, 0, 0, 0),
        datetime.datetime(2022, 2, 1, 0, 0, 0),
        "100000",
        None,
        True,
        pytest.raises(NotImplementedError),
    ),
    # list of times
    (
        None,
        None,
        None,
        [(x * 86400) + 694267200.0 for x in range(40)],
        False,
        does_not_raise(),
    ),
    (
        None,
        None,
        None,
        [(x * 86400) + 694267200.0 for x in range(40)],
        True,
        does_not_raise(),
    ),
    (
        None,
        None,
        None,
        [(x * 86400) + 694267200.0 for x in range(130)],
        False,
        pytest.raises(ValueError),
    ),
    (
        None,
        None,
        None,
        [(x * 86400) + 694267200.0 for x in range(130)],
        True,
        does_not_raise(),
    ),
]

user_input_short = [
    (
        datetime.datetime(2022, 1, 1, 0, 0, 0),
        datetime.datetime(2022, 2, 1, 0, 0, 0),
        "5d",
        None,
        False,
        does_not_raise(),
    ),
    (
        None,
        None,
        None,
        [(x * 86400) + 694267200.0 for x in range(90)],
        True,
        does_not_raise(),
    ),
]


@pytest.mark.parametrize("query_id,location,frame_orientation,frame_origin", targets)
@pytest.mark.parametrize(
    "epoch_start,epoch_end,epoch_step,epoch_list,extended_query,expectation", user_input
)
def test_JPL_user_input(
    query_id,
    location,
    frame_orientation,
    frame_origin,
    epoch_start,
    epoch_end,
    epoch_step,
    epoch_list,
    extended_query,
    expectation,
):
    """Check if input handling is working as expected"""
    with expectation as exc_info:
        query = HorizonsQuery(
            query_id=query_id,
            location=location,
            epoch_start=epoch_start,
            epoch_end=epoch_end,
            epoch_step=epoch_step,
            epoch_list=epoch_list,
            extended_query=extended_query,
        )


@pytest.mark.parametrize("query_id,location,frame_orientation,frame_origin", targets)
@pytest.mark.parametrize(
    "epoch_start,epoch_end,epoch_step,epoch_list,extended_query,expectation",
    user_input_short,
)
def test_JPL_methods(
    query_id,
    location,
    frame_orientation,
    frame_origin,
    epoch_start,
    epoch_end,
    epoch_step,
    epoch_list,
    extended_query,
    expectation,
):
    """Check if the HorizonsQuery methods and parameters are working"""
    with expectation as exc_info:
        query = HorizonsQuery(
            query_id=query_id,
            location=location,
            epoch_start=epoch_start,
            epoch_end=epoch_end,
            epoch_step=epoch_step,
            epoch_list=epoch_list,
            extended_query=extended_query,
        )

        query.cartesian(frame_orientation=frame_orientation)
        jpl_horizons_from_query(
            query,
            frame_orientation=frame_orientation,
            frame_origin=frame_origin,
        )
        query.interpolated_observations(reference_system=frame_orientation)

        query.MPC_number
        query.designation
        query.name


@pytest.mark.parametrize("query_id,location,frame_orientation,frame_origin", targets)
def test_hybrid_function_wrapper(query_id, location, frame_orientation, frame_origin):
    """Coverage test for the hybrid function, no comparisons just seeing if it runs"""
    juice_eph_settings = jpl_horizons(
        horizons_query=query_id,
        horizons_location=location,
        frame_origin=frame_origin,  # tudat frame origin and orientation
        frame_orientation=frame_orientation,
        epoch_start=datetime.datetime(2020, 1, 1, 0, 0, 0),
        epoch_end=datetime.datetime(2023, 1, 1, 0, 0, 0),
        epoch_step="1d",
        extended_query=True,
    )

    bodies_to_create = [
        "Earth",
    ]
    global_frame_origin = frame_origin
    global_frame_orientation = frame_orientation
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation
    )

    body_settings.add_empty_settings("Eros")

    body_settings.get("Eros").ephemeris_settings = juice_eph_settings


@pytest.mark.parametrize("query_id,location,frame_orientation,frame_origin", targets)
@pytest.mark.parametrize(
    "epoch_start,epoch_end,epoch_step,epoch_list,extended_query,expectation",
    user_input_short,
)
def test_horizons_batch(
    query_id,
    location,
    frame_orientation,
    frame_origin,
    epoch_start,
    epoch_end,
    epoch_step,
    epoch_list,
    extended_query,
    expectation,
):
    batch = HorizonsBatch(
        query_id_list=[query_id, "123", "456"],
        location=location,
        epoch_list=epoch_list,
        epoch_start=epoch_start,
        epoch_end=epoch_end,
        epoch_step=epoch_step,
        extended_query=extended_query,
    )

    bodies_to_create = [
        "Earth",
    ]
    global_frame_origin = frame_origin
    global_frame_orientation = frame_orientation
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation
    )

    add_horizons_batch_ephemerides(
        batch,
        body_settings=body_settings,
        frame_origin=frame_origin,
        frame_orientation=frame_orientation,
    )

    batch.names


@pytest.mark.skipif(
    "vector_table" not in inspect.signature(Horizons.vectors_async).parameters,
    reason="Requires astroquery >= 0.4.12",
)
def test_cartesian_uncertainties():
    """Checks that the trace of the position and velocity covariance is the same in
    all uncertainty frames, and that bodies without covariance give NaN"""
    query = HorizonsQuery(
        query_id="1;",
        location="500@10",
        epoch_start=datetime.datetime(2020, 1, 1),
        epoch_end=datetime.datetime(2020, 1, 10),
        epoch_step="1d",
    )
    sigmas = {
        frame: query.cartesian_uncertainties(uncertainty_frame=frame)
        for frame in ["xyz", "acn", "rtn"]
    }
    assert np.all(sigmas["xyz"][:, 1:] > 0)

    position_trace = np.sum(sigmas["xyz"][:, 1:4] ** 2, axis=1)
    velocity_trace = np.sum(sigmas["xyz"][:, 4:7] ** 2, axis=1)
    for frame in ["acn", "rtn"]:
        np.testing.assert_array_equal(sigmas[frame][:, 0], sigmas["xyz"][:, 0])
        np.testing.assert_allclose(
            np.sum(sigmas[frame][:, 1:4] ** 2, axis=1), position_trace, rtol=1e-6
        )
        np.testing.assert_allclose(
            np.sum(sigmas[frame][:, 4:7] ** 2, axis=1), velocity_trace, rtol=1e-6
        )

    query = HorizonsQuery(query_id="499", location="500@10", epoch_list=[0.0])
    assert np.all(np.isnan(query.cartesian_uncertainties()[:, 1:]))
