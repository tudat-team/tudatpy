from tudatpy.data.horizons import HorizonsQuery
from tudatpy.numerical_simulation.environment_setup.ephemeris import jpl_horizons
from tudatpy.numerical_simulation import environment_setup


from tudatpy.interface import spice
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
import datetime


spice.load_standard_kernels()

# 81% test coverage, Horizons Batch not tested


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
        state_spice = spice.get_body_cartesian_state_at_epoch(
            "Venus", "Earth", ref, "NONE", time
        )
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
        datetime.datetime(2022, 1, 1),
        datetime.datetime(2022, 2, 1),
        "5d",
        None,
        False,
        does_not_raise(),
    ),
    (
        datetime.datetime(2022, 1, 1),
        datetime.datetime(2022, 5, 1),
        "1m",
        None,
        True,
        does_not_raise(),
    ),
    (
        694267200.000,
        datetime.datetime(2022, 5, 1),
        "1m",
        None,
        True,
        does_not_raise(),
    ),
    (
        694267200.000,
        datetime.datetime(2022, 2, 1),
        "5d",
        None,
        False,
        does_not_raise(),
    ),
    (
        datetime.datetime(2022, 1, 1),
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
        datetime.datetime(2022, 1, 1),
        datetime.datetime(2022, 2, 1),
        "20",
        None,
        False,
        does_not_raise(),
    ),
    (
        datetime.datetime(2022, 1, 1),
        datetime.datetime(2022, 2, 1),
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
        datetime.datetime(2022, 1, 1),
        datetime.datetime(2022, 2, 1),
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
        query.create_ephemeris_tabulated(
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
        epoch_start=datetime.datetime(2020, 1, 1),
        epoch_end=datetime.datetime(2023, 1, 1),
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
