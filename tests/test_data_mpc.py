from tudatpy.data.mpc import BatchMPC
from tudatpy.data.horizons import HorizonsQuery

from tudatpy.numerical_simulation import environment_setup
from tudatpy.interface import spice

import numpy as np
import pytest
import datetime

spice.load_standard_kernels()

# TESTS DO NOT COVER:
# properties
# plot_observations_temporal()
# plot_observations_sky()
# summary()

# TESTS DO NOT CHECK:
# positions of observatories.


get_observations_input = [
    ([999, 222], {"999", "222"}),
    ([222, "C/2012 S1"], {"222", "2012 S1"}),
]
get_observations_input2 = [
    (222, ValueError, "MPCcodes parameter must be list of integers/strings"),
    (
        [222, 1.0],
        ValueError,
        "All codes in the MPCcodes parameter must be integers or string",
    ),
]


@pytest.mark.parametrize("inp,expected", get_observations_input)
def test_BatchMPC_getobservations(inp, expected):
    """Check if get observations works"""
    query = BatchMPC()
    query.get_observations(inp)
    assert set(query.MPC_objects) == expected


@pytest.mark.parametrize("inp,errtype,errvalue", get_observations_input2)
def test_BatchMPC_getobservations2(inp, errtype, errvalue):
    """Check if get observations works, this time with errors"""
    query = BatchMPC()
    with pytest.raises(Exception) as exc_info:
        query.get_observations(inp)

    assert exc_info.type is errtype
    assert str(exc_info.value) == errvalue


filter_test_input = [
    (
        999,
        datetime.datetime(2022, 1, 1),
        datetime.datetime(2023, 1, 1),
        ["C51"],
        ["T08", "T05", "U55"],
        684,
        264,
        241,
        141,
    ),
    (
        222,
        datetime.datetime(2022, 1, 1),
        datetime.datetime(2023, 1, 1),
        ["C51"],
        ["T08", "T05", "U55"],
        575,
        214,
        209,
        7,
    ),
]


@pytest.mark.parametrize(
    "mpc_code,datestart,datestop,obs_exlude,obs_include,exp_size_start,exp_size_stop,exp_size_excl,exp_size_incl",
    filter_test_input,
)
def test_BatchMPC_filter(
    mpc_code,
    datestart,
    datestop,
    obs_exlude,
    obs_include,
    exp_size_start,
    exp_size_stop,
    exp_size_excl,
    exp_size_incl,
):
    """Test if BatchMPC filtering works correctly"""
    b = BatchMPC()

    b.get_observations([mpc_code])
    b.filter(epoch_start=datestart)
    assert b.size == exp_size_start
    assert len(b._table) == exp_size_start

    b.filter(epoch_end=datestop)
    assert b.size == exp_size_stop
    assert len(b._table) == exp_size_stop

    b.filter(observatories_exclude=obs_exlude)
    assert b.size == exp_size_excl
    assert len(b._table) == exp_size_excl

    b.filter(observatories=obs_include)
    assert b.size == exp_size_incl
    assert len(b._table) == exp_size_incl

    # test not in place:
    b.filter(observatories=[], in_place=False)
    assert b.size == exp_size_incl
    assert len(b._table) == exp_size_incl


to_tudat_inp = [222, 999]


@pytest.mark.parametrize("mpc_code", to_tudat_inp)
def test_BatchMPC_to_tudat(mpc_code):
    """Check if observatory table matches observation_collection"""
    query = BatchMPC()
    query.get_observations([mpc_code])
    query.filter(observatories=["T05", "T08"])

    # table values are sorted for easier comparison
    query._table = query._table.sort_values(["observatory", "epochJ2000secondsTDB"])

    RADEC = query.table.loc[:, ["RA", "DEC"]].to_numpy().T
    times = query.table.loc[:, ["epochJ2000secondsTDB"]].to_numpy().T[0]
    times = np.array([times, times])  # concat times are doubled due to RA + DEC

    # to_tudat needs a system of bodies with earth in it as input
    bodies_to_create = [
        "Earth",
    ]
    global_frame_origin = "SSB"
    global_frame_orientation = "J2000"
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation
    )
    bodies = environment_setup.create_system_of_bodies(body_settings)

    observation_collection = query.to_tudat(
        bodies=bodies,
        included_satellites=None,
    )

    # reshape to [2, ...] where 2 is RA + DEC
    obscol_RADEC = (np.array(observation_collection.concatenated_observations)).reshape(
        2, -1, order="F"
    )
    obscol_times = (np.array(observation_collection.concatenated_times)).reshape(
        2, -1, order="F"
    )

    # max error between the two should be zero
    assert (np.max(obscol_times - times)) == pytest.approx(0.00)
    assert (np.max(obscol_RADEC - RADEC)) == pytest.approx(0.00)


to_tudat_inp2 = [222, 999]


@pytest.mark.parametrize("mpc_code", to_tudat_inp2)
def test_BatchMPC_to_tudat_with_satelite(mpc_code):
    """Check if observatory table matches observation_collection"""
    query = BatchMPC()
    query.get_observations([mpc_code])
    query.filter(observatories=["C51"])

    # table values are sorted for easier comparison
    query._table = query._table.sort_values(["observatory", "epochJ2000secondsTDB"])

    RADEC = query.table.loc[:, ["RA", "DEC"]].to_numpy().T
    times = query.table.loc[:, ["epochJ2000secondsTDB"]].to_numpy().T[0]
    times = np.array([times, times])  # concat times are doubled due to RA + DEC

    # to_tudat needs a system of bodies with earth in it as input
    bodies_to_create = [
        "Earth",
    ]
    global_frame_origin = "SSB"
    global_frame_orientation = "J2000"
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation
    )
    bodies = environment_setup.create_system_of_bodies(body_settings)
    bodies.create_empty_body("Wise")

    observation_collection = query.to_tudat(
        bodies=bodies,
        included_satellites={"C51": "Wise"},
    )

    # reshape to [2, ...] where 2 is RA + DEC
    obscol_RADEC = (np.array(observation_collection.concatenated_observations)).reshape(
        2, -1, order="F"
    )
    obscol_times = (np.array(observation_collection.concatenated_times)).reshape(
        2, -1, order="F"
    )

    # max error between the two should be zero
    assert (np.max(obscol_times - times)) == pytest.approx(0.00)
    assert (np.max(obscol_RADEC - RADEC)) == pytest.approx(0.00)


def test_compare_mpc_horizons_eph():
    """Compares true observations from BatchMPC to interpolated simulated RA/DEC from JPL Horizons"""
    batch = BatchMPC()
    batch.get_observations([433])
    batch.filter(
        epoch_start=datetime.datetime(2017, 1, 1),
        epoch_end=datetime.datetime(2022, 1, 1),
        observatories=["T08"],
    )
    batch_times = batch.table.epochJ2000secondsTDB.to_list()

    eros = HorizonsQuery(
        query_id="433;", location="T08@399", epoch_list=batch_times, extended_query=True
    )

    radec_horizons = eros.interpolated_observations(degrees=False)
    radec_mpc = batch.table.loc[:, ["epochJ2000secondsTDB", "RA", "DEC"]].reset_index(
        drop=True
    )

    diff = (radec_horizons - radec_mpc).to_numpy()
    diff = np.abs(diff).max(axis=0)

    time_diff = diff[0]
    RA_diff = diff[1]
    DEC_diff = diff[2]

    assert time_diff < 1e-3
    assert RA_diff < 1e-5
    assert DEC_diff < 1e-5
