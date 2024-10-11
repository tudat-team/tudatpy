import sys
sys.path.insert(0, "/home/dominic/Tudat/tudat-bundle/tudat-bundle/cmake-build-debug/tudatpy")


from tudatpy.data.mpc import BatchMPC
from tudatpy.data.horizons import HorizonsQuery

from tudatpy.numerical_simulation import environment_setup
from tudatpy.interface import spice

import numpy as np
import pytest
import datetime

from astroquery.mpc import MPC as astroquery_MPC


spice.load_standard_kernels()

# coverage = 88%
# TESTS DO NOT CHECK/VALIDATE:
# positions of observatories.

# Parameterised inputs
mpc_codes_test = [222, 999]
mpc_codes_test2 = [3]

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


# for the weights tests 
observatory_set_single = ["M22"]
observatory_set_multi = ["K19", "D67", "089", "706"]
weights_test_combinations = [
    (observatory_set_single, True),  # just one obs
    (observatory_set_single, False),
    (observatory_set_multi, False),
    (None, False),  # all data
]


#@pytest.mark.parametrize("inp,expected", get_observations_input)
#def test_BatchMPC_getobservations(inp, expected):
#    query = BatchMPC()
#    query.get_observations(inp)
#    assert set(query.MPC_objects) == expected


#@pytest.mark.parametrize("inp,errtype,errvalue", get_observations_input2)
#def test_BatchMPC_getobservations2(inp, errtype, errvalue):
#    query = BatchMPC()
#    with pytest.raises(Exception) as exc_info:
#        query.get_observations(inp)
#
#    assert exc_info.type is errtype
#    assert str(exc_info.value) == errvalue


@pytest.mark.parametrize("mpc_code", mpc_codes_test)
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
        apply_star_catalog_debias=False
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


@pytest.mark.parametrize("mpc_code", mpc_codes_test)
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
        apply_star_catalog_debias=False
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


@pytest.mark.parametrize("mpc_code", mpc_codes_test2)
def test_mpc_coverage(mpc_code):
    batch_base = BatchMPC()
    batch_base.get_observations([mpc_code])
    batch_base.filter(
        epoch_start=datetime.datetime(2021, 1, 1),
        epoch_end=datetime.datetime(2022, 1, 1),
    )

    # properties
    batch_base.table
    batch_base.observatories
    batch_base.space_telescopes
    batch_base.MPC_objects
    batch_base.size
    batch_base.bands
    batch_base.epoch_start
    batch_base.epoch_end
    len(batch_base)

    # addition
    batch2 = BatchMPC()
    batch2.get_observations([1])
    batch2.filter(
        epoch_start=datetime.datetime(2021, 1, 1),
        epoch_end=datetime.datetime(2022, 1, 1),
    )
    batch3 = batch_base + batch2

    # copy
    batch3copy = batch3.copy()

    # from_pandas + from_astropy
    batch4 = BatchMPC()
    batch5 = BatchMPC()

    batch4.from_astropy(astroquery_MPC.get_observations(mpc_code))
    batch5.from_pandas(batch_base._table)  # type: ignore

    # plotting
    batch_base.plot_observations_temporal()
    batch_base.plot_observations_sky()
    batch_base.plot_observations_sky(projection="hammer")
    batch_base.plot_observations_sky(projection="mollweide")
    batch_base.plot_observations_sky(projection="lambert")

    # obs_table
    batch_base.observatories_table(only_in_batch=False)
    batch_base.observatories_table(only_space_telescopes=True)
    batch_base.observatories_table(exclude_space_telescopes=True)
    batch_base.observatories_table(include_positions=True)

    # summary
    batch_base.summary()
