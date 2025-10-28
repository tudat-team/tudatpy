
from tudatpy.data.mpc import BatchMPC, MPC80ColsParser
from tudatpy.data.horizons import HorizonsQuery

from tudatpy.dynamics import environment_setup
from tudatpy.interface import spice

import numpy as np
import pytest
import datetime
import pytest
from tudatpy.astro.time_representation import DateTime
from astroquery.mpc import MPC
import pandas as pd

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

    # batch.filter takes python datetimes in UTC!
    batch.filter(
        epoch_start=datetime.datetime(2017, 1, 1),
        epoch_end=datetime.datetime(2022, 1, 1),
        observatories=["T08"],
    )

    # Horizons Query wants batch_times (or start_epoch, end_epoch) in UTC!!!
    utc_datetimes = batch.table.epochUTC
    batch_times = [DateTime.to_epoch(DateTime.from_python_datetime(t)) for t in utc_datetimes]
    eros = HorizonsQuery(
        query_id="433;", location="T08@399", epoch_list=batch_times, extended_query=True
    )

    # interpolated_observations returns times in TDB!!!
    radec_horizons = eros.interpolated_observations(degrees=False)

    # the retrieved batch.table has time columns: epoch [julian days in UTC], epochUTC [UTC datetime], epochJ2000secondsTDB [TDB seconds]
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


# COMMENTED DUE TO REGULAR TIMEOUT AND FAILURE ON AZURE
# @pytest.mark.parametrize("mpc_code", mpc_codes_test2)
# def test_mpc_coverage(mpc_code):
#     batch_base = BatchMPC()
#     batch_base.get_observations([mpc_code])
#     batch_base.filter(
#         epoch_start=datetime.datetime(2021, 1, 1),
#         epoch_end=datetime.datetime(2022, 1, 1),
#     )
#
#     # properties
#     batch_base.table
#     batch_base.observatories
#     batch_base.space_telescopes
#     batch_base.MPC_objects
#     batch_base.size
#     batch_base.bands
#     batch_base.epoch_start
#     batch_base.epoch_end
#     len(batch_base)
#
#     # addition
#     batch2 = BatchMPC()
#     batch2.get_observations([1])
#     batch2.filter(
#         epoch_start=datetime.datetime(2021, 1, 1),
#         epoch_end=datetime.datetime(2022, 1, 1),
#     )
#     batch3 = batch_base + batch2
#
#     # copy
#     batch3copy = batch3.copy()
#
#     # from_pandas + from_astropy
#     batch4 = BatchMPC()
#     batch5 = BatchMPC()
#
#     batch4.from_astropy(astroquery_MPC.get_observations(mpc_code))
#     batch5.from_pandas(batch_base._table)  # type: ignore
#
#     # plotting
#     batch_base.plot_observations_temporal()
#     batch_base.plot_observations_sky()
#     batch_base.plot_observations_sky(projection="hammer")
#     batch_base.plot_observations_sky(projection="mollweide")
#     batch_base.plot_observations_sky(projection="lambert")
#
#     # obs_table
#     batch_base.observatories_table(only_in_batch=False)
#     batch_base.observatories_table(only_space_telescopes=True)
#     batch_base.observatories_table(exclude_space_telescopes=True)
#     batch_base.observatories_table(include_positions=True)
#
#     # summary
#     batch_base.summary()

def test_80cols_line_parser():

    batch = BatchMPC()
    batch.get_observations(['3I', 433, 134341, '2025 FA22'],  id_types = ['comet_number', 'asteroid_number', 'asteroid_number', 'asteroid_designation'])

    # Observation Lines are taken from astroquery.MPC.get_observations with the 'get_mpcformat = True' flag.
    line_atlas = '0003I         S2025 05 08.51765919 12 35.590-18 42 21.35         21.57VVER063C5' # Interstellar Comet
    line_eros = '00433         A1893 10 29.4132  06 08 59.32 +53 39 04.2                 HA053802' # Asteroid/Minor Planet
    line_charon = 'D4341J79M00A*4A1979 06 25.66181 20 27 06.64 -15 37 11.5          19.0   M4986413' # Natural Satellite
    line_2025FA22 = '     K25F22A  C2025 10 13.24277 00 20 45.76 +25 53 06.1          18.3 RrET147718'
    MPC_parser = MPC80ColsParser()
    parsed_line_atlas = MPC_parser.parse_80cols_identification_fields(line_atlas)
    parsed_line_eros = MPC_parser.parse_80cols_identification_fields(line_eros)
    parsed_line_charon = MPC_parser.parse_80cols_identification_fields(line_charon)
    parsed_line_2025FA22 = MPC_parser.parse_80cols_identification_fields(line_2025FA22)

    assert(parsed_line_2025FA22['desig'] == batch.MPC_objects[-1]) # at the time of writing, 2025FA22 does not have a number. We test the designation.
    assert [int(parsed_line_atlas['number']), int(parsed_line_eros['number']), int(parsed_line_charon['number'])] == [int(x) for x in batch.MPC_objects[:-1]]

def test_parse_80cols_file():
    batch = BatchMPC()
    batch.get_observations([433])
    batch.filter(epoch_start = datetime.datetime(2021, 6, 7, 00, 4), epoch_end =  datetime.datetime(2021, 6, 7, 16, 4,2))
    MPC_parser = MPC80ColsParser()
    file_path = '/Users/lgisolfi/CLionProjects/tudatpy_examples/estimation/data/eros_obs.txt'
    table_output = MPC_parser.parse_80cols_file(file_path)

    epochs1 = pd.to_datetime(table_output['epoch_utc']).to_numpy()
    epochs2 = batch.table['epochUTC'].to_numpy()
    # Get difference in seconds
    diff = np.sort(epochs1) - np.sort(epochs2)
    diff_seconds = diff / np.timedelta64(1, 's')

    tol = 5e-5 # not completely sure why some are zero and some are not.
    assert not (diff_seconds > tol).any()

test_compare_mpc_horizons_eph()