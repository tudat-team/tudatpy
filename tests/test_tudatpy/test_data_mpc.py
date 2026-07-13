from tudatpy.data_access.tracking.mpc import BatchMPC, filter_augmented_optical_table
from tudatpy.data_access.environment.horizons import HorizonsQuery
from tudatpy.estimation.observations import create_observation_collection
from tudatpy.astro import time_representation
from tudatpy.dynamics import environment_setup
from tudatpy.dynamics.environment_setup import ground_station
from tudatpy.interface import spice
import numpy as np
import datetime
import pytest

spice.load_standard_kernels()
_TIME_SCALE_CONVERTER = time_representation.default_time_scale_converter()

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


def _utc_seconds_to_tdb(epoch_seconds_utc):
    return np.array(
        [
            _TIME_SCALE_CONVERTER.convert_time(
                input_scale=time_representation.utc_scale,
                output_scale=time_representation.tdb_scale,
                input_value=float(epoch),
            )
            for epoch in epoch_seconds_utc
        ]
    )


def _batch_from_augmented_table(table) -> BatchMPC:
    batch = BatchMPC()
    batch.from_pandas(table, in_degrees=False)
    return batch


def _sorted_table(batch: BatchMPC):
    return batch.table.sort_values(["observatory", "epoch_seconds_UTC"]).reset_index(drop=True)


def _create_earth_bodies():
    body_settings = environment_setup.get_default_body_settings(["Earth"], "SSB", "J2000")
    body_settings.get("Earth").ground_station_settings = ground_station.optical_telescope_stations()
    return environment_setup.create_system_of_bodies(body_settings)


def _create_observation_collection_from_batch(batch: BatchMPC, bodies):
    tracking_data, supplementary_data = batch.to_tracking_dataset(
        add_star_catalog_corrections=False
    )
    assert supplementary_data == []
    return create_observation_collection(tracking_data, bodies)


def _assert_tracking_dataset_matches_batch(batch: BatchMPC) -> None:
    table = _sorted_table(batch)
    batch._table = table
    tracking_data, supplementary_data = batch.to_tracking_dataset(
        add_star_catalog_corrections=False
    )

    groups = list(table.groupby(["number", "observatory"], sort=True))
    assert supplementary_data == []
    assert len(tracking_data) == len(groups)

    for data_object, ((target, observatory), group) in zip(tracking_data, groups):
        assert data_object.observable_type == "AngularPosition"
        assert data_object.reference_link_end == "receiver"
        assert data_object.time_scale == "UTC"
        assert data_object.link_ends == [
            ((str(target), ""), "transmitter"),
            (("Earth", str(observatory)), "receiver"),
        ]

        expected_observations = group.loc[:, ["RA", "DEC"]].to_numpy()
        expected_epochs = group["epoch_seconds_UTC"].to_numpy()
        actual_observations = np.array(data_object.observations)
        actual_epochs = np.array([epoch.to_float() for epoch in data_object.epochs])

        assert np.max(np.abs(actual_observations - expected_observations)) == pytest.approx(0.0)
        assert np.max(np.abs(actual_epochs - expected_epochs)) == pytest.approx(0.0)


# @pytest.mark.parametrize("inp,expected", get_observations_input)
# def test_BatchMPC_getobservations(inp, expected):
#    query = BatchMPC()
#    query.get_observations(inp)
#    assert set(query.MPC_objects) == expected


# @pytest.mark.parametrize("inp,errtype,errvalue", get_observations_input2)
# def test_BatchMPC_getobservations2(inp, errtype, errvalue):
#    query = BatchMPC()
#    with pytest.raises(Exception) as exc_info:
#        query.get_observations(inp)
#
#    assert exc_info.type is errtype
#    assert str(exc_info.value) == errvalue


@pytest.mark.parametrize("mpc_code", mpc_codes_test)
def test_to_tracking_dataset_preserves_full_mpc_table(mpc_code):
    """Check if the MPC table is preserved when converted to TrackingData."""
    query = BatchMPC()
    query.get_observations([mpc_code])
    _assert_tracking_dataset_matches_batch(query)


@pytest.mark.parametrize("mpc_code", mpc_codes_test)
def test_mpc_custom_name_metadata(mpc_code):
    batch = BatchMPC()
    # Use a real MPC code but provide a custom name
    custom_name = "Death_Star"
    batch.get_observations([mpc_code], custom_name=custom_name)

    # Verify the custom_name column exists and has the correct value
    assert "custom_name" in batch.table.columns
    assert (batch.table["custom_name"] == custom_name).all()

    # Verify MPC_objects returns the custom name instead of mpc code
    assert custom_name in batch.MPC_objects
    assert mpc_code not in batch.MPC_objects


@pytest.mark.parametrize("mpc_code", mpc_codes_test)
def test_to_tracking_dataset_preserves_filtered_observatories(mpc_code):
    """Check if dataframe-filtered MPC observations survive TrackingData conversion."""
    query = BatchMPC()
    query.get_observations([mpc_code])
    query = _batch_from_augmented_table(
        filter_augmented_optical_table(query.table, observatories=["T05", "T08"])
    )
    assert set(query.table["observatory"].unique()) <= {"T05", "T08"}
    _assert_tracking_dataset_matches_batch(query)


@pytest.mark.parametrize("mpc_code", mpc_codes_test)
def test_to_tracking_dataset_handles_alphanumeric_observatory_codes(mpc_code):
    """Check if alphanumeric MPC observatory codes survive TrackingData conversion."""
    query = BatchMPC()
    query.get_observations([mpc_code])
    query = _batch_from_augmented_table(
        filter_augmented_optical_table(query.table, observatories=["F51"])
    )
    assert set(query.table["observatory"].unique()) == {"F51"}
    _assert_tracking_dataset_matches_batch(query)


def test_tracking_dataset_can_create_observation_collection():
    """Check the integration boundary from UTC TrackingData to ObservationCollection."""
    query = BatchMPC()
    query.get_observations([222])
    query = _batch_from_augmented_table(
        filter_augmented_optical_table(query.table, observatories=["T05"])
    )
    query._table = _sorted_table(query)

    observation_collection = _create_observation_collection_from_batch(
        query, _create_earth_bodies()
    )

    expected_observations = query.table.loc[:, ["RA", "DEC"]].to_numpy().T
    expected_times = np.array(
        [
            _utc_seconds_to_tdb(query.table["epoch_seconds_UTC"]),
            _utc_seconds_to_tdb(query.table["epoch_seconds_UTC"]),
        ]
    )
    actual_observations = np.array(observation_collection.concatenated_observations).reshape(
        2, -1, order="F"
    )
    actual_times = np.array(observation_collection.concatenated_times).reshape(2, -1, order="F")

    assert np.max(np.abs(actual_observations - expected_observations)) == pytest.approx(0.0)
    assert np.max(np.abs(actual_times - expected_times)) < 1.0e-5


def test_compare_mpc_horizons_eph():
    """Compares true observations from BatchMPC to interpolated simulated RA/DEC from JPL Horizons"""
    batch = BatchMPC()
    batch.get_observations([433])

    batch = _batch_from_augmented_table(
        filter_augmented_optical_table(
            batch.table,
            epoch_start=datetime.datetime(2017, 1, 1),
            epoch_end=datetime.datetime(2022, 1, 1),
            observatories=["T08"],
        )
    )

    # Horizons Query wants batch_times (or start_epoch, end_epoch) in UTC!!!
    batch_times = batch.table.epoch_seconds_UTC.to_numpy()
    eros = HorizonsQuery(
        query_id="433;", location="T08@399", epoch_list=batch_times, extended_query=True
    )

    radec_horizons = eros.interpolated_observations(degrees=False)

    radec_mpc = batch.table.loc[:, ["RA", "DEC"]].reset_index(drop=True)

    diff = radec_horizons[:, 1:3] - radec_mpc.to_numpy()
    diff = np.abs(diff).max(axis=0)
    RA_diff = diff[0]
    DEC_diff = diff[1]

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
