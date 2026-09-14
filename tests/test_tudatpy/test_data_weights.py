# tests for data weights functionality
from tudatpy.dynamics import environment_setup
from tudatpy.dynamics.environment_setup import ground_station
from tudatpy.estimation.observations import create_observation_collection_from_tracking_data
from tudatpy.data_input.tracking_data.mpc import BatchMPC
from tudatpy.data_input.tracking_data.optical_utilities import (
    create_augmented_optical_table,
    filter_augmented_optical_table,
)
from tudatpy.astro import time_representation
from tudatpy.data_input.environment_data import spice

import numpy as np
import datetime

import pytest

spice.load_standard_kernels()
_TIME_SCALE_CONVERTER = time_representation.default_time_scale_converter()

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
    batch._table = create_augmented_optical_table(table, in_degrees=False)
    batch._refresh_metadata()
    return batch


@pytest.mark.parametrize(
    "observatories_to_filter,use_single_observation", weights_test_combinations
)
@pytest.mark.parametrize("use_dummy_weights", [True, False])
@pytest.mark.remote_data
def test_MPC_weights_to_ObsCol(observatories_to_filter, use_dummy_weights, use_single_observation):
    """Test if the weights are transfered correctly to observation collection"""
    if use_dummy_weights:
        pytest.skip("Custom per-observation MPC weights are not part of the current BatchMPC API.")

    target_mpc_code = "433"
    mpc_codes = [target_mpc_code]

    observations_start = datetime.datetime(2023, 1, 1)
    observations_end = datetime.datetime(2024, 1, 1)
    global_frame_origin = "SSB"
    global_frame_orientation = "J2000"

    # Create system of bodies
    bodies_to_create = ["Sun", "Earth"]
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation
    )
    body_settings.get("Earth").ground_station_settings = ground_station.optical_telescope_stations()
    bodies = environment_setup.create_system_of_bodies(body_settings)

    batch = BatchMPC()
    batch.get_observations(mpc_codes)
    batch = _batch_from_augmented_table(
        filter_augmented_optical_table(
            batch.table,
            epoch_start=observations_start,
            epoch_end=observations_end,
            observatories=observatories_to_filter,
            observatories_exclude=(["C51"] if observatories_to_filter is None else None),
        )
    )

    if use_single_observation:
        # gets the first item and remakes the batch from this 1 item dataframe
        single_observation_table = batch.table.iloc[0:1].copy()
        batch = _batch_from_augmented_table(single_observation_table)

    tracking_data, supplementary_data = batch.to_tracking_dataset(
        add_star_catalog_corrections=True,
        add_weights=True,
    )
    assert supplementary_data == []
    assert all(data.weighing_scheme == "VFCC17" for data in tracking_data)
    observation_collection = create_observation_collection_from_tracking_data(tracking_data, bodies)

    # tudat's observationcollection sorts by observatory then time
    temp_table = batch._table.sort_values(["observatory", "epoch_seconds_UTC"], ascending=True)

    # concatted values go [RA1, DEC1, RA2, DEC2, ...]
    batch_times = np.ravel(2 * [_utc_seconds_to_tdb(temp_table.epoch_seconds_UTC)], "F")

    collection_weights = np.array(observation_collection.concatenated_weights)
    time_difference = batch_times - np.array(observation_collection.concatenated_times)

    assert len(collection_weights) == 2 * len(temp_table)
    assert np.all(np.isfinite(collection_weights))
    assert np.all(collection_weights > 0.0)
    assert np.max(np.abs(time_difference)) < 1.0e-5
