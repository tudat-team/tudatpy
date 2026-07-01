# tests for data weights functionality
from tudatpy.dynamics import environment_setup
from tudatpy.data.mpc import BatchMPC
from tudatpy.interface import spice

import numpy as np
import datetime

import pytest

spice.load_standard_kernels()

observatory_set_single = ["M22"]
observatory_set_multi = ["K19", "D67", "089", "706"]
weights_test_combinations = [
    (observatory_set_single, True),  # just one obs
    (observatory_set_single, False),
    (observatory_set_multi, False),
    (None, False),  # all data
]


@pytest.mark.parametrize(
    "observatories_to_filter,use_single_observation", weights_test_combinations
)
@pytest.mark.parametrize("use_dummy_weights", [True, False])
def test_MPC_weights_to_observation_dataset(
    observatories_to_filter, use_dummy_weights, use_single_observation
):
    """Test if the weights are transferred correctly to the observation dataset."""
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
    bodies = environment_setup.create_system_of_bodies(body_settings)

    batch = BatchMPC()
    batch.get_observations(mpc_codes)
    batch.filter(
        epoch_start=observations_start,
        epoch_end=observations_end,
        observatories=observatories_to_filter,
        observatories_exclude=(["C51"] if observatories_to_filter is None else None),
    )

    if use_single_observation:
        # gets the first item and remakes the batch from this 1 item dataframe
        batch.from_pandas(batch.table.iloc[0:1])

    if use_dummy_weights:
        # sets the weights to be a list in ascending order from 0, 1, 2,...
        batch.set_weights(np.array(list(range(0, batch.size))))

    observation_dataset = batch.to_tudat(
        bodies=bodies,
        included_satellites=None,
        apply_star_catalog_debias=True,
        apply_weights_VFCC17=True,
    )

    # Tudat's ordered flattened dataset data sorts by observatory then time.
    temp_table = batch._table.query("observatory != @batch.space_telescopes").sort_values(
        ["observatory", "epoch_seconds_TDB"], ascending=True
    )

    # concatted weights goes [RA1, DEC1, RA2, DEC2, ...]
    batch_weights = np.ravel(2 * [temp_table.weight.to_numpy()], "F")
    batch_times = np.ravel(2 * [temp_table.epoch_seconds_TDB.to_numpy()], "F")

    flattened_data = observation_dataset.ordered_flattened_observation_data()

    # The flattened dataset must preserve every expected weight in ordered-output order.
    assert len(batch_weights) == len(flattened_data.weight_vector)
    np.testing.assert_allclose(np.array(flattened_data.weight_vector), batch_weights)

    # Times are checked entry-by-entry so an ordering regression cannot be hidden by cancellation.
    np.testing.assert_allclose(np.array(flattened_data.times), batch_times)
