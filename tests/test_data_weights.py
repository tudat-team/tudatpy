
# tests for data weights functionality
from tudatpy.numerical_simulation import environment_setup
from tudatpy.numerical_simulation import estimation
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
@pytest.mark.parametrize("use_dummy_weights", [(True,), (False,)])
def test_MPC_weights_to_ObsCol(
    observatories_to_filter, use_dummy_weights, use_single_observation
):
    """Test if the weights are transfered correctly to observation collection"""
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

    observation_collection = batch.to_tudat(
        bodies=bodies,
        included_satellites=None,
        apply_star_catalog_debias=True,
        apply_weights_VFCC17=True,
    )

    # tudat's observationcollection sorts by observatory then time
    temp_table = batch._table.sort_values(
        ["observatory", "epochJ2000secondsTDB"], ascending=True
    )
    # concatted weights goes [RA1, DEC1, RA2, DEC2, ...]
    batch_weights = np.ravel(2 * [temp_table.weight.to_numpy()], "F")
    batch_times = np.ravel(2 * [temp_table.epochJ2000secondsTDB.to_numpy()], "F")

    # check if lengths match and if the difference is zero
    print(len(batch_weights), len(observation_collection.concatenated_weights))
    assert len(batch_weights) == len(observation_collection.concatenated_weights)
    total_diff = np.sum(
        batch_weights - np.array(observation_collection.concatenated_weights)
    )
    total_diff_time = np.sum(
        batch_times - np.array(observation_collection.concatenated_times)
    )

    assert total_diff_time == 0
    assert total_diff == 0

    # test pod_input
    # provide the observation collection as input, and limit number of iterations for estimation.
    pod_input = estimation.EstimationInput(
        observations_and_times=observation_collection,
        convergence_checker=estimation.estimation_convergence_checker(
            maximum_iterations=1,
        ),
    )
    pod_input.set_weights_from_observation_collection()

