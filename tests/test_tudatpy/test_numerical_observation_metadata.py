"""Numerical metadata follows observation identities rather than transient row indices."""

import pickle
from pathlib import Path

import numpy as np
import pytest

from tudatpy.astro.time_representation import Time
from tudatpy.data_input.tracking_data import ObservationWeightSettings, TrackingData
from tudatpy.data_input.tracking_data.gaia import GaiaAstrometry
from tudatpy.dynamics import environment_setup
from tudatpy.estimation import observations
from tudatpy.estimation.observations import observations_processing


def tracking_data(epochs=(3.0, 1.0, 1.0, 2.0)):
    return TrackingData(
        observable_type="AngularPosition",
        link_ends=[(("673", ""), "transmitter"), (("Gaia", ""), "receiver")],
        observations=[[0.1 + 0.01 * i, 0.2] for i in range(len(epochs))],
        epochs=list(epochs),
        reference_link_end="receiver",
        time_scale="TDB",
    )


def dataset_from_tracking(sources):
    bodies = environment_setup.create_system_of_bodies(
        environment_setup.get_default_body_settings([], "SSB", "J2000")
    )
    return observations.create_observation_dataset_from_tracking_data(sources, bodies)


def test_tracking_metadata_length_removal_and_detached_getters():
    source = tracking_data()
    source.add_numerical_observation_metadata("user supplied key", [30, 10, 11, 20])
    source.add_numerical_observation_metadata("another key", [1, 2, 3, 4])
    with pytest.raises(RuntimeError, match="one entry per observation"):
        source.add_numerical_observation_metadata("user supplied key", [0, 0])
    source.remove_single_observation_entry(1)
    assert source.get_numerical_observation_metadata("user supplied key") == [30, 11, 20]
    assert source.get_numerical_observation_metadata("another key") == [1, 3, 4]
    snapshot = source.numerical_observation_metadata
    snapshot["user supplied key"][0] = -100
    assert source.get_numerical_observation_metadata("user supplied key")[0] == 30


@pytest.mark.parametrize("weighted", [False, True])
def test_conversion_sorts_duplicate_epochs_without_changing_weights(weighted):
    source = tracking_data()
    source.add_numerical_observation_metadata("along_scan_angle", [0.3, 0.1, 0.11, 0.2])
    weight_matrix = np.diag(np.arange(1.0, 9.0)) + 0.01 * np.ones((8, 8))
    if weighted:
        source.set_observation_weight_settings(ObservationWeightSettings.set_block(weight_matrix))
    dataset = dataset_from_tracking([source])
    ids = dataset.observation_ids_for_set(0)
    np.testing.assert_array_equal(
        dataset.get_numerical_observation_metadata("along_scan_angle", ids),
        [0.1, 0.11, 0.2, 0.3],
    )
    np.testing.assert_allclose(
        np.asarray(dataset.observations_for_set(0))[:, 0], [0.11, 0.12, 0.13, 0.1]
    )
    scalar_order = [2, 3, 4, 5, 6, 7, 0, 1]
    expected_weights = weight_matrix[np.ix_(scalar_order, scalar_order)] if weighted else np.eye(8)
    np.testing.assert_array_equal(dataset.get_weight_matrix().toarray(), expected_weights)
    assert source.get_numerical_observation_metadata("along_scan_angle") == [0.3, 0.1, 0.11, 0.2]


def test_native_precision_sorting_does_not_match_metadata_by_float_epoch():
    first = Time(1.0e9) + Time(1.0e-12)
    last = first + Time(2.0e-12)
    assert first < last and float(first) == float(last)
    source = tracking_data([last, first, first])
    source.add_numerical_observation_metadata("angle", [3, 1, 2])
    dataset = dataset_from_tracking([source])
    assert dataset.get_numerical_observation_metadata(
        "angle", dataset.observation_ids_for_set(0)
    ) == [1, 2, 3]


def test_covariance_accessor_includes_rejected_rows_of_full_correlated_set():
    source = tracking_data()
    weights = np.diag(np.arange(1.0, 9.0)) + 0.1 * np.ones((8, 8))
    source.set_observation_weight_settings(ObservationWeightSettings.set_block(weights))
    dataset = dataset_from_tracking([source])
    dataset.reject_observations(observations.observation_query.time == 1)
    snapshot = dataset.observation_vector_data(include_rejected=False)
    scalar_order = [2, 3, 4, 5, 6, 7, 0, 1]
    covariance = np.linalg.inv(weights[np.ix_(scalar_order, scalar_order)])
    for observation_id in dataset.observation_ids_for_set(0):
        start = 2 * observation_id
        np.testing.assert_allclose(
            snapshot.inverse_weight_matrix_for_observation(observation_id),
            covariance[start : start + 2, start : start + 2],
            rtol=1.0e-14,
        )
    detached = snapshot.inverse_weight_matrix_for_observation(2)
    detached[0, 0] = -1
    assert snapshot.inverse_weight_matrix_for_observation(2)[0, 0] > 0


def test_rejection_copy_removal_and_appended_rows_preserve_ids():
    source = tracking_data()
    source.add_numerical_observation_metadata("angle", [30, 10, 11, 20])
    dataset = dataset_from_tracking([source])
    before = dataset.numerical_observation_metadata
    query = observations.observation_query
    dataset.reject_observations(query.time == 1, "test rejection")
    assert dataset.numerical_observation_metadata == before
    subset = dataset.create_new_and_keep(query.time > 1)
    assert subset.numerical_observation_metadata["angle"] == {2: 20, 3: 30}
    subset.add_numerical_observation_metadata("angle", [2], [200])
    assert dataset.numerical_observation_metadata == before
    dataset.restore_observations(query.time == 1)
    dataset.remove_observations(query.time == 1)
    assert dataset.numerical_observation_metadata["angle"] == {2: 20, 3: 30}
    dataset.add_observations_to_set(0, [[0.1, 0.2]], [0.0], sort_observations=True)
    assert dataset.observation_ids_for_set(0) == [4, 2, 3]
    assert dataset.get_numerical_observation_metadata("angle", [3, 2]) == [30, 20]
    with pytest.raises(IndexError):
        dataset.get_numerical_observation_metadata("angle", [4])


def test_user_metadata_updates_validate_ids_and_do_not_invent_missing_values():
    dataset = dataset_from_tracking([tracking_data()])
    dataset.add_numerical_observation_metadata("custom", [0, 2], [0.0, 2.0])
    before = dataset.numerical_observation_metadata
    with pytest.raises(IndexError):
        dataset.add_numerical_observation_metadata("custom", [0, 99], [10, 20])
    with pytest.raises(RuntimeError):
        dataset.add_numerical_observation_metadata("custom", [0, 2], [10])
    assert dataset.numerical_observation_metadata == before
    assert dataset.get_numerical_observation_metadata("custom", [2, 0]) == [2, 0]
    with pytest.raises(IndexError):
        dataset.get_numerical_observation_metadata("custom", [1])
    with pytest.raises(IndexError):
        dataset.get_numerical_observation_metadata("missing key", [0])
    snapshot = dataset.numerical_observation_metadata
    snapshot["custom"][0] = -1
    assert dataset.numerical_observation_metadata == before


def test_copying_sets_remaps_metadata_to_new_ids_including_within_same_dataset():
    source = tracking_data()
    source.add_numerical_observation_metadata("angle", [30, 10, 11, 20])
    dataset = dataset_from_tracking([source])
    target = dataset_from_tracking([tracking_data([0.0])])
    new_set = target.add_observation_set_from_dataset(dataset, 0)
    assert target.observation_ids_for_set(new_set) == [1, 2, 3, 4]
    assert target.get_numerical_observation_metadata("angle", [1, 2, 3, 4]) == [10, 11, 20, 30]
    same_dataset_set = dataset.add_observation_set_from_dataset(dataset, 0)
    assert dataset.get_numerical_observation_metadata(
        "angle", dataset.observation_ids_for_set(same_dataset_set)
    ) == [10, 11, 20, 30]


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
def test_transfer_to_filtered_set_and_serialization_preserve_numerical_metadata():
    source = tracking_data([1.0, 2.0, 3.0])
    source.add_numerical_observation_metadata("angle", [10, 20, 30])
    dataset = dataset_from_tracking([source])
    observation_set = observations.create_single_observation_set_from_dataset(dataset, 0)
    observation_set.filter_observations(
        observations_processing.observation_filter(
            observations_processing.ObservationFilterType.time_bounds_filtering,
            1.5,
            2.5,
            use_opposite_condition=True,
        )
    )
    assert dataset.numerical_observation_metadata["angle"] == {1: 20}
    restored = pickle.loads(pickle.dumps(observation_set))
    filtered = observations.create_observation_dataset_from_single_observation_set(
        restored.filtered_observation_set
    )
    assert filtered.get_numerical_observation_metadata(
        "angle", filtered.observation_ids_for_set(0)
    ) == [10, 30]


def test_gaia_angles_match_ccd_rows_after_filtering_and_dataset_conversion():
    archive = Path(__file__).with_name("gaia_astro_archive_for_tests.parquet")
    gaia = GaiaAstrometry.load_from_local_archive(archive, 673)
    gaia.apply_filters()
    sources, _ = gaia.to_tracking_data()
    table = gaia.table
    expected = []
    for source, (_, rows) in zip(sources, table.groupby("transit_id", sort=False)):
        angles = rows["position_angle_scan"].tolist()
        assert source.get_numerical_observation_metadata("along_scan_angle") == angles
        expected.append(angles)
    dataset = dataset_from_tracking(sources)
    snapshot = dataset.observation_vector_data()
    for set_id, (angles, covariance) in enumerate(
        zip(expected, gaia._get_observation_covariance(673))
    ):
        ids = dataset.observation_ids_for_set(set_id)
        assert dataset.get_numerical_observation_metadata("along_scan_angle", ids) == angles
        for index, observation_id in enumerate(ids):
            start = 2 * index
            np.testing.assert_allclose(
                snapshot.inverse_weight_matrix_for_observation(observation_id),
                covariance[start : start + 2, start : start + 2],
                rtol=1.0e-6,
                atol=1.0e-25,
            )
