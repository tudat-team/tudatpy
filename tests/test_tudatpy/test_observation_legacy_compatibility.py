import warnings

import numpy as np
import pytest

from tudatpy.dynamics import environment_setup
from tudatpy.estimation import observations
from tudatpy.estimation.observations import observations_processing


def _link_ends(receiver_body):
    return {
        observations.transmitter: observations.LinkEndId("Probe", ""),
        observations.receiver: observations.LinkEndId(receiver_body, ""),
    }


def _new_dataset_single_set(observable_type, receiver_body, observation_values, times):
    dataset = observations.ObservationDataset()
    dataset.add_observation_set(
        observable_type,
        observations.LinkDefinition(_link_ends(receiver_body)),
        [np.asarray(value, dtype=float) for value in observation_values],
        times,
        observations.receiver,
    )
    return dataset


@pytest.fixture
def sample_dataset():
    dataset = _new_dataset_single_set(
        observations.one_way_range,
        "Earth",
        [[10.0], [20.0], [30.0]],
        [1.0, 2.0, 4.0],
    )
    dataset.set_residuals_for_set(
        0,
        [np.array([0.1]), np.array([5.0]), np.array([-0.2])],
    )

    angular_dataset = _new_dataset_single_set(
        observations.angular_position,
        "Mars",
        [[1.0, 2.0], [3.0, 4.0]],
        [3.0, 5.0],
    )
    assert dataset.add_observation_set_from_dataset(angular_dataset, 0) == 1

    dataset.set_residuals_for_set(
        1,
        [np.array([0.01, 0.02]), np.array([10.0, 0.04])],
    )
    dataset.set_weight_vector_for_set(1, np.array([2.0, 3.0, 4.0, 5.0]))

    return dataset


def _legacy_collection(dataset):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return observations.create_observation_collection_from_dataset(dataset)


def _link_definition_signature(link_definition):
    signature = []
    for link_end_type, link_end_id in link_definition.link_ends.items():
        signature.append(
            (
                str(link_end_type),
                link_end_id.body_name,
                link_end_id.reference_point,
            )
        )
    return tuple(sorted(signature))


def _link_definition_map_signature(link_definition_map):
    return {
        key: _link_definition_signature(observations.LinkDefinition(link_ends))
        for key, link_ends in link_definition_map.items()
    }


def _link_ends_per_observable_signature(link_ends_per_observable_type):
    return {
        observable_type: [
            _link_definition_signature(observations.LinkDefinition(link_ends))
            for link_ends in link_ends_list
        ]
        for observable_type, link_ends_list in link_ends_per_observable_type.items()
    }


def _to_dense_matrix(matrix):
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


def test_legacy_single_observation_set_conversion_matches_dataset():
    """Verify legacy single-set conversion preserves dataset values and metadata."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        observation_set = observations.create_single_observation_set(
            observations.angular_position,
            _link_ends("Mars"),
            [np.array([1.0, 2.0]), np.array([3.0, 4.0])],
            [3.0, 5.0],
            observations.receiver,
        )
        dataset = observations.create_observation_dataset_from_single_observation_set(
            observation_set
        )

    np.testing.assert_allclose(dataset.observation_vector_for_set(0), [1.0, 2.0, 3.0, 4.0])
    np.testing.assert_allclose(
        [float(time) for time in dataset.observation_times_for_set(0)],
        [3.0, 5.0],
    )


def test_legacy_saved_filtered_set_tracks_link_definition_updates():
    """Verify repeated legacy filtering after a link-definition update."""
    with warnings.catch_warnings(record=True) as warning_records:
        warnings.simplefilter("always")
        observation_set = observations.create_single_observation_set(
            observations.one_way_range,
            _link_ends("Earth"),
            [np.array([10.0]), np.array([20.0]), np.array([30.0])],
            [1.0, 2.0, 4.0],
            observations.receiver,
        )

        # Save observations outside the time window in the filtered set.
        time_filter = observations_processing.observation_filter(
            observations_processing.ObservationFilterType.time_bounds_filtering,
            1.5,
            3.5,
            use_opposite_condition=True,
        )
        observation_set.filter_observations(time_filter)

        # Update the active and saved observations before filtering again.
        updated_link_ends = _link_ends("Earth")
        updated_link_ends[observations.receiver] = observations.LinkEndId("Earth", "Station")
        observation_set.link_definition = observations.LinkDefinition(updated_link_ends)
        observation_set.set_residuals([np.array([5.0])])
        residual_filter = observations_processing.observation_filter(
            observations_processing.ObservationFilterType.residual_filtering,
            0.1,
        )
        observation_set.filter_observations(residual_filter)
        active_observation_count = observation_set.number_of_observables
        filtered_observation_count = observation_set.filtered_observation_set.number_of_observables

    # Both legacy operations remain usable and explain their modern replacements.
    messages = [str(record.message) for record in warning_records]
    assert any("create_single_observation_set is deprecated" in message for message in messages)
    filter_messages = [
        message for message in messages if "SingleObservationSet.filter_observations" in message
    ]
    assert len(filter_messages) >= 2
    assert all("ObservationDataset" in message for message in filter_messages)
    assert active_observation_count == 0
    assert filtered_observation_count == 3


def test_legacy_reference_point_keeps_live_dataset_and_filtered_metadata():
    """Verify reference-point changes preserve live data and repeated filtering."""
    dataset = _new_dataset_single_set(
        observations.one_way_range,
        "Probe",
        [[10.0], [20.0], [30.0]],
        [1.0, 2.0, 4.0],
    )
    dataset.add_observation_set(
        observations.one_way_range,
        observations.LinkDefinition(_link_ends("Mars")),
        [],
        [],
        observations.receiver,
    )
    collection = _legacy_collection(dataset)
    bodies = environment_setup.create_system_of_bodies(
        environment_setup.BodyListSettings("SSB", "J2000")
    )
    bodies.create_empty_body("Probe")

    with warnings.catch_warnings(record=True) as warning_records:
        warnings.simplefilter("always")

        # Save two rows, then change the reference point of both active and saved data.
        time_filter = observations_processing.observation_filter(
            observations_processing.ObservationFilterType.time_bounds_filtering,
            1.5,
            3.5,
            use_opposite_condition=True,
        )
        collection.filter_observations(time_filter)
        collection.remove_empty_observation_sets()
        collection.set_reference_point(
            bodies,
            np.array([1.0, 2.0, 3.0]),
            "Antenna",
            "Probe",
            observations.receiver,
        )

        # Mutations through the returned dataset must remain visible to the collection.
        live_dataset = observations.create_observation_dataset_from_collection(collection)
        live_dataset.set_residuals_for_set(0, [np.array([5.0])])
        residuals = np.asarray(collection.get_concatenated_residuals()).reshape(-1)

        # A second filter must move the active row into the compatible saved set.
        residual_filter = observations_processing.observation_filter(
            observations_processing.ObservationFilterType.residual_filtering,
            0.1,
        )
        collection.filter_observations(residual_filter)
        observation_set = collection.get_single_observation_sets()[0]
        active_observation_count = observation_set.number_of_observables
        filtered_observation_count = observation_set.filtered_observation_set.number_of_observables

    messages = [str(record.message) for record in warning_records]
    assert any("set_reference_point is deprecated" in message for message in messages)
    assert any(
        "create_observation_dataset_from_collection is deprecated" in message
        for message in messages
    )
    np.testing.assert_allclose(residuals, [5.0])
    assert active_observation_count == 0
    assert filtered_observation_count == 3


def test_legacy_reference_point_switch_history_ignores_empty_sets():
    """Verify reference-point switch histories tolerate empty observation sets."""
    dataset = observations.ObservationDataset()
    dataset.add_observation_set(
        observations.one_way_range,
        observations.LinkDefinition(_link_ends("Earth")),
        [],
        [],
        observations.receiver,
    )
    dataset.add_observation_set(
        observations.one_way_range,
        observations.LinkDefinition(_link_ends("Mars")),
        [np.array([10.0]), np.array([20.0]), np.array([30.0])],
        [1.0, 2.0, 4.0],
        observations.receiver,
    )
    collection = _legacy_collection(dataset)
    bodies = environment_setup.create_system_of_bodies(
        environment_setup.BodyListSettings("SSB", "J2000")
    )
    bodies.create_empty_body("Probe")

    with warnings.catch_warnings(record=True) as warning_records:
        warnings.simplefilter("always")

        # Split the populated set at the antenna switch while leaving the empty set alone.
        collection.set_reference_points(
            bodies,
            {
                0.0: np.array([1.0, 0.0, 0.0]),
                3.0: np.array([2.0, 0.0, 0.0]),
                5.0: np.array([2.0, 0.0, 0.0]),
            },
            "Probe",
            observations.transmitter,
        )
        observation_sets = collection.get_single_observation_sets()
        populated_sets = [
            observation_set
            for observation_set in observation_sets
            if observation_set.number_of_observables > 0
        ]
        populated_reference_points = {
            observation_set.link_definition.link_ends[observations.transmitter].reference_point
            for observation_set in populated_sets
        }

        # Residual updates through the converted dataset must remain live.
        live_dataset = observations.create_observation_dataset_from_collection(collection)
        live_dataset.set_residuals_for_set(0, [np.array([5.0]), np.array([5.0])])
        live_dataset.set_residuals_for_set(1, [np.array([5.0])])
        residuals = np.asarray(collection.get_concatenated_residuals()).reshape(-1)

    messages = [str(record.message) for record in warning_records]
    assert any("set_reference_points is deprecated" in message for message in messages)
    assert len(observation_sets) == 2
    assert len(populated_sets) == 2
    assert populated_reference_points == {"Antenna1", "Antenna2"}
    np.testing.assert_allclose(residuals, [5.0, 5.0, 5.0])


def test_legacy_ephemeris_reference_point_adopts_grouped_collection_dataset():
    """Verify an ephemeris reference point keeps grouped collections live."""
    with warnings.catch_warnings(record=True) as warning_records:
        warnings.simplefilter("always")
        observation_set = observations.create_single_observation_set(
            observations.one_way_range,
            _link_ends("Earth"),
            [np.array([10.0]), np.array([20.0])],
            [1.0, 2.0],
            observations.receiver,
        )
        collection = observations.ObservationCollection([observation_set])
        bodies = environment_setup.create_system_of_bodies(
            environment_setup.BodyListSettings("SSB", "J2000")
        )
        bodies.create_empty_body("Probe")
        antenna_ephemeris = environment_setup.ephemeris.create_ephemeris(
            environment_setup.ephemeris.constant(
                np.zeros(6),
                "SSB",
                "J2000",
            ),
            "Antenna",
        )

        # This overload is used by the older MEX residual example.
        collection.set_reference_point(
            bodies,
            antenna_ephemeris,
            "Antenna",
            "Probe",
            observations.transmitter,
        )
        live_dataset = observations.create_observation_dataset_from_collection(collection)
        live_dataset.set_residuals_for_set(0, [np.array([5.0]), np.array([6.0])])
        residuals = np.asarray(collection.get_concatenated_residuals()).reshape(-1)

    messages = [str(record.message) for record in warning_records]
    assert any("set_reference_point is deprecated" in message for message in messages)
    assert any(
        "create_observation_dataset_from_collection is deprecated" in message
        for message in messages
    )
    np.testing.assert_allclose(residuals, [5.0, 6.0])


def test_dataset_add_observation_set_accepts_single_observation_set_object():
    """Verify a legacy single-set object can be added directly to a dataset."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        observation_set = observations.create_single_observation_set(
            observations.angular_position,
            _link_ends("Mars"),
            [np.array([1.0, 2.0]), np.array([3.0, 4.0])],
            [3.0, 5.0],
            observations.receiver,
        )

    dataset = observations.ObservationDataset()
    assert dataset.add_observation_set(observation_set) == 0

    np.testing.assert_allclose(dataset.observation_vector_for_set(0), [1.0, 2.0, 3.0, 4.0])
    np.testing.assert_allclose(
        [float(time) for time in dataset.observation_times_for_set(0)],
        [3.0, 5.0],
    )


def test_dataset_add_observation_set_component_shape_matches_keyword_construction():
    """Verify positional and keyword dataset construction preserve component shape."""
    link_definition = observations.LinkDefinition(_link_ends("Earth"))
    observation_values = [np.array([20.0]), np.array([10.0])]
    observation_times = [2.0, 1.0]

    positional_dataset = observations.ObservationDataset()
    assert (
        positional_dataset.add_observation_set(
            observations.one_way_range,
            link_definition,
            observation_values,
            observation_times,
            observations.receiver,
        )
        == 0
    )

    keyword_dataset = observations.ObservationDataset()
    assert (
        keyword_dataset.add_observation_set(
            observable_type=observations.one_way_range,
            link_definition=link_definition,
            observation_values=observation_values,
            observation_times=observation_times,
            reference_link_end=observations.receiver,
        )
        == 0
    )

    np.testing.assert_allclose(
        positional_dataset.ordered_observation_vector_data().observation_vector,
        keyword_dataset.ordered_observation_vector_data().observation_vector,
    )
    np.testing.assert_allclose(
        positional_dataset.ordered_observation_vector_data().times,
        keyword_dataset.ordered_observation_vector_data().times,
    )
    np.testing.assert_allclose(
        [float(time) for time in positional_dataset.observation_times_for_set(0)],
        observation_times,
    )


def test_dataset_add_observation_set_wrong_shape_names_accepted_signatures():
    """Verify invalid add-observation-set calls report both supported signatures."""
    with pytest.raises(TypeError) as exception_info:
        observations.ObservationDataset().add_observation_set(object())

    error_message = str(exception_info.value)
    assert "single_observation_set" in error_message
    assert "observable_type" in error_message
    assert "link_definition" in error_message


def test_dataset_exposes_legacy_collection_vector_properties(sample_dataset):
    """Verify dataset legacy vector properties match the collection facade."""
    legacy_collection = _legacy_collection(sample_dataset)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.testing.assert_allclose(
            sample_dataset.concatenated_times,
            legacy_collection.concatenated_times,
        )
        np.testing.assert_allclose(
            sample_dataset.concatenated_observations,
            legacy_collection.concatenated_observations,
        )
        np.testing.assert_allclose(
            sample_dataset.concatenated_weights,
            legacy_collection.concatenated_weights,
        )
        assert sample_dataset.concatenated_link_definition_ids == (
            legacy_collection.concatenated_link_definition_ids
        )
        assert sample_dataset.observation_vector_size == sample_dataset.total_scalar_size
        assert sample_dataset.observation_vector_size == (legacy_collection.observation_vector_size)
        assert sample_dataset.time_bounds == legacy_collection.time_bounds


def test_dataset_exposes_legacy_collection_metadata_properties(sample_dataset):
    """Verify dataset legacy metadata properties match the collection facade."""
    legacy_collection = _legacy_collection(sample_dataset)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert sample_dataset.observable_type_start_index_and_size == (
            legacy_collection.observable_type_start_index_and_size
        )
        assert sample_dataset.observation_set_start_index_and_size == (
            legacy_collection.observation_set_start_index_and_size
        )
        assert _link_definition_map_signature(sample_dataset.link_definition_ids) == (
            _link_definition_map_signature(legacy_collection.link_definition_ids)
        )
        assert _link_ends_per_observable_signature(
            sample_dataset.link_ends_per_observable_type
        ) == _link_ends_per_observable_signature(legacy_collection.link_ends_per_observable_type)
        assert {
            observable_type: [
                _link_definition_signature(link_definition) for link_definition in link_definitions
            ]
            for observable_type, link_definitions in (
                sample_dataset.link_definitions_per_observable.items()
            )
        } == {
            observable_type: [
                _link_definition_signature(link_definition) for link_definition in link_definitions
            ]
            for observable_type, link_definitions in (
                legacy_collection.link_definitions_per_observable.items()
            )
        }
        assert sample_dataset.sorted_per_set_time_bounds == (
            legacy_collection.sorted_per_set_time_bounds
        )
        assert sample_dataset.sorted_observation_sets.keys() == (
            legacy_collection.sorted_observation_sets.keys()
        )


def test_dataset_exposes_legacy_collection_lookup_methods(sample_dataset):
    """Verify dataset legacy lookup methods match the collection facade."""
    legacy_collection = _legacy_collection(sample_dataset)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert sample_dataset.get_observable_types() == legacy_collection.get_observable_types()
        assert sample_dataset.get_bodies_in_link_ends() == (
            legacy_collection.get_bodies_in_link_ends()
        )
        assert [
            _link_definition_signature(link_definition)
            for link_definition in sample_dataset.get_link_definitions_for_observables(
                observations.angular_position
            )
        ] == [
            _link_definition_signature(link_definition)
            for link_definition in legacy_collection.get_link_definitions_for_observables(
                observations.angular_position
            )
        ]

        angular_link = sample_dataset.get_link_definitions_for_observables(
            observations.angular_position
        )[0]
        dataset_single_link = sample_dataset.get_single_link_and_type_observations(
            observations.angular_position,
            angular_link,
        )
        legacy_single_link = legacy_collection.get_single_link_and_type_observations(
            observations.angular_position,
            angular_link,
        )
        assert len(dataset_single_link) == len(legacy_single_link)
        np.testing.assert_allclose(
            dataset_single_link[0].concatenated_observations,
            legacy_single_link[0].concatenated_observations,
        )


def test_dataset_delegates_remaining_legacy_collection_method_names(sample_dataset):
    """Verify retained legacy method names delegate to dataset behavior."""
    legacy_collection = _legacy_collection(sample_dataset)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.testing.assert_allclose(
            sample_dataset.get_concatenated_observations(),
            legacy_collection.get_concatenated_observations(),
        )
        np.testing.assert_allclose(
            sample_dataset.get_concatenated_weights(),
            legacy_collection.get_concatenated_weights(),
        )
        assert sample_dataset.get_time_bounds_list() == legacy_collection.get_time_bounds_list()

        sample_dataset.set_constant_weight(12.0)
        np.testing.assert_allclose(
            sample_dataset.ordered_observation_vector_data().weight_vector,
            np.full(sample_dataset.total_scalar_size, 12.0),
        )


def test_dataset_legacy_full_vector_setters_update_dataset(sample_dataset):
    """Verify legacy full-vector setters update the dataset backend."""
    vector_data = sample_dataset.ordered_observation_vector_data()
    new_observations = np.asarray(vector_data.observation_vector, dtype=float) + 100.0
    new_residuals = np.asarray(vector_data.residual_vector, dtype=float) - 0.25

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sample_dataset.set_observations(new_observations)
        sample_dataset.set_residuals(new_residuals)

    updated = sample_dataset.ordered_observation_vector_data()
    np.testing.assert_allclose(updated.observation_vector, new_observations)
    np.testing.assert_allclose(updated.residual_vector, new_residuals)


def test_legacy_weight_setters_match_dataset(sample_dataset):
    """Verify legacy collection weight setters match dataset weights."""
    legacy_collection = _legacy_collection(sample_dataset)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        range_parser = observations_processing.observation_parser(observations.one_way_range)
        angular_parser = observations_processing.observation_parser(observations.angular_position)

    def assert_weights_match(expected_weights):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            np.testing.assert_allclose(
                np.asarray(legacy_collection.concatenated_weights),
                expected_weights,
            )
        vector_data = sample_dataset.ordered_observation_vector_data()
        np.testing.assert_allclose(vector_data.weight_vector, expected_weights)
        np.testing.assert_allclose(
            _to_dense_matrix(vector_data.sparse_weight_matrix),
            np.diag(expected_weights),
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            if hasattr(legacy_collection, "get_full_weights_matrix"):
                np.testing.assert_allclose(
                    np.asarray(legacy_collection.get_full_weights_matrix()),
                    np.diag(expected_weights),
                )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        legacy_collection.set_constant_weight_per_observation_parser(
            {range_parser: 4.0, angular_parser: 9.0}
        )
    assert_weights_match(np.array([4.0, 4.0, 4.0, 9.0, 9.0, 9.0, 9.0]))

    tabulated_weights = {
        range_parser: np.array([1.1, 1.2, 1.3]),
        angular_parser: np.array([2.1, 2.2, 2.3, 2.4]),
    }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        legacy_collection.set_tabulated_weights(tabulated_weights)
    assert_weights_match(np.array([1.1, 1.2, 1.3, 2.1, 2.2, 2.3, 2.4]))

    source_snapshot = sample_dataset.get_observations(observations.observation_query.active)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        legacy_collection.remove_single_observation_sets(range_parser)
    # Removing a collection member changes its grouping, not its dataset owner.
    assert len(source_snapshot) == 5
    assert sample_dataset.number_of_observation_sets == 2
    np.testing.assert_allclose(
        sample_dataset.ordered_observation_vector_data().weight_vector,
        [1.1, 1.2, 1.3, 2.1, 2.2, 2.3, 2.4],
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.testing.assert_allclose(
            np.asarray(legacy_collection.concatenated_weights),
            [2.1, 2.2, 2.3, 2.4],
        )


@pytest.mark.parametrize("container", [list, lambda values: np.array(values, dtype=object)])
def test_dataset_preserves_precise_time_objects(container):
    """Verify dataset input preserves precise time objects in common containers."""
    from tudatpy.astro.time_representation import Time

    epoch = Time(1000000, 0.000000001)
    dataset = _new_dataset_single_set(
        observations.one_way_range, "Earth", [[10.0]], container([epoch])
    )
    stored = dataset.observation_times_for_set(0)[0]
    assert stored == epoch
    assert float(stored - Time(1000000, 0.0)) == pytest.approx(1.0e-9, rel=1.0e-12)


def test_scalar_setters_accept_one_dimensional_arrays(sample_dataset):
    """Verify scalar observation and residual setters accept one-dimensional arrays."""
    sample_dataset.set_observations_for_set(0, np.array([101.0, 102.0, 103.0]))
    sample_dataset.set_residuals_for_set(0, np.array([0.5, 1.5, 2.5]))
    np.testing.assert_array_equal(
        sample_dataset.observation_vector_for_set(0), [101.0, 102.0, 103.0]
    )
    np.testing.assert_array_equal(sample_dataset.residual_vector_for_set(0), [0.5, 1.5, 2.5])
