import numpy as np
import pytest

from tudatpy.dynamics import environment_setup
from tudatpy.estimation import estimation_analysis, observations
from tudatpy.estimation.observable_models_setup import model_settings
from tudatpy.estimation.observations_setup import observations_dependent_variables as dependent
from tudatpy.estimation.observations_setup import observations_simulation_settings as simulation


def test_legacy_residual_computation_updates_original_sets():
    body_settings = environment_setup.BodyListSettings("SSB", "J2000")
    for name, x in [("Probe", 1000.0), ("Earth", 0.0)]:
        body_settings.add_empty_settings(name)
        body_settings.get(name).ephemeris_settings = environment_setup.ephemeris.constant(
            np.array([x, 0, 0, 0, 0, 0]), "SSB", "J2000"
        )
    bodies = environment_setup.create_system_of_bodies(body_settings)
    definition = observations.LinkDefinition(_link_ends("Earth"))
    settings = simulation.tabulated_simulation_settings(
        observations.one_way_range, definition, [1.0, 2.0]
    )
    settings.add_dependent_variables(
        [
            dependent.target_range_between_link_ends_dependent_variable(
                observations.transmitter, observations.receiver
            )
        ]
    )
    sets = [
        observations.SingleObservationSet(
            observations.one_way_range,
            definition,
            [[1005.0 + i]],
            [1.0 + i],
            observations.receiver,
            dependent_variable_bookkeeping=settings.dependent_variable_bookkeeping,
        )
        for i in range(2)
    ]
    collection = observations.ObservationCollection(sets)
    original_sets = collection.get_single_observation_sets()
    simulators = simulation.create_observation_simulators(
        [model_settings.one_way_range(definition)], bodies
    )
    for expected_residuals in ([5.0, 6.0], [10.0, 20.0]):
        collection.set_observations(1000.0 + np.array(expected_residuals))
        observations.compute_residuals_and_dependent_variables(collection, simulators, bodies)
        np.testing.assert_allclose(
            collection.get_concatenated_residuals(), expected_residuals, atol=1e-9
        )
        for original, current, residual in zip(
            original_sets, collection.get_single_observation_sets(), expected_residuals
        ):
            assert current is original
            np.testing.assert_allclose(
                np.asarray(original.residuals).ravel(), [residual], atol=1e-9
            )
            np.testing.assert_allclose(np.asarray(original.dependent_variables).ravel(), [1000.0])


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
    return observations.create_observation_collection_from_dataset(dataset)


@pytest.mark.parametrize(
    "input_class",
    [estimation_analysis.CovarianceAnalysisInput, estimation_analysis.EstimationInput],
)
@pytest.mark.parametrize("source_type", ["dataset", "collection"])
def test_retained_input_weight_alias_updates_source_and_diagonal_is_readonly(
    sample_dataset, input_class, source_type
):
    source = sample_dataset if source_type == "dataset" else _legacy_collection(sample_dataset)
    input_settings = input_class(source)
    input_settings.set_constant_weight(7.0)
    if source_type == "dataset":
        np.testing.assert_allclose(source.get_weight_diagonal(), 7.0)
    else:
        np.testing.assert_allclose(source.concatenated_weights, 7.0)
    with pytest.raises(AttributeError):
        input_settings.weight_matrix_diagonal = np.ones(sample_dataset.total_scalar_size)


def test_older_deprecated_observation_and_input_apis_remain_removed():
    assert not hasattr(observations, "single_observation_set")
    assert not hasattr(observations.SingleObservationSet, "ancilliary_settings")
    for input_class in (
        estimation_analysis.CovarianceAnalysisInput,
        estimation_analysis.EstimationInput,
    ):
        for name in (
            "set_weights_from_observation_collection",
            "set_constant_single_observable_weight",
            "set_constant_single_observable_vector_weight",
            "set_constant_single_observable_and_link_end_weight",
            "set_constant_single_observable_and_link_end_vector_weight",
            "set_total_single_observable_and_link_end_vector_weight",
            "set_constant_weight_per_observable",
            "set_constant_vector_weight_per_observable",
        ):
            assert not hasattr(input_class, name)


def test_legacy_single_observation_set_conversion_matches_dataset():
    """Verify legacy single-set conversion preserves dataset values and metadata."""
    observation_set = observations.create_single_observation_set(
        observations.angular_position,
        _link_ends("Mars"),
        [np.array([1.0, 2.0]), np.array([3.0, 4.0])],
        [3.0, 5.0],
        observations.receiver,
    )
    dataset = observations.create_observation_dataset_from_single_observation_set(observation_set)

    np.testing.assert_allclose(dataset.observation_vector_for_set(0), [1.0, 2.0, 3.0, 4.0])
    np.testing.assert_allclose(
        [float(time) for time in dataset.observation_times_for_set(0)],
        [3.0, 5.0],
    )


def test_dataset_add_observation_set_accepts_single_observation_set_object():
    """Verify a legacy single-set object can be added directly to a dataset."""
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
            observations=observation_values,
            times=observation_times,
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


@pytest.mark.parametrize("positional", [False, True])
def test_estimation_input_combines_outlier_and_apriori_settings(sample_dataset, positional):
    """The dataset interface retains the seventh outlier argument when enabling a priori deviations."""
    outliers = estimation_analysis.simple_outlier_rejection_settings(10.0)
    if positional:
        estimation_input = estimation_analysis.EstimationInput(
            sample_dataset,
            np.empty((0, 0)),
            estimation_analysis.estimation_convergence_checker(),
            np.empty((0, 0)),
            np.empty(0),
            True,
            outliers,
            True,
        )
    else:
        estimation_input = estimation_analysis.EstimationInput(
            sample_dataset,
            outlier_rejection_settings=outliers,
            apply_apriori_parameter_deviation=True,
        )
    assert estimation_input.outlier_rejection_settings is outliers
    assert estimation_input.observation_dataset is sample_dataset


def test_collection_estimation_input_does_not_support_outlier_rejection(sample_dataset):
    """Estimation input created from an ObservationCollection keeps its original signature, without outlier rejection."""
    estimation_input = estimation_analysis.EstimationInput(
        _legacy_collection(sample_dataset),
        np.empty((0, 0)),
        estimation_analysis.estimation_convergence_checker(),
        np.empty((0, 0)),
        np.empty(0),
        True,
        True,
    )
    with pytest.raises(RuntimeError):
        estimation_input.outlier_rejection_settings = (
            estimation_analysis.simple_outlier_rejection_settings(10.0)
        )
