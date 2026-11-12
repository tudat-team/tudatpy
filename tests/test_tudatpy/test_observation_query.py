import numpy as np
import pytest

from tudatpy.estimation import observations
from tudatpy.estimation.observable_models_setup import links, model_settings


def _link_ends(receiver_body):
    """Create a minimal two-way link-end map used by the query fixtures."""
    return {
        links.transmitter: links.LinkEndId("Probe", ""),
        links.receiver: links.LinkEndId(receiver_body, ""),
    }


def _new_dataset_single_set(observable_type, receiver_body, observation_values, times):
    """Create a one-set ObservationDataset using the dataset-native API."""
    dataset = observations.ObservationDataset()
    dataset.add_observation_set(
        observable_type,
        links.LinkDefinition(_link_ends(receiver_body)),
        [np.asarray(value, dtype=float) for value in observation_values],
        times,
        links.receiver,
    )
    return dataset


@pytest.fixture
def sample_dataset():
    """Build a deterministic two-set dataset with scalar and vector observables."""
    dataset = _new_dataset_single_set(
        model_settings.one_way_range_type,
        "Earth",
        [[10.0], [20.0], [30.0]],
        [1.0, 2.0, 4.0],
    )
    dataset.set_residuals_for_set(
        0,
        [np.array([0.1]), np.array([5.0]), np.array([-0.2])],
    )

    angular_dataset = _new_dataset_single_set(
        model_settings.angular_position_type,
        "Mars",
        [[1.0, 2.0], [3.0, 4.0]],
        [3.0, 5.0],
    )

    # The second source dataset should be appended as set id 1; many row-order
    # assertions below depend on this fixed two-set layout.
    assert dataset.add_observation_set_from_dataset(angular_dataset, 0) == 1

    dataset.set_residuals_for_set(
        1,
        [np.array([0.01, 0.02]), np.array([10.0, 0.04])],
    )
    dataset.set_weight_vector_for_set(1, np.array([2.0, 3.0, 4.0, 5.0]))

    return dataset


def _to_dense_matrix(matrix):
    """Convert a scipy sparse or numpy-like matrix returned by pybind to a dense array."""
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


def test_query_conditions_select_expected_dataset_rows(sample_dataset):
    """Check primitive conditions against actual ObservationDataset rows."""
    observation_query = observations.observation_query
    dataset = sample_dataset

    # The fixture should contain one scalar range set and one vector angular set.
    assert dataset.number_of_observation_sets == 2
    # The fixture should expose five observation rows/events.
    assert dataset.number_of_observations == 5
    # The fixture should expose seven scalar components after vector expansion.
    assert dataset.total_scalar_size == 7

    # set_id == 0 should select the three Earth range rows.
    assert dataset.get_observation_ids(observation_query.set_id == 0) == [0, 1, 2]
    # set_id == 1 should select the two Mars angular-position rows.
    assert dataset.get_observation_ids(observation_query.set_id == 1) == [3, 4]
    # one_way_range should select only the scalar range observation set.
    assert dataset.get_observation_ids(
        observation_query.observable_type == model_settings.one_way_range_type
    ) == [0, 1, 2]
    # angular_position should select only the two vector-valued observations.
    assert dataset.get_observation_ids(
        observation_query.observable_type == model_settings.angular_position_type
    ) == [3, 4]
    # Receiver Earth should match the range set link definition.
    assert dataset.get_observation_ids(
        observation_query.receiver == links.LinkEndId("Earth", "")
    ) == [0, 1, 2]
    # Receiver Mars should match the angular-position set link definition.
    assert dataset.get_observation_ids(
        observation_query.receiver == links.LinkEndId("Mars", "")
    ) == [3, 4]
    # A link end type should select the rows of every set with a link end of that type.
    assert dataset.get_observation_ids(
        observations.ObservationSelectionCondition.link_end_type(links.receiver)
    ) == [0, 1, 2, 3, 4]
    assert (
        dataset.get_observation_ids(
            observations.ObservationSelectionCondition.link_end_type(links.observer)
        )
        == []
    )
    # The full link definition of the Mars set should select only its rows.
    assert dataset.get_observation_ids(
        observation_query.link_definition == links.LinkDefinition(_link_ends("Mars"))
    ) == [3, 4]
    # Inclusive time bounds [2, 4] should select rows at times 2, 3 and 4.
    assert dataset.get_observation_ids(observation_query.time.between(2.0, 4.0)) == [
        1,
        2,
        3,
    ]
    # Observation-value thresholds should select rows where any component exceeds 25.
    assert dataset.get_observation_ids(observation_query.observation.abs_greater_than(25.0)) == [2]
    # Residual thresholds should select the intentionally large residual rows.
    assert dataset.get_observation_ids(observation_query.residual.abs_greater_than(1.0)) == [1, 4]


def test_query_conditions_compose_when_selecting_dataset_rows(sample_dataset):
    """Check AND, OR and NOT composition against real dataset filtering."""
    observation_query = observations.observation_query

    # AND should intersect set-id and time conditions, leaving rows 1 and 2.
    assert sample_dataset.get_observation_ids(
        (observation_query.set_id == 0)
        & (observation_query.time >= 2.0)
        & (observation_query.time <= 4.0)
    ) == [1, 2]
    # OR should union all Earth rows with the angular row containing value 4.
    assert sample_dataset.get_observation_ids(
        (observation_query.receiver == links.LinkEndId("Earth", ""))
        | (observation_query.observation.abs_greater_than(np.array([3.5])))
    ) == [0, 1, 2, 4]
    # NOT should remove the angular row whose first residual component is large.
    assert sample_dataset.get_observation_ids(
        (observation_query.observable_type == model_settings.angular_position_type)
        & ~(observation_query.residual.abs_greater_than(np.array([1.0, 1.0])))
    ) == [3]


def test_query_conditions_drive_snapshot_getters(sample_dataset):
    """Verify query conditions select aligned fields in detached snapshots."""
    query = observations.observation_query
    selection = (query.set_id == 1) | (query.time < 2.0)
    data = sample_dataset.get_data(
        selection,
        fields=("observations", "observation_ids", "set_ids", "rows", "scalar_components"),
        ordering="estimation",
    )
    assert data["observation_ids"] == [0, 3, 4]
    assert data["set_ids"] == [0, 1, 1]
    assert data["rows"][1]["set_id"] == 1
    np.testing.assert_array_equal(data["observations"][1], [1.0, 2.0])
    np.testing.assert_array_equal(np.concatenate(data["observations"]), [10, 1, 2, 3, 4])
    assert data["scalar_components"] == [(0, 0), (3, 0), (3, 1), (4, 0), (4, 1)]


def test_python_sparse_weight_block_binding_materializes_off_diagonal_weights(sample_dataset):
    """Check Python bindings for sparse/off-diagonal dataset weight blocks."""
    angular_ids = sample_dataset.observation_ids_for_set(1)
    cross_weight_block = np.array([[0.2, 0.3], [0.4, 0.5]])

    sample_dataset.set_weight_block(
        row_observation_ids=[angular_ids[0]],
        column_observation_ids=[angular_ids[1]],
        weight_block=cross_weight_block,
    )

    # The dataset-level flag verifies that the Python set_weight_block call stored advanced blocks.
    assert sample_dataset.has_extra_weight_blocks is True
    vector_data = sample_dataset.observation_vector_data(True)
    # A cross-observation block must mark the vector data as non-diagonal.
    assert vector_data.has_off_diagonal_weights is True

    dense_weight_matrix = _to_dense_matrix(vector_data.sparse_weight_matrix)
    first_row_start = vector_data.vector_row(angular_ids[0], 0)
    second_row_start = vector_data.vector_row(angular_ids[1], 0)

    # The requested block must appear at the vector rows of the selected observations.
    np.testing.assert_allclose(
        dense_weight_matrix[
            first_row_start : first_row_start + 2, second_row_start : second_row_start + 2
        ],
        cross_weight_block,
    )
    # Symmetric insertion must materialize the transposed block in the opposite quadrant.
    np.testing.assert_allclose(
        dense_weight_matrix[
            second_row_start : second_row_start + 2, first_row_start : first_row_start + 2
        ],
        cross_weight_block.T,
    )


def test_python_snapshot_estimation_order():
    """Check snapshot estimation order against the observation vector data."""
    angular_dataset = _new_dataset_single_set(
        model_settings.angular_position_type,
        "Mars",
        [[1.0, 2.0], [3.0, 4.0]],
        [3.0, 4.0],
    )
    range_dataset = _new_dataset_single_set(
        model_settings.one_way_range_type,
        "Earth",
        [[10.0], [20.0]],
        [1.0, 2.0],
    )
    assert angular_dataset.add_observation_set_from_dataset(range_dataset, 0) == 1

    assert angular_dataset.get_observation_ids() == [0, 1, 2, 3]
    assert angular_dataset.get_observation_ids(ordering="estimation") == [2, 3, 0, 1]
    values = angular_dataset.get_observations(ordering="estimation")
    np.testing.assert_array_equal(np.concatenate(values), [10, 20, 1, 2, 3, 4])
    vector_data = angular_dataset.observation_vector_data()
    np.testing.assert_array_equal(np.concatenate(values), vector_data.observation_vector)

    # Scalar-aligned link IDs follow the same estimation order as the observation vector.
    assert vector_data.link_definition_ids == [1, 1, 0, 0, 0, 0]
    # Dataset bounds and observable-specific links remain available without a legacy collection.
    assert [bound.to_float() for bound in angular_dataset.observation_time_bounds] == [
        1.0,
        4.0,
    ]
    assert angular_dataset.link_definitions_for_observable(model_settings.one_way_range_type) == [
        range_dataset.link_definition(0)
    ]


def test_query_conditions_drive_rejection_restoration_and_filtered_datasets(
    sample_dataset,
):
    """Check active/rejected query semantics and reduced dataset creation."""
    observation_query = observations.observation_query

    sample_dataset.reject_observations(observation_query.residual.abs_greater_than(1.0), "large")
    # Rejected query should select exactly the rows with intentionally large residuals.
    assert sample_dataset.get_observation_ids(observation_query.rejected) == [1, 4]
    # Active query should select all rows not rejected by the previous operation.
    assert sample_dataset.get_observation_ids(observation_query.active) == [0, 2, 3]
    # Row metadata should expose the inactive state after rejection.
    assert sample_dataset.observation_row(1).is_active is False
    # Row metadata should preserve the rejection reason for diagnostics.
    assert sample_dataset.observation_row(1).rejection_reason == "large"

    # Estimation vector data should exclude rejected rows by default.
    assert sample_dataset.observation_vector_data().observation_ids == [
        0,
        2,
        3,
        3,
    ]
    # include_rejected=True should keep all scalar components, including rejected rows.
    assert sample_dataset.observation_vector_data(True).observation_ids == [
        0,
        1,
        2,
        3,
        3,
        4,
        4,
    ]

    kept = sample_dataset.create_new_and_keep(observation_query.active)
    # create_new_and_keep(observation_query.active) should physically contain only active rows.
    assert kept.number_of_observations == 3
    # The kept range set should retain rows 0 and 2 after dropping rejected row 1.
    np.testing.assert_allclose(kept.observation_vector_for_set(0), [10.0, 30.0])
    # The kept angular set should retain only row 3 after dropping rejected row 4.
    np.testing.assert_allclose(kept.observation_vector_for_set(1), [1.0, 2.0])

    dropped = sample_dataset.create_new_and_drop(observation_query.rejected)
    # Dropping rejected rows should produce the same active-only row count.
    assert dropped.number_of_observations == 3
    # The dropped range set should match the active-only range values.
    np.testing.assert_allclose(dropped.observation_vector_for_set(0), [10.0, 30.0])
    # The dropped angular set should match the active-only angular values.
    np.testing.assert_allclose(dropped.observation_vector_for_set(1), [1.0, 2.0])

    sample_dataset.restore_observations(observation_query.set_id == 0)
    # Restoring set 0 should leave only the rejected angular row 4.
    assert sample_dataset.get_observation_ids(observation_query.rejected) == [4]
    # Restored row 1 should become active again.
    assert sample_dataset.observation_row(1).is_active is True
    # Physical deletion preserves all surviving row identities.
    sample_dataset.delete_rejected_observations()
    # After physical removal, no rejected rows should remain selectable.
    assert sample_dataset.get_observation_ids(observation_query.rejected) == []
    # The dataset should now contain the four rows that were not rejected at removal time.
    assert sample_dataset.number_of_observations == 4
    assert sample_dataset.get_observation_ids(~observation_query.rejected) == [
        0,
        1,
        2,
        3,
    ]


def test_query_conditions_can_be_used_by_weight_api(sample_dataset):
    """Check that query-built conditions drive the canonical weight API."""
    observation_query = observations.observation_query

    sample_dataset.set_constant_single_observation_scalar_weight(
        condition=observation_query.set_id == 0,
        weight=9.0,
    )
    for observation_id in sample_dataset.observation_ids_for_set(0):
        # A set-id condition should select every row in the range set.
        np.testing.assert_allclose(sample_dataset.weight_value(observation_id), [9.0])

    sample_dataset.set_constant_single_observation_diagonal_weight(
        condition=observation_query.observable_type == model_settings.angular_position_type,
        weight=np.array([7.0, 8.0]),
    )
    for observation_id in sample_dataset.observation_ids_for_set(1):
        # An observable-type condition should select both angular rows and apply a diagonal vector.
        np.testing.assert_allclose(
            sample_dataset.weight_value(observation_id),
            [7.0, 8.0],
        )

    sample_dataset.set_constant_single_observation_scalar_weight(
        weight=6.0,
        condition=(observation_query.receiver == links.LinkEndId("Earth", ""))
        & (observation_query.time >= 2.0),
    )
    # Row 0 should keep its previous weight because it fails the time condition.
    np.testing.assert_allclose(sample_dataset.weight_value(0), [9.0])
    # Row 1 should receive the condition-selected scalar weight.
    np.testing.assert_allclose(sample_dataset.weight_value(1), [6.0])
    # Row 2 should receive the condition-selected scalar weight.
    np.testing.assert_allclose(sample_dataset.weight_value(2), [6.0])


def test_dataset_remove_observations_with_residual_condition(sample_dataset):
    """Check physical removal with the residual condition used by MPC-style examples."""
    observations_before = sample_dataset.number_of_observations
    scalar_size_before = sample_dataset.total_scalar_size

    sample_dataset.remove_observations(
        observations.ObservationSelectionCondition.residual_absolute_value_greater_than(
            np.array([1.0])
        )
    )

    assert sample_dataset.number_of_observations == observations_before - 2
    assert sample_dataset.total_scalar_size == scalar_size_before - 3
    assert (
        sample_dataset.get_observation_ids(
            observations.observation_query.residual.abs_greater_than(np.array([1.0]))
        )
        == []
    )


@pytest.mark.parametrize("epoch, expected_ids", [(1.0, [0]), (3.0, [3]), (7.0, [])])
def test_time_equality_and_inequality_are_complementary(sample_dataset, epoch, expected_ids):
    """Verify time equality and inequality queries partition dataset rows."""
    query = observations.observation_query
    equal = sample_dataset.get_observation_ids(query.time == epoch)
    unequal = sample_dataset.get_observation_ids(query.time != epoch)
    assert equal == expected_ids
    assert unequal == [row for row in range(5) if row not in expected_ids]
    assert sample_dataset.get_observation_ids(~(query.time == epoch)) == unequal


def test_python_row_and_metadata_descriptions_survive_dataset_removal(sample_dataset):
    """Verify Python row and metadata snapshots survive dataset removal."""
    dataset = sample_dataset
    row = dataset.observation_rows[0]
    metadata = dataset.observation_set_metadata[0]
    original_id = row.observation_id
    original_time = row.time
    original_type = metadata.observable_type
    dataset.remove_observations(
        observations.observation_query.active | observations.observation_query.rejected
    )
    del dataset
    assert row.observation_id == original_id
    assert row.time == original_time
    assert metadata.observable_type == original_type
