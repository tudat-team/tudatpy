import warnings

import numpy as np
import pytest

from tudatpy.estimation import observations


def _link_ends(receiver_body):
    """Create a minimal two-way link-end map used by the query fixtures."""
    return {
        observations.transmitter: observations.LinkEndId("Probe", ""),
        observations.receiver: observations.LinkEndId(receiver_body, ""),
    }


def _converted_single_set(observable_type, receiver_body, observation_values, times):
    """Create a dataset through the legacy facade to exercise compatibility paths."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        observation_set = observations.create_single_observation_set(
            observable_type,
            _link_ends(receiver_body),
            [np.asarray(value, dtype=float) for value in observation_values],
            times,
            observations.receiver,
        )
        return observations.create_observation_dataset_from_single_observation_set(observation_set)


@pytest.fixture
def sample_dataset():
    """Build a deterministic two-set dataset with scalar and vector observables."""
    dataset = _converted_single_set(
        observations.one_way_range,
        "Earth",
        [[10.0], [20.0], [30.0]],
        [1.0, 2.0, 4.0],
    )
    dataset.set_residuals_for_set(
        0,
        [np.array([0.1]), np.array([5.0]), np.array([-0.2])],
    )

    angular_dataset = _converted_single_set(
        observations.angular_position,
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


def _assert_link_end_id(link_end_id, body_name, reference_point=""):
    """Assert both parts of a link-end id so selector tests do not compare identity."""
    # The body name verifies that the query stored the intended physical body.
    assert link_end_id.body_name == body_name
    # The reference point verifies that station-level ids survive round-tripping.
    assert link_end_id.reference_point == reference_point


def _to_dense_matrix(matrix):
    """Convert a scipy sparse or numpy-like matrix returned by pybind to a dense array."""
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


def test_observation_query_builds_basic_backend_conditions():
    """Check that simple observation_query selectors create inspectable backend conditions."""
    observation_query = observations.observation_query

    condition = observation_query.set_id == 4
    # The query helper must produce the C++ condition type consumed by datasets.
    assert isinstance(condition, observations.ObservationSelectionCondition)
    # Equality on observation_query.set_id should map to the dedicated set-id condition node.
    assert condition.condition_type == observations.ObservationSelectionConditionType.set_id
    # The set-id value must be preserved so downstream filtering selects set 4.
    assert condition.set_id_value == 4

    negated = observation_query.set_id != 4
    # Inequality should be represented as logical NOT, not as a separate ad hoc type.
    assert negated.condition_type == observations.ObservationSelectionConditionType.not_condition
    # The NOT child should be the same primitive set-id selector used by equality.
    assert negated.child_conditions[0].condition_type == (
        observations.ObservationSelectionConditionType.set_id
    )
    # The child condition must retain the rejected set id after negation.
    assert negated.child_conditions[0].set_id_value == 4

    observable_condition = observation_query.observable_type == observations.one_way_range
    # Observable-type equality should build the canonical observable-type node.
    assert observable_condition.condition_type == (
        observations.ObservationSelectionConditionType.observable_type
    )
    # The enum value must round-trip because observable type drives set selection.
    assert observable_condition.observable_type_value == observations.one_way_range

    link_definition = observations.LinkDefinition(_link_ends("Earth"))
    link_definition_condition = observation_query.link_definition == link_definition
    # Link-definition equality should build the canonical link-definition node.
    assert link_definition_condition.condition_type == (
        observations.ObservationSelectionConditionType.link_definition
    )
    # The receiver id inside the stored link definition verifies deep value storage.
    _assert_link_end_id(
        link_definition_condition.link_definition_value.link_end_id(observations.receiver),
        "Earth",
    )


@pytest.mark.parametrize(
    ("condition", "condition_type", "time_value"),
    [
        (
            observations.observation_query.time >= 1.0,
            observations.ObservationSelectionConditionType.time_greater_equal,
            1.0,
        ),
        (
            observations.observation_query.time > 2.0,
            observations.ObservationSelectionConditionType.time_greater_than,
            2.0,
        ),
        (
            observations.observation_query.time <= 3.0,
            observations.ObservationSelectionConditionType.time_less_equal,
            3.0,
        ),
        (
            observations.observation_query.time < 4.0,
            observations.ObservationSelectionConditionType.time_less_than,
            4.0,
        ),
    ],
)
def test_observation_query_builds_one_sided_time_conditions(condition, condition_type, time_value):
    """Check that each Python comparison operator maps to the right time node."""
    # The condition type distinguishes inclusive and exclusive time comparisons.
    assert condition.condition_type == condition_type
    # The scalar time value must be stored unchanged for later row filtering.
    assert condition.time_value == time_value


def test_observation_query_builds_time_bounds_and_logical_trees():
    """Check bounded-time queries and logical composition tree structure."""
    observation_query = observations.observation_query

    bounded = observation_query.time.between(1.0, 2.0)
    # between() should use the explicit time-bounds node, not two comparisons.
    assert bounded.condition_type == observations.ObservationSelectionConditionType.time_bounds
    # Both bounds must be stored in order because the C++ side applies them directly.
    assert bounded.time_bounds_value == (1.0, 2.0)

    composed = (
        (observation_query.set_id == 4) & (observation_query.time < 10.0)
    ) | ~observation_query.rejected
    # Top-level OR verifies that Python | creates a composable backend tree.
    assert composed.condition_type == observations.ObservationSelectionConditionType.or_condition
    # The OR node should preserve exactly the two operands supplied by the user.
    assert len(composed.child_conditions) == 2
    # The first child should be the parenthesized AND expression.
    assert composed.child_conditions[0].condition_type == (
        observations.ObservationSelectionConditionType.and_condition
    )
    # The second child should be the explicit negation of observation_query.rejected.
    assert composed.child_conditions[1].condition_type == (
        observations.ObservationSelectionConditionType.not_condition
    )


@pytest.mark.parametrize(
    ("selector_name", "link_end_type"),
    [
        ("receiver", observations.receiver),
        ("transmitter", observations.transmitter),
        ("observer", observations.observer),
        ("observed_body", observations.observed_body),
    ],
)
def test_observation_query_link_end_aliases(selector_name, link_end_type):
    """Check that named link-end aliases are equivalent to explicit link_end()."""
    station = observations.LinkEndId("Earth", "DSS-63")

    alias_condition = getattr(observations.observation_query, selector_name) == station
    explicit_condition = observations.ObservationSelectionCondition.link_end(link_end_type, station)

    # Alias selectors should create the same condition node type as link_end().
    assert alias_condition.condition_type == explicit_condition.condition_type
    # The alias must map to the intended LinkEndType enum.
    assert alias_condition.link_end_type_value == explicit_condition.link_end_type_value
    # The station id should be stored by value, including reference point.
    _assert_link_end_id(alias_condition.link_end_id_value, "Earth", "DSS-63")

    negated_alias = getattr(observations.observation_query, selector_name) != station
    # Inequality on aliases should use normal NOT composition.
    assert negated_alias.condition_type == (
        observations.ObservationSelectionConditionType.not_condition
    )
    # The negated child must still use the correct link-end role.
    assert negated_alias.child_conditions[0].link_end_type_value == link_end_type


def test_observation_query_value_threshold_conditions_store_scalar_and_vector_limits():
    """Check threshold selectors for scalar and vector-valued limits."""
    observation_query = observations.observation_query

    residual_condition = observation_query.residual.abs_greater_than(3.0)
    # Residual thresholds should build the residual absolute-value condition.
    assert residual_condition.condition_type == (
        observations.ObservationSelectionConditionType.residual_absolute_value_greater_than
    )
    # Scalar thresholds are stored as a length-one vector for the C++ evaluator.
    np.testing.assert_allclose(residual_condition.vector_limit, [3.0])

    observation_condition = observation_query.observation.abs_greater_than(np.array([1.0, 2.0]))
    # Observation thresholds should build the observed-value absolute-value condition.
    assert observation_condition.condition_type == (
        observations.ObservationSelectionConditionType.observation_absolute_value_greater_than
    )
    # Vector thresholds must preserve all components for vector observables.
    np.testing.assert_allclose(observation_condition.vector_limit, [1.0, 2.0])


@pytest.mark.parametrize(
    "selector_expression",
    [
        lambda observation_query: observation_query,
        lambda observation_query: observation_query.observable_type,
        lambda observation_query: observation_query.set_id,
        lambda observation_query: observation_query.link_definition,
        lambda observation_query: observation_query.link_end(observations.receiver),
        lambda observation_query: observation_query.receiver,
        lambda observation_query: observation_query.time,
        lambda observation_query: observation_query.residual,
        lambda observation_query: observation_query.observation,
    ],
)
def test_observation_query_rejects_incomplete_selector_bool_contexts(
    selector_expression,
):
    """Check that unfinished selector objects cannot be used as booleans."""
    # Incomplete selectors must fail loudly so users write comparisons/calls first.
    with pytest.raises(TypeError, match="Incomplete observation query selector"):
        bool(selector_expression(observations.observation_query))


def test_observation_query_rejects_condition_bool_contexts():
    """Check that complete conditions reject Python and/or/not semantics."""
    observation_query = observations.observation_query

    # Direct bool conversion must fail because truthiness would be ambiguous.
    with pytest.raises(
        TypeError, match="ObservationSelectionCondition objects cannot be converted to bool"
    ):
        bool(observation_query.set_id == 1)

    # Python and would call bool(); this guards the documented use of & instead.
    with pytest.raises(
        TypeError, match="ObservationSelectionCondition objects cannot be converted to bool"
    ):
        (observation_query.set_id == 1) and (observation_query.time < 2.0)


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
    assert dataset.observation_ids_matching_condition(observation_query.set_id == 0) == [0, 1, 2]
    # set_id == 1 should select the two Mars angular-position rows.
    assert dataset.observation_ids_matching_condition(observation_query.set_id == 1) == [3, 4]
    # one_way_range should select only the scalar range observation set.
    assert dataset.observation_ids_matching_condition(
        observation_query.observable_type == observations.one_way_range
    ) == [0, 1, 2]
    # angular_position should select only the two vector-valued observations.
    assert dataset.observation_ids_matching_condition(
        observation_query.observable_type == observations.angular_position
    ) == [3, 4]
    # Receiver Earth should match the range set link definition.
    assert dataset.observation_ids_matching_condition(
        observation_query.receiver == observations.LinkEndId("Earth", "")
    ) == [0, 1, 2]
    # Receiver Mars should match the angular-position set link definition.
    assert dataset.observation_ids_matching_condition(
        observation_query.receiver == observations.LinkEndId("Mars", "")
    ) == [3, 4]
    # Inclusive time bounds [2, 4] should select rows at times 2, 3 and 4.
    assert dataset.observation_ids_matching_condition(observation_query.time.between(2.0, 4.0)) == [
        1,
        2,
        3,
    ]
    # Observation-value thresholds should select rows where any component exceeds 25.
    assert dataset.observation_ids_matching_condition(
        observation_query.observation.abs_greater_than(25.0)
    ) == [2]
    # Residual thresholds should select the intentionally large residual rows.
    assert dataset.observation_ids_matching_condition(
        observation_query.residual.abs_greater_than(1.0)
    ) == [1, 4]


def test_query_conditions_compose_when_selecting_dataset_rows(sample_dataset):
    """Check AND, OR and NOT composition against real dataset filtering."""
    observation_query = observations.observation_query

    # AND should intersect set-id and time conditions, leaving rows 1 and 2.
    assert sample_dataset.observation_ids_matching_condition(
        (observation_query.set_id == 0)
        & (observation_query.time >= 2.0)
        & (observation_query.time <= 4.0)
    ) == [1, 2]
    # OR should union all Earth rows with the angular row containing value 4.
    assert sample_dataset.observation_ids_matching_condition(
        (observation_query.receiver == observations.LinkEndId("Earth", ""))
        | (observation_query.observation.abs_greater_than(np.array([3.5])))
    ) == [0, 1, 2, 4]
    # NOT should remove the angular row whose first residual component is large.
    assert sample_dataset.observation_ids_matching_condition(
        (observation_query.observable_type == observations.angular_position)
        & ~(observation_query.residual.abs_greater_than(np.array([1.0, 1.0])))
    ) == [3]


def test_query_conditions_drive_viewers_and_flattened_data(sample_dataset):
    """Check that query results drive viewers and flattened provenance correctly."""
    observation_query = observations.observation_query

    viewer = sample_dataset.create_viewer(
        (observation_query.set_id == 1) | (observation_query.time < 2.0)
    )
    # The viewer should contain row 0 plus both angular-position rows.
    assert viewer.number_of_observations == 3
    # The stored row ids verify that viewer ordering follows dataset row order.
    assert viewer.observation_ids == [0, 3, 4]

    # Viewer index 1 should refer to dataset row 3 from angular set id 1.
    assert viewer.observation_row(1).set_id == 1
    # The same viewer index should return the first angular-position value.
    np.testing.assert_allclose(viewer.observation_value(1), [1.0, 2.0])

    flattened = viewer.estimation_flattened_observation_data()
    # Flattened values should expand vector observations into scalar components.
    np.testing.assert_allclose(
        flattened.observation_vector,
        [10.0, 1.0, 2.0, 3.0, 4.0],
    )
    # Each flattened scalar should point back to the source observation row.
    assert flattened.observation_ids == [0, 3, 3, 4, 4]
    # Each flattened scalar should point back to the source observation set.
    assert flattened.set_ids == [0, 1, 1, 1, 1]
    # The per-set ordering should record first appearance in flattened data.
    assert flattened.set_ids_in_row_order == [0, 1]
    # Unique row ids for set 1 should collapse vector components back to rows.
    assert flattened.unique_observation_ids_for_set(1) == [3, 4]
    # Row lookup should map observation row 3 component 1 to flattened row 2.
    assert flattened.flattened_row(3, 1) == 2


def test_python_sparse_weight_block_binding_materializes_off_diagonal_weights(sample_dataset):
    """Check Python bindings for sparse/off-diagonal dataset weight blocks."""
    angular_ids = sample_dataset.observation_ids_for_set(1)
    cross_weight_block = np.array([[0.2, 0.3], [0.4, 0.5]])

    sample_dataset.set_weight_block(
        row_observation_ids=[angular_ids[0]],
        column_observation_ids=[angular_ids[1]],
        weight_block=cross_weight_block,
        make_symmetric=True,
    )

    # The dataset-level flag verifies that the Python set_weight_block call stored advanced blocks.
    assert sample_dataset.has_extra_weight_blocks is True
    # Symmetric insertion should store the requested block and its transpose.
    assert len(sample_dataset.extra_weight_blocks) == 2
    # The first exposed block should preserve the row scalar-component ids for the first angular row.
    assert sample_dataset.extra_weight_blocks[0].row_scalar_component_ids == [3, 4]
    # The first exposed block should preserve the column scalar-component ids for the second angular row.
    assert sample_dataset.extra_weight_blocks[0].column_scalar_component_ids == [5, 6]
    # The dense block value must round-trip through the Python ObservationWeightBlock binding.
    np.testing.assert_allclose(
        sample_dataset.extra_weight_blocks[0].weight_block, cross_weight_block
    )

    flattened = sample_dataset.estimation_flattened_observation_data(True)
    # A cross-observation block must mark the flattened data as non-diagonal.
    assert flattened.has_off_diagonal_weights is True
    # The inverse flag should also reflect the presence of off-diagonal weights.
    assert flattened.is_diagonal_weight_only is False

    dense_weight_matrix = _to_dense_matrix(flattened.sparse_weight_matrix)
    first_row_start = flattened.flattened_row(angular_ids[0], 0)
    second_row_start = flattened.flattened_row(angular_ids[1], 0)

    # The requested block must materialize at the flattened rows of the selected observations.
    np.testing.assert_allclose(
        dense_weight_matrix[
            first_row_start : first_row_start + 2, second_row_start : second_row_start + 2
        ],
        cross_weight_block,
    )
    # make_symmetric=True must materialize the transposed block in the opposite quadrant.
    np.testing.assert_allclose(
        dense_weight_matrix[
            second_row_start : second_row_start + 2, first_row_start : first_row_start + 2
        ],
        cross_weight_block.T,
    )


def test_python_viewer_ordered_flattening_reorders_selected_rows():
    """Check Python viewer ordered flattening against dataset-row flattening."""
    angular_dataset = _converted_single_set(
        observations.angular_position,
        "Mars",
        [[1.0, 2.0], [3.0, 4.0]],
        [3.0, 4.0],
    )
    range_dataset = _converted_single_set(
        observations.one_way_range,
        "Earth",
        [[10.0], [20.0]],
        [1.0, 2.0],
    )
    assert angular_dataset.add_observation_set_from_dataset(range_dataset, 0) == 1

    viewer = angular_dataset.create_viewer(observations.ObservationSelectionCondition.all())
    # The viewer itself follows dataset row insertion order: angular rows before range rows.
    assert viewer.observation_ids == [0, 1, 2, 3]

    estimation_flattened = viewer.estimation_flattened_observation_data()
    # Estimation flattening expands vector rows but keeps the selected dataset row order.
    assert estimation_flattened.observation_ids == [0, 0, 1, 1, 2, 3]

    ordered_flattened = viewer.ordered_flattened_observation_data()
    # Ordered flattening should reorder selected rows into Tudat's ordered-output convention.
    assert ordered_flattened.observation_ids == [2, 3, 0, 0, 1, 1]


def test_query_conditions_drive_rejection_restoration_and_filtered_datasets(
    sample_dataset,
):
    """Check active/rejected query semantics and reduced dataset creation."""
    observation_query = observations.observation_query

    sample_dataset.reject_observations(observation_query.residual.abs_greater_than(1.0), "large")
    # Rejected query should select exactly the rows with intentionally large residuals.
    assert sample_dataset.observation_ids_matching_condition(observation_query.rejected) == [1, 4]
    # Active query should select all rows not rejected by the previous operation.
    assert sample_dataset.observation_ids_matching_condition(observation_query.active) == [0, 2, 3]
    # Row metadata should expose the inactive state after rejection.
    assert sample_dataset.observation_row(1).is_active is False
    # Row metadata should preserve the rejection reason for diagnostics.
    assert sample_dataset.observation_row(1).rejection_reason == "large"

    # Estimation flattening should exclude rejected rows by default.
    assert sample_dataset.estimation_flattened_observation_data().observation_ids == [
        0,
        2,
        3,
        3,
    ]
    # include_rejected=True should keep all scalar components, including rejected rows.
    assert sample_dataset.estimation_flattened_observation_data(True).observation_ids == [
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
    assert sample_dataset.observation_ids_matching_condition(observation_query.rejected) == [4]
    # Restored row 1 should become active again.
    assert sample_dataset.observation_row(1).is_active is True
    # The renamed removal API should physically delete all rows that are still rejected.
    sample_dataset.remove_rejected_observations()
    # After physical removal, no rejected rows should remain selectable.
    assert sample_dataset.observation_ids_matching_condition(observation_query.rejected) == []
    # The dataset should now contain the four rows that were not rejected at removal time.
    assert sample_dataset.number_of_observations == 4
    # The old delete_* spelling should not remain on the primary ObservationDataset API.
    assert not hasattr(sample_dataset, "delete_rejected_observations")


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
        condition=observation_query.observable_type == observations.angular_position,
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
        condition=(observation_query.receiver == observations.LinkEndId("Earth", ""))
        & (observation_query.time >= 2.0),
    )
    # Row 0 should keep its previous weight because it fails the time condition.
    np.testing.assert_allclose(sample_dataset.weight_value(0), [9.0])
    # Row 1 should receive the condition-selected scalar weight.
    np.testing.assert_allclose(sample_dataset.weight_value(1), [6.0])
    # Row 2 should receive the condition-selected scalar weight.
    np.testing.assert_allclose(sample_dataset.weight_value(2), [6.0])


def test_redundant_weight_aliases_and_legacy_dataset_methods_are_not_public(sample_dataset):
    """Check that the primary ObservationDataset Python surface contains no removed aliases."""

    # Short ambiguous aliases should not be available beside the explicit scalar/diagonal/matrix names.
    assert not hasattr(sample_dataset, "set_constant_weight")
    # The diagonal alias should not be available beside set_constant_single_observation_diagonal_weight.
    assert not hasattr(sample_dataset, "set_constant_diagonal_weight")
    # The matrix alias should not be available beside set_constant_single_observation_matrix_weight.
    assert not hasattr(sample_dataset, "set_constant_matrix_weight")

    # Keyword selector overloads were removed; selection should go through ObservationSelectionCondition objects.
    with pytest.raises(TypeError):
        sample_dataset.set_constant_single_observation_scalar_weight(
            weight=1.0,
            set_id=0,
        )

    # Legacy flat mutation helpers should not be public on the primary dataset API.
    assert not hasattr(sample_dataset, "set_observation_vector_for_set")
    # Legacy flat residual mutation should likewise remain out of the primary Python surface.
    assert not hasattr(sample_dataset, "set_residual_vector_for_set")
    # Parser/filter compatibility helpers should not be exposed directly on ObservationDataset.
    assert not hasattr(sample_dataset, "filtered_observation_indices")
    # Duplicate-erasure is a legacy set-processing helper and should not be public on ObservationDataset.
    assert not hasattr(sample_dataset, "erase_duplicate_observations_from_set")
