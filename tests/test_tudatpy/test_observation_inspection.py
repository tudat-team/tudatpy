"""Snapshot ownership, precise times, and shared estimator ordering at the Python boundary."""

import gc

import numpy as np
import pytest

from tudatpy.astro.time_representation import Time
from tudatpy.estimation import observations as obs
from tudatpy.estimation.observations_setup import ancillary_settings as ancillary
from tudatpy.estimation.observations_setup import observations_dependent_variables as dependent
from tudatpy.estimation.observations_setup import observations_simulation_settings as simulation

FIELDS = (
    "times",
    "observations",
    "residuals",
    "observation_ids",
    "set_ids",
    "rows",
    "dependent_variables",
    "scalar_components",
    "weight_diagonal",
    "weight_matrix",
    "metadata",
)


def link(station):
    return obs.LinkDefinition(
        {
            obs.transmitter: obs.LinkEndId("Earth", station),
            obs.receiver: obs.LinkEndId("Vehicle", ""),
        }
    )


def dataset_with_mixed_rows():
    dataset = obs.ObservationDataset()
    dataset.add_observation_set(
        obs.angular_position, link("B"), [[10, 11], [20, 21]], [4, 4], obs.receiver
    )
    dataset.add_observation_set(obs.one_way_range, link("B"), [[30], [40]], [3, 1], obs.receiver)
    dataset.add_observation_set(obs.one_way_range, link("A"), [[50]], [4], obs.receiver)
    dataset.add_observations_to_set(0, [[60, 61]], [2], sort_observations=True)
    for set_id in range(3):
        values = dataset.observations_for_set(set_id)
        dataset.set_residuals_for_set(set_id, [0.01 * value for value in values])
        dataset.set_dependent_variables_for_set(
            set_id,
            [[100 + observation_id] for observation_id in dataset.observation_ids_for_set(set_id)],
        )
    weights = np.diag(np.arange(10.0, 19.0))
    for i in range(9):
        for j in range(i + 1, 9):
            weights[i, j] = weights[j, i] = 0.001 * (10 * i + j)
    dataset.set_weight_block(list(range(6)), list(range(6)), weights)
    dataset.reject_observations(obs.observation_query.time == 4, "equal epochs")
    return dataset, weights


@pytest.mark.parametrize("ordering", ["internal", "estimation"])
@pytest.mark.parametrize("selection_kind", ["all", "active", "noncontiguous", "empty"])
def test_ordering_membership_alignment_and_weight_axes(ordering, selection_kind):
    dataset, weights = dataset_with_mixed_rows()
    query = obs.observation_query
    condition = {
        "all": obs.ObservationSelectionCondition.all(),
        "active": query.active,
        "noncontiguous": query.time != 1,
        "empty": query.time < 0,
    }[selection_kind]
    # A noncontiguous event selection independent of output ordering.
    if selection_kind == "noncontiguous":
        condition = (query.time != 1) & ~(
            query.observation.abs_greater_than(19) & ~query.observation.abs_greater_than(22)
        )
    original = dataset.get_data(fields=FIELDS)
    data = dataset.get_data(condition, fields=FIELDS, ordering=ordering)
    ids = data["observation_ids"]
    expected_sequence = [0, 1, 2, 3, 4, 5] if ordering == "internal" else [4, 2, 3, 5, 0, 1]
    selected_ids = dataset.observation_ids_matching_condition(condition)
    assert ids == [event for event in expected_sequence if event in selected_ids]
    assert sorted(ids) == sorted(dataset.get_observation_ids(condition, ordering="internal"))
    assert set(data["metadata"]) == set(data["set_ids"])
    scalar_indices = []
    expected_components = []
    for index, event in enumerate(ids):
        row = original["rows"][event]
        assert data["times"][index] == original["times"][event]
        assert data["rows"][index]["is_active"] == row["is_active"]
        assert data["set_ids"][index] == row["set_id"]
        for field in ("observations", "residuals", "dependent_variables"):
            np.testing.assert_array_equal(data[field][index], original[field][event])
        scalar_indices.extend(
            range(row["first_scalar_component"], row["first_scalar_component"] + row["scalar_size"])
        )
        expected_components.extend((event, component) for component in range(row["scalar_size"]))
    assert data["scalar_components"] == expected_components
    np.testing.assert_array_equal(data["weight_diagonal"], np.diag(weights)[scalar_indices])
    np.testing.assert_array_equal(
        data["weight_matrix"].toarray(), weights[np.ix_(scalar_indices, scalar_indices)]
    )
    for field in FIELDS:
        single = getattr(dataset, "get_" + field)(condition, ordering=ordering)
        if field == "weight_matrix":
            np.testing.assert_array_equal(single.toarray(), data[field].toarray())
        elif field == "weight_diagonal":
            np.testing.assert_array_equal(single, data[field])
        elif field in ("observations", "residuals", "dependent_variables"):
            for a, b in zip(single, data[field]):
                np.testing.assert_array_equal(a, b)
        elif field not in ("rows", "metadata"):
            assert single == data[field]
    if ordering == "estimation":
        projection = dataset.create_new_and_keep(condition).create_estimation_projection(
            include_rejected=True
        )
        for field, projection_field in (
            ("observations", "observation_vector"),
            ("residuals", "residual_vector"),
        ):
            actual = np.concatenate(data[field]) if ids else np.empty(0)
            np.testing.assert_array_equal(actual, getattr(projection, projection_field))
        assert [event for event, _ in expected_components] == projection.observation_ids
        np.testing.assert_array_equal(
            data["weight_matrix"].toarray(), projection.sparse_weight_matrix.toarray()
        )
    assert dataset.get_observation_ids() == original["observation_ids"]
    for a, b in zip(dataset.get_observations(), original["observations"]):
        np.testing.assert_array_equal(a, b)


@pytest.mark.parametrize("ordering", ["internal", "estimation"])
@pytest.mark.parametrize("multi", [False, True])
def test_snapshots_survive_all_mutations_and_dataset_destruction(ordering, multi):
    dataset, _ = dataset_with_mixed_rows()
    if multi:
        saved = dataset.get_data(fields=FIELDS, ordering=ordering)
    else:
        saved = {field: getattr(dataset, "get_" + field)(ordering=ordering) for field in FIELDS}
    ids = saved["observation_ids"].copy()
    first = ids.index(0)
    saved_matrix = saved["weight_matrix"].toarray().copy()
    saved_diagonal = saved["weight_diagonal"].copy()
    dataset.set_observations_for_set(0, [[160, 161], [110, 111], [120, 121]])
    dataset.set_residuals_for_set(0, [[6, 6.1], [1, 1.1], [2, 2.1]])
    dataset.set_dependent_variables_for_set(0, [[160], [110], [120]])
    dataset.set_constant_single_observation_diagonal_weight_for_set(0, [40, 50])
    dataset.set_weight_block([0], [2], [[0.7], [0.8]])
    dataset.set_link_end_reference_point("Earth", "Changed", obs.transmitter)
    dataset.restore_observations(obs.observation_query.rejected)
    dataset.reject_observations(obs.observation_query.active, "new rejection")
    dataset.restore_observations(obs.observation_query.rejected)
    dataset.add_observations_to_set(0, [[70, 71]], [0], sort_observations=True)
    dataset.remove_observations(obs.observation_query.time == 1)
    assert 3 not in dataset.get_observation_ids()
    assert 6 in dataset.get_observation_ids()
    np.testing.assert_array_equal(dataset.observation_value(0), [110, 111])
    assert saved["observation_ids"] == ids
    np.testing.assert_array_equal(saved["observations"][first], [10, 11])
    np.testing.assert_array_equal(saved["residuals"][first], [0.1, 0.11])
    np.testing.assert_array_equal(saved["dependent_variables"][first], [100])
    assert saved["rows"][first]["rejection_reason"] == "equal epochs"
    assert not saved["rows"][first]["is_active"]
    assert saved["times"][first] == Time(4)
    np.testing.assert_array_equal(saved["weight_matrix"].toarray(), saved_matrix)
    np.testing.assert_array_equal(saved["weight_diagonal"], saved_diagonal)
    # Editing nested snapshots must also leave the dataset untouched.
    saved["observations"][first][0] = -10
    saved["residuals"][first][0] = -20
    saved["dependent_variables"][first][0] = -30
    saved["rows"][first]["dependent_variables"][0] = -40
    saved["rows"][first]["rejection_reason"] = "edited"
    saved["weight_diagonal"][0] = -50
    saved["weight_matrix"].data[:] = -60
    saved["metadata"][0]["observable_size"] = 99
    np.testing.assert_array_equal(dataset.observation_value(0), [110, 111])
    np.testing.assert_array_equal(dataset.residual_value(0), [1, 1.1])
    np.testing.assert_array_equal(dataset.dependent_variables(0), [110])
    np.testing.assert_array_equal(dataset.weight_value(0), [40, 50])
    assert dataset.get_metadata()[0]["observable_size"] == 2
    del dataset
    gc.collect()
    assert saved["observations"][first][0] == -10
    assert saved["rows"][first]["dependent_variables"][0] == -40
    assert saved["times"][first] == Time(4)
    assert (
        saved["metadata"][0]["link_definition"].link_end_id(obs.transmitter).reference_point == "B"
    )
    assert saved["weight_matrix"].toarray().shape == saved_matrix.shape


def test_nested_metadata_is_detached():
    definition = link("Station")
    settings = simulation.tabulated_simulation_settings(obs.one_way_range, definition, [1.0])
    setting = dependent.elevation_angle_dependent_variable(
        obs.transmitter, obs.LinkEndId("Earth", "Station")
    )
    settings.add_dependent_variables([setting])
    ancillary_value = ancillary.ObservationAncillarySimulationSettings()
    ancillary_value.set_float_settings(ancillary.doppler_integration_time, 10)
    dataset = obs.ObservationDataset()
    dataset.add_observation_set(
        obs.one_way_range,
        definition,
        [[1]],
        [1],
        obs.receiver,
        dependent_variables=[[0.5]],
        dependent_variable_bookkeeping=settings.dependent_variable_bookkeeping,
        ancillary_settings=ancillary_value,
    )
    metadata = dataset.get_metadata()
    captured = metadata[0]["ancillary_settings"]
    source = dataset.ancillary_settings_for_set(0)
    source.set_float_settings(ancillary.doppler_integration_time, 20)
    assert captured.get_float_settings(ancillary.doppler_integration_time) == 10
    captured.set_float_settings(ancillary.doppler_integration_time, 30)
    assert source.get_float_settings(ancillary.doppler_integration_time) == 20
    layout = metadata[0]["dependent_variable_layout"]
    assert list(layout) == [(0, 1)]
    assert layout[(0, 1)] == setting
    assert layout[(0, 1)] is not setting
    layout.clear()
    assert list(dataset.get_metadata()[0]["dependent_variable_layout"]) == [(0, 1)]
    del dataset, source, settings, setting
    gc.collect()
    assert captured.get_float_settings(ancillary.doppler_integration_time) == 30


@pytest.mark.parametrize("ordering", ["internal", "estimation"])
def test_precise_times_and_empty_missing_fields(ordering):
    dataset = obs.ObservationDataset()
    data = dataset.get_data(fields=FIELDS, ordering=ordering)
    assert all(
        data[field] == []
        for field in FIELDS
        if field not in ("metadata", "weight_diagonal", "weight_matrix")
    )
    assert data["metadata"] == {}
    assert data["weight_matrix"].shape == (0, 0)
    assert data["weight_diagonal"].size == 0
    t0, t1 = Time(1000000, 1e-9), Time(1000000, 2e-9)
    assert float(t0) == float(t1)
    assert t0 != t1
    dataset.add_observation_set(obs.one_way_range, link("A"), [[1], [2]], [t1, t0], obs.receiver)
    saved = dataset.get_data(fields=FIELDS, ordering=ordering)
    assert saved["times"] == [t1, t0]
    assert all(isinstance(t, Time) for t in saved["times"])
    assert saved["rows"][0]["time"] == t1
    assert saved["rows"][1]["time"] == t0
    np.testing.assert_array_equal(np.concatenate(saved["residuals"]), [0, 0])
    assert all(value.size == 0 for value in saved["dependent_variables"])
    assert saved["metadata"][0]["ancillary_settings"] is None
    assert saved["metadata"][0]["dependent_variable_layout"] == {}
    assert dataset.get_times(obs.observation_query.time < 0, ordering=ordering) == []
    dataset.remove_observations(obs.observation_query.time == t0)
    assert dataset.get_times(ordering=ordering) == [t1]
    del dataset
    gc.collect()
    assert saved["times"] == [t1, t0]


def test_invalid_ordering_fields_and_removed_api():
    dataset = obs.ObservationDataset()
    for field in FIELDS:
        with pytest.raises(ValueError, match="ordering"):
            getattr(dataset, "get_" + field)(ordering="time")
    with pytest.raises(ValueError, match="ordering"):
        dataset.get_data(ordering="time")
    with pytest.raises(ValueError, match="Unknown"):
        dataset.get_data(fields=("unknown",))
    with pytest.raises(ValueError, match="Duplicate"):
        dataset.get_data(fields=("times", "times"))
    assert dataset.get_data(fields=()) == {}
    assert not hasattr(dataset, "create_viewer")
    assert not hasattr(obs, "ObservationDatasetViewer")
