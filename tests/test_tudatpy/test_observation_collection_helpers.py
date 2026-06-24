from types import SimpleNamespace

import numpy as np
import pytest

from tudatpy.estimation import observations
from tudatpy.estimation.observations import _observation_collection_helpers as helpers

DOPPLER = "doppler"
PIXEL = "pixel"


class FakeParser:
    def __init__(self, observable_type, use_opposite_condition=False):
        self.observable_type = observable_type
        self.use_opposite_condition = use_opposite_condition


class FakeFilter:
    def __init__(self, filter_type, filter_value, filter_out=True, use_opposite_condition=False):
        self.filter_type = filter_type
        self.filter_value = filter_value
        self.filter_out = filter_out
        self.use_opposite_condition = use_opposite_condition


class FakeObservationCollection:
    def __init__(self):
        self.observable_type_start_index_and_size = {
            DOPPLER: (0, 3),
            PIXEL: (3, 4),
        }
        self.observations = np.array([100.0, 101.0, 102.0, 1.0, 2.0, 3.0, 4.0])
        self.residuals = np.array([0.1, 0.2, 0.3, 3.0, 4.0, 5.0, 6.0])
        self.weights = np.ones(7)
        self.concatenated_times = np.array([10.0, 20.0, 30.0, 40.0, 40.0, 50.0, 50.0])
        self.observation_times = {
            DOPPLER: np.array([10.0, 20.0, 30.0]),
            PIXEL: np.array([40.0, 50.0]),
        }
        self.filter_calls = []
        self.weight_calls = []

    def _slice_for_parser(self, parser):
        start, size = self.observable_type_start_index_and_size[parser.observable_type]
        return slice(start, start + size)

    def get_concatenated_observations(self, parser):
        return self.observations[self._slice_for_parser(parser)]

    def get_concatenated_residuals(self, parser):
        return self.residuals[self._slice_for_parser(parser)]

    def get_concatenated_weights(self, parser):
        return self.weights[self._slice_for_parser(parser)]

    def get_concatenated_observation_times(self, parser):
        return self.observation_times[parser.observable_type]

    def set_constant_weight(self, weight, parser):
        self.weight_calls.append((weight, parser))

    def filter_observations(self, filter_settings, parser, save_filtered_observations=True):
        self.filter_calls.append((filter_settings, parser, save_filtered_observations))


@pytest.fixture(autouse=True)
def patch_tudat_helpers(monkeypatch):
    monkeypatch.setattr(
        helpers,
        "observation_parser",
        lambda observable_type, use_opposite_condition=False: FakeParser(
            observable_type, use_opposite_condition
        ),
    )
    monkeypatch.setattr(helpers, "observation_filter", FakeFilter)
    monkeypatch.setattr(
        helpers,
        "ObservationFilterType",
        SimpleNamespace(residual_filtering="residual_filtering"),
    )
    monkeypatch.setattr(
        helpers,
        "get_observable_size",
        lambda observable_type: {DOPPLER: 1, PIXEL: 2}[observable_type],
    )


def test_helpers_are_reexported_from_observations_package():
    assert observations.filter_residuals_by_type is helpers.filter_residuals_by_type
    assert (
        observations.get_observable_type_residual_norms
        is helpers.get_observable_type_residual_norms
    )


def test_filter_residuals_by_type_routes_to_only_that_observable_type():
    collection = FakeObservationCollection()

    helpers.filter_residuals_by_type(
        collection,
        DOPPLER,
        10.0,
        save_filtered_observations=False,
    )

    assert len(collection.filter_calls) == 1
    filter_settings, parser, save_filtered = collection.filter_calls[0]
    assert parser.observable_type == DOPPLER
    assert filter_settings.filter_type == "residual_filtering"
    assert filter_settings.filter_value == 10.0
    assert save_filtered is False


def test_set_constant_weight_from_sigma_by_type_accepts_component_sigmas():
    collection = FakeObservationCollection()

    helpers.set_constant_weight_from_sigma_by_type(collection, PIXEL, np.array([2.0, 4.0]))

    assert len(collection.weight_calls) == 1
    weight, parser = collection.weight_calls[0]
    assert parser.observable_type == PIXEL
    np.testing.assert_allclose(weight, np.array([0.25, 0.0625]))


def test_observable_type_vector_extraction_uses_collection_start_and_size():
    collection = FakeObservationCollection()
    full_vector = np.arange(7.0)

    np.testing.assert_array_equal(
        helpers.get_observable_type_residuals(collection, PIXEL, full_vector),
        np.array([3.0, 4.0, 5.0, 6.0]),
    )

    residual_history = np.column_stack((full_vector, full_vector + 10.0))
    np.testing.assert_array_equal(
        helpers.get_observable_type_residual_history(collection, PIXEL, residual_history),
        np.array([[3.0, 13.0], [4.0, 14.0], [5.0, 15.0], [6.0, 16.0]]),
    )


def test_vector_observable_residuals_can_be_reported_as_rows_norms_and_times():
    collection = FakeObservationCollection()

    np.testing.assert_array_equal(
        helpers.get_observable_type_residual_matrix(collection, PIXEL),
        np.array([[3.0, 4.0], [5.0, 6.0]]),
    )
    np.testing.assert_allclose(
        helpers.get_observable_type_residual_norms(collection, PIXEL),
        np.array([5.0, np.sqrt(61.0)]),
    )
    np.testing.assert_array_equal(
        helpers.get_observable_type_times(collection, PIXEL),
        np.array([40.0, 40.0, 50.0, 50.0]),
    )
    np.testing.assert_array_equal(
        helpers.get_observable_type_times_per_observation(collection, PIXEL),
        np.array([40.0, 50.0]),
    )
    np.testing.assert_array_equal(
        helpers.get_observable_type_times_per_observation(
            collection, PIXEL, collection.concatenated_times
        ),
        np.array([40.0, 50.0]),
    )


def test_split_vector_by_observable_type_returns_independent_slices():
    collection = FakeObservationCollection()
    split = helpers.split_vector_by_observable_type(collection, np.arange(7.0))

    np.testing.assert_array_equal(split[DOPPLER], np.array([0.0, 1.0, 2.0]))
    np.testing.assert_array_equal(split[PIXEL], np.array([3.0, 4.0, 5.0, 6.0]))
