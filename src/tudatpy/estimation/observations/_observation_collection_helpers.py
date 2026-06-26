"""Convenience helpers for mixed-observable observation collections.

The Tudat ``ObservationCollection`` stores observations by observable type and
link ends, while many estimation outputs expose one concatenated vector.  These
helpers keep observable-specific operations explicit when a collection contains
mixed units, for example Doppler residuals in Hz and pixel residuals in px.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from tudatpy.kernel.estimation.observations import (
    create_filtered_observation_collection as _create_filtered_observation_collection,
)
from tudatpy.estimation.observable_models_setup.model_settings import get_observable_size
from tudatpy.estimation.observations.observations_processing import (
    ObservationFilterType,
    observation_filter,
    observation_parser,
)


def _as_1d_array(values: Any, name: str) -> np.ndarray:
    """Return ``values`` as a one-dimensional float array."""
    array = np.asarray(values, dtype=float)
    if array.ndim == 0:
        return array.reshape(1)
    return array.reshape(-1)


def _as_filter_value(value: Any) -> float | np.ndarray:
    """Return a scalar or vector value suitable for Tudat filter overloads."""
    if np.isscalar(value):
        return float(value)
    return _as_1d_array(value, "filter value")


def _as_weight_value(value: Any) -> float | np.ndarray:
    """Return a scalar or vector value suitable for Tudat weight overloads."""
    if np.isscalar(value):
        return float(value)
    return _as_1d_array(value, "weight")


def _available_observable_types(observation_collection) -> list:
    return list(observation_collection.observable_type_start_index_and_size.keys())


def observable_type_parser(observable_type, use_opposite_condition: bool = False):
    """Create an observation parser selecting a single observable type."""
    return observation_parser(observable_type, use_opposite_condition)


def get_observable_type_slice(observation_collection, observable_type) -> slice:
    """Return the slice occupied by ``observable_type`` in full concatenated vectors.

    The returned slice applies to vectors with the same ordering as
    ``observation_collection.concatenated_observations`` and to estimation
    residual vectors that were produced from this collection.
    """
    start_and_size_by_type = observation_collection.observable_type_start_index_and_size
    if observable_type not in start_and_size_by_type:
        raise KeyError(
            f"Observable type {observable_type!r} is not present in the observation collection. "
            f"Available observable types: {_available_observable_types(observation_collection)}"
        )
    start_and_size = start_and_size_by_type[observable_type]
    start = int(start_and_size[0])
    size = int(start_and_size[1])
    return slice(start, start + size)


def extract_observable_type_vector(observation_collection, observable_type, vector) -> np.ndarray:
    """Extract one observable type from a full concatenated vector."""
    values = np.asarray(vector)
    observable_slice = get_observable_type_slice(observation_collection, observable_type)
    if values.ndim == 1:
        if observable_slice.stop > values.shape[0]:
            raise ValueError(
                "Input vector is shorter than the requested observable-type slice: "
                f"stop={observable_slice.stop}, length={values.shape[0]}"
            )
        return values[observable_slice]
    if values.ndim == 2:
        if observable_slice.stop > values.shape[0]:
            raise ValueError(
                "Input matrix has too few rows for the requested observable-type slice: "
                f"stop={observable_slice.stop}, rows={values.shape[0]}"
            )
        return values[observable_slice, :]
    raise ValueError("Input vector must be one- or two-dimensional")


def get_observable_type_observations(
    observation_collection,
    observable_type,
    observation_vector: Any | None = None,
) -> np.ndarray:
    """Return observations for one observable type.

    If ``observation_vector`` is supplied, it is interpreted as a full
    concatenated vector and sliced using collection metadata.  Otherwise the
    values are retrieved from the collection through an observable-type parser.
    """
    if observation_vector is not None:
        return extract_observable_type_vector(
            observation_collection, observable_type, observation_vector
        )
    parser = observable_type_parser(observable_type)
    return np.asarray(observation_collection.get_concatenated_observations(parser))


def get_observable_type_residuals(
    observation_collection,
    observable_type,
    residual_vector: Any | None = None,
) -> np.ndarray:
    """Return residuals for one observable type.

    ``residual_vector`` may be a full residual-history column from an estimation
    output.  If omitted, residuals stored in the collection are returned.
    """
    if residual_vector is not None:
        return extract_observable_type_vector(
            observation_collection, observable_type, residual_vector
        )
    parser = observable_type_parser(observable_type)
    return np.asarray(observation_collection.get_concatenated_residuals(parser))


def get_observable_type_residual_history(
    observation_collection,
    observable_type,
    residual_history,
) -> np.ndarray:
    """Return all residual-history rows for one observable type."""
    return extract_observable_type_vector(observation_collection, observable_type, residual_history)


def get_observable_type_weights(
    observation_collection,
    observable_type,
    weights_vector: Any | None = None,
) -> np.ndarray:
    """Return weights for one observable type."""
    if weights_vector is not None:
        return extract_observable_type_vector(
            observation_collection, observable_type, weights_vector
        )
    parser = observable_type_parser(observable_type)
    return np.asarray(observation_collection.get_concatenated_weights(parser))


def get_observable_type_times(
    observation_collection,
    observable_type,
    times_vector: Any | None = None,
) -> np.ndarray:
    """Return scalar-aligned concatenated times for one observable type.

    The returned vector has the same length and ordering as the scalar entries
    in ``observation_collection.concatenated_observations``.  For vector
    observables, the observation epoch is repeated once per component.
    """
    if times_vector is not None:
        return extract_observable_type_vector(observation_collection, observable_type, times_vector)
    return extract_observable_type_vector(
        observation_collection, observable_type, observation_collection.concatenated_times
    ).astype(float)


def get_observable_type_observation_times(
    observation_collection,
    observable_type,
) -> np.ndarray:
    """Return one observation epoch per scalar or vector observation."""
    parser = observable_type_parser(observable_type)
    return np.asarray(
        observation_collection.get_concatenated_observation_times(parser), dtype=float
    )


def reshape_observable_components(
    values, observable_type, observable_size: int | None = None
) -> np.ndarray:
    """Reshape a flat observable vector to ``(n_observations, observable_size)``."""
    if observable_size is None:
        observable_size = int(get_observable_size(observable_type))
    if observable_size <= 0:
        raise ValueError(f"Observable size must be positive, got {observable_size}")
    flat_values = _as_1d_array(values, "observable values")
    if flat_values.size % observable_size != 0:
        raise ValueError(
            f"Cannot reshape {flat_values.size} values into observable size {observable_size}"
        )
    return flat_values.reshape((-1, observable_size))


def get_observable_type_residual_matrix(
    observation_collection,
    observable_type,
    residual_vector: Any | None = None,
) -> np.ndarray:
    """Return residuals as one row per observation and one column per component."""
    residuals = get_observable_type_residuals(
        observation_collection, observable_type, residual_vector
    )
    return reshape_observable_components(residuals, observable_type)


def get_observable_type_residual_norms(
    observation_collection,
    observable_type,
    residual_vector: Any | None = None,
    norm_order: int | float | None = None,
) -> np.ndarray:
    """Return vector-norm residuals for one observable type."""
    residual_matrix = get_observable_type_residual_matrix(
        observation_collection, observable_type, residual_vector
    )
    return np.linalg.norm(residual_matrix, ord=norm_order, axis=1)


def get_observable_type_times_per_observation(
    observation_collection,
    observable_type,
    times_vector: Any | None = None,
    validate_component_times: bool = True,
) -> np.ndarray:
    """Return one time per observation for scalar or vector observables."""
    if times_vector is None:
        return get_observable_type_observation_times(observation_collection, observable_type)

    times = get_observable_type_times(observation_collection, observable_type, times_vector)
    observable_size = int(get_observable_size(observable_type))
    if observable_size <= 1:
        return _as_1d_array(times, "observation times")

    flat_times = _as_1d_array(times, "observation times")
    if flat_times.size % observable_size != 0:
        raise ValueError(
            f"Cannot group {flat_times.size} time entries into observable size {observable_size}"
        )
    time_matrix = flat_times.reshape((-1, observable_size))
    if validate_component_times and not np.all(time_matrix == time_matrix[:, [0]]):
        raise ValueError("Component times for at least one vector observation are not identical")
    return time_matrix[:, 0]


def split_vector_by_observable_type(observation_collection, vector) -> dict:
    """Split a full concatenated vector into a dictionary by observable type."""
    return {
        observable_type: extract_observable_type_vector(
            observation_collection, observable_type, vector
        )
        for observable_type in _available_observable_types(observation_collection)
    }


def set_constant_weight_by_type(observation_collection, observable_type, weight):
    """Set a constant scalar or component-wise vector weight for one observable type."""
    parser = observable_type_parser(observable_type)
    observation_collection.set_constant_weight(_as_weight_value(weight), parser)
    return observation_collection


def set_constant_weight_from_sigma_by_type(observation_collection, observable_type, sigma):
    """Set constant weights ``1 / sigma**2`` for one observable type."""
    sigma_array = _as_weight_value(sigma)
    if np.isscalar(sigma_array):
        if sigma_array <= 0.0:
            raise ValueError("sigma must be positive")
        weight = 1.0 / sigma_array**2
    else:
        if np.any(sigma_array <= 0.0):
            raise ValueError("all sigma values must be positive")
        weight = 1.0 / sigma_array**2
    return set_constant_weight_by_type(observation_collection, observable_type, weight)


def filter_observations_by_type(
    observation_collection,
    observable_type,
    filter_settings,
    save_filtered_observations: bool = True,
):
    """Apply an observation filter in-place to one observable type."""
    parser = observable_type_parser(observable_type)
    observation_collection.filter_observations(filter_settings, parser, save_filtered_observations)
    return observation_collection


def create_filtered_observation_collection_by_type(
    observation_collection,
    observable_type,
    filter_settings,
):
    """Create a filtered copy by applying a filter to one observable type."""
    parser = observable_type_parser(observable_type)
    return _create_filtered_observation_collection(observation_collection, filter_settings, parser)


def filter_residuals_by_type(
    observation_collection,
    observable_type,
    residual_threshold,
    filter_out: bool = True,
    use_opposite_condition: bool = False,
    save_filtered_observations: bool = True,
):
    """Apply residual filtering in-place to one observable type only."""
    filter_settings = observation_filter(
        ObservationFilterType.residual_filtering,
        _as_filter_value(residual_threshold),
        filter_out,
        use_opposite_condition,
    )
    return filter_observations_by_type(
        observation_collection,
        observable_type,
        filter_settings,
        save_filtered_observations,
    )


def create_residual_filtered_observation_collection_by_type(
    observation_collection,
    observable_type,
    residual_threshold,
    filter_out: bool = True,
    use_opposite_condition: bool = False,
):
    """Create a filtered copy using residual filtering on one observable type only."""
    filter_settings = observation_filter(
        ObservationFilterType.residual_filtering,
        _as_filter_value(residual_threshold),
        filter_out,
        use_opposite_condition,
    )
    return create_filtered_observation_collection_by_type(
        observation_collection,
        observable_type,
        filter_settings,
    )


__all__ = [
    "observable_type_parser",
    "get_observable_type_slice",
    "extract_observable_type_vector",
    "get_observable_type_observations",
    "get_observable_type_residuals",
    "get_observable_type_residual_history",
    "get_observable_type_weights",
    "get_observable_type_times",
    "get_observable_type_observation_times",
    "reshape_observable_components",
    "get_observable_type_residual_matrix",
    "get_observable_type_residual_norms",
    "get_observable_type_times_per_observation",
    "split_vector_by_observable_type",
    "set_constant_weight_by_type",
    "set_constant_weight_from_sigma_by_type",
    "filter_observations_by_type",
    "create_filtered_observation_collection_by_type",
    "filter_residuals_by_type",
    "create_residual_filtered_observation_collection_by_type",
]
