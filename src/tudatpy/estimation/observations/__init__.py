from functools import wraps

from tudatpy._deprecation import deprecation_warning, property_deprecation
from tudatpy.estimation.observable_models_setup.links import (
    LinkDefinition,
    LinkEndId,
    LinkEndType,
    observed_body,
    observer,
    receiver,
    transmitter,
)
from tudatpy.estimation.observable_models_setup.model_settings import (
    ObservableType,
    angular_position_type as angular_position,
    angular_position_type,
    azimuth_elevation_type as azimuth_elevation,
    azimuth_elevation_type,
    one_way_instantaneous_doppler_type as one_way_doppler,
    one_way_instantaneous_doppler_type,
    one_way_range_type as one_way_range,
    one_way_range_type,
)
from tudatpy.kernel.estimation.observations import *

from ._query import observation_query

for _name, _object in list(globals().items()):
    if getattr(_object, "__module__", None) == "tudatpy.kernel.estimation.observations":
        _object.__module__ = "tudatpy.estimation.observations"

SingleObservationSet.ancilliary_settings = property_deprecation(
    "SingleObservationSet.ancilliary_settings", "SingleObservationSet.ancillary_settings"
)(SingleObservationSet.ancillary_settings)

_native_add_observation_set = ObservationDataset.add_observation_set


@wraps(_native_add_observation_set)
def _add_observation_set(self, *args, **kwargs):
    # Keep previously used keyword spellings; pybind validates the actual signatures.
    for alias, canonical in (
        ("observation_values", "observations"),
        ("observation_times", "times"),
        ("observation_set", "single_observation_set"),
    ):
        if alias in kwargs:
            if canonical in kwargs:
                raise TypeError(f"Both {alias!r} and {canonical!r} were supplied")
            kwargs[canonical] = kwargs.pop(alias)
    return _native_add_observation_set(self, *args, **kwargs)


ObservationDataset.add_observation_set = _add_observation_set


def _legacy_collection(dataset):
    return create_observation_collection_from_dataset(dataset)


def _dataset_property_deprecation(old_name, new_name, getter):
    def wrapped(self):
        deprecation_warning(f"ObservationDataset.{old_name}", new_name)
        return getter(self)

    return property(wrapped)


def _dataset_object_deprecation(old_name, new_name, method):
    def wrapped(self, *args, **kwargs):
        deprecation_warning(f"ObservationDataset.{old_name}", new_name)
        return method(self, *args, **kwargs)

    wrapped.__name__ = old_name
    return wrapped


# Keep the legacy aliases explicit: membership changes must not act on a temporary collection.
for _legacy_name, _replacement in {
    "concatenated_times": "ObservationDataset.get_times",
    "concatenated_times_objects": "ObservationDataset.get_times",
    "concatenated_weights": "ObservationDataset.get_weight_diagonal",
    "concatenated_observations": "ObservationDataset.get_observations",
    "concatenated_link_definition_ids": "ObservationDataset.get_set_ids",
    "link_definition_ids": "ObservationDataset.link_definition",
    "observable_type_start_index_and_size": "ObservationDataset.get_data",
    "observation_set_start_index_and_size": "ObservationDataset.get_data",
    "sorted_observation_sets": "ObservationDataset.observation_set_metadata",
    "link_ends_per_observable_type": "ObservationDataset.observation_set_metadata",
    "link_definitions_per_observable": "ObservationDataset.observation_set_metadata",
    "time_bounds": "ObservationDataset.get_times",
    "time_bounds_time_object": "ObservationDataset.get_times",
    "sorted_per_set_time_bounds": "ObservationDataset.observation_set_metadata",
}.items():
    setattr(
        ObservationDataset,
        _legacy_name,
        _dataset_property_deprecation(
            _legacy_name,
            _replacement,
            lambda dataset, _member=_legacy_name: getattr(_legacy_collection(dataset), _member),
        ),
    )

ObservationDataset.observation_vector_size = _dataset_property_deprecation(
    "observation_vector_size",
    "ObservationDataset.total_scalar_size",
    lambda dataset: dataset.total_scalar_size,
)

for _legacy_name, _replacement in {
    "set_observations": "ObservationDataset.set_observations_for_set",
    "set_residuals": "ObservationDataset.set_residuals_for_set",
    "get_link_definitions_for_observables": "ObservationDataset.observation_set_metadata",
    "get_single_link_and_type_observations": "ObservationDataset.get_data",
    "get_observable_types": "ObservationDataset.observation_set_metadata",
    "get_bodies_in_link_ends": "ObservationDataset.observation_set_metadata",
    "get_concatenated_observations": "ObservationDataset",
    "get_concatenated_weights": "ObservationDataset",
    "get_concatenated_residuals": "ObservationDataset",
    "get_concatenated_computed_observations": "ObservationDataset",
    "get_concatenated_observation_times": "ObservationDataset",
    "get_concatenated_observation_times_objects": "ObservationDataset",
    "get_concatenated_observations_and_times": "ObservationDataset",
    "get_concatenated_observations_and_times_objects": "ObservationDataset",
    "get_concatenated_link_definition_ids": "ObservationDataset",
    "get_time_bounds_list": "ObservationDataset",
    "set_constant_weight": "ObservationDataset",
    "set_tabulated_weights": "ObservationDataset",
}.items():
    setattr(
        ObservationDataset,
        _legacy_name,
        _dataset_object_deprecation(
            _legacy_name,
            _replacement,
            lambda dataset, *args, _member=_legacy_name, **kwargs: getattr(
                _legacy_collection(dataset), _member
            )(*args, **kwargs),
        ),
    )

del _name, _object, _legacy_name, _replacement
