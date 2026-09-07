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


ObservationDataset.concatenated_times = _dataset_property_deprecation(
    "concatenated_times",
    "ObservationDataset.ordered_flattened_observation_data().times",
    lambda dataset: _legacy_collection(dataset).concatenated_times,
)
ObservationDataset.concatenated_times_objects = _dataset_property_deprecation(
    "concatenated_times_objects",
    "ObservationDataset.ordered_flattened_observation_data().times",
    lambda dataset: _legacy_collection(dataset).concatenated_times_objects,
)
ObservationDataset.concatenated_weights = _dataset_property_deprecation(
    "concatenated_weights",
    "ObservationDataset.ordered_flattened_observation_data().weight_vector",
    lambda dataset: _legacy_collection(dataset).concatenated_weights,
)
ObservationDataset.concatenated_observations = _dataset_property_deprecation(
    "concatenated_observations",
    "ObservationDataset.ordered_flattened_observation_data().observation_vector",
    lambda dataset: _legacy_collection(dataset).concatenated_observations,
)
ObservationDataset.concatenated_link_definition_ids = _dataset_property_deprecation(
    "concatenated_link_definition_ids",
    "ObservationDataset.ordered_flattened_observation_data().set_ids",
    lambda dataset: _legacy_collection(dataset).concatenated_link_definition_ids,
)
ObservationDataset.link_definition_ids = _dataset_property_deprecation(
    "link_definition_ids",
    "ObservationDataset.link_definition",
    lambda dataset: _legacy_collection(dataset).link_definition_ids,
)
ObservationDataset.observable_type_start_index_and_size = _dataset_property_deprecation(
    "observable_type_start_index_and_size",
    "ObservationDataset.ordered_flattened_observation_data()",
    lambda dataset: _legacy_collection(dataset).observable_type_start_index_and_size,
)
ObservationDataset.observation_set_start_index_and_size = _dataset_property_deprecation(
    "observation_set_start_index_and_size",
    "ObservationDataset.ordered_flattened_observation_data()",
    lambda dataset: _legacy_collection(dataset).observation_set_start_index_and_size,
)
ObservationDataset.observation_vector_size = _dataset_property_deprecation(
    "observation_vector_size",
    "ObservationDataset.total_scalar_size",
    lambda dataset: dataset.total_scalar_size,
)
ObservationDataset.sorted_observation_sets = _dataset_property_deprecation(
    "sorted_observation_sets",
    "ObservationDataset.observation_set_metadata",
    lambda dataset: _legacy_collection(dataset).sorted_observation_sets,
)
ObservationDataset.link_ends_per_observable_type = _dataset_property_deprecation(
    "link_ends_per_observable_type",
    "ObservationDataset.observation_set_metadata",
    lambda dataset: _legacy_collection(dataset).link_ends_per_observable_type,
)
ObservationDataset.link_definitions_per_observable = _dataset_property_deprecation(
    "link_definitions_per_observable",
    "ObservationDataset.observation_set_metadata",
    lambda dataset: _legacy_collection(dataset).link_definitions_per_observable,
)
ObservationDataset.time_bounds = _dataset_property_deprecation(
    "time_bounds",
    "ObservationDataset.ordered_flattened_observation_data().times",
    lambda dataset: _legacy_collection(dataset).time_bounds,
)
ObservationDataset.time_bounds_time_object = _dataset_property_deprecation(
    "time_bounds_time_object",
    "ObservationDataset.ordered_flattened_observation_data().times",
    lambda dataset: _legacy_collection(dataset).time_bounds_time_object,
)
ObservationDataset.sorted_per_set_time_bounds = _dataset_property_deprecation(
    "sorted_per_set_time_bounds",
    "ObservationDataset.observation_set_metadata",
    lambda dataset: _legacy_collection(dataset).sorted_per_set_time_bounds,
)

ObservationDataset.set_observations = _dataset_object_deprecation(
    "set_observations",
    "ObservationDataset.set_observations_for_set",
    lambda dataset, *args, **kwargs: _legacy_collection(dataset).set_observations(*args, **kwargs),
)
ObservationDataset.set_residuals = _dataset_object_deprecation(
    "set_residuals",
    "ObservationDataset.set_residuals_for_set",
    lambda dataset, *args, **kwargs: _legacy_collection(dataset).set_residuals(*args, **kwargs),
)
ObservationDataset.get_link_definitions_for_observables = _dataset_object_deprecation(
    "get_link_definitions_for_observables",
    "ObservationDataset.observation_set_metadata",
    lambda dataset, *args, **kwargs: _legacy_collection(
        dataset
    ).get_link_definitions_for_observables(*args, **kwargs),
)
ObservationDataset.get_single_link_and_type_observations = _dataset_object_deprecation(
    "get_single_link_and_type_observations",
    "ObservationDataset.create_viewer",
    lambda dataset, *args, **kwargs: _legacy_collection(
        dataset
    ).get_single_link_and_type_observations(*args, **kwargs),
)
ObservationDataset.get_observable_types = _dataset_object_deprecation(
    "get_observable_types",
    "ObservationDataset.observation_set_metadata",
    lambda dataset, *args, **kwargs: _legacy_collection(dataset).get_observable_types(
        *args, **kwargs
    ),
)
ObservationDataset.get_bodies_in_link_ends = _dataset_object_deprecation(
    "get_bodies_in_link_ends",
    "ObservationDataset.observation_set_metadata",
    lambda dataset, *args, **kwargs: _legacy_collection(dataset).get_bodies_in_link_ends(
        *args, **kwargs
    ),
)


# Explicit compatibility aliases. Collection membership operations belong on the
# legacy facade; automatically attaching them to Dataset could silently modify
# a temporary grouping instead of the dataset.
for _legacy_method_name in (
    "get_concatenated_observations",
    "get_concatenated_weights",
    "get_concatenated_residuals",
    "get_concatenated_computed_observations",
    "get_concatenated_observation_times",
    "get_concatenated_observation_times_objects",
    "get_concatenated_observations_and_times",
    "get_concatenated_observations_and_times_objects",
    "get_concatenated_link_definition_ids",
    "get_time_bounds_list",
    "set_constant_weight",
    "set_tabulated_weights",
):
    setattr(
        ObservationDataset,
        _legacy_method_name,
        _dataset_object_deprecation(
            _legacy_method_name,
            "ObservationDataset",
            lambda dataset, *args, _method=_legacy_method_name, **kwargs: getattr(
                _legacy_collection(dataset), _method
            )(*args, **kwargs),
        ),
    )

del _name, _object, _legacy_method_name
