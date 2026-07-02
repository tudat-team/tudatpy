import warnings

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

_native_add_observation_set_object_deprecation = ObservationDataset.add_observation_set


def _object_deprecation_add_observation_set_python_compatibility(self, *args, **kwargs):
    try:
        return _native_add_observation_set_object_deprecation(self, *args, **kwargs)
    except TypeError:
        if len(args) < 5:
            raise

        observable_type = args[0]
        link_definition = args[1]
        observation_values = args[2]
        observation_times = args[3]
        reference_link_end = args[4]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            observation_set = create_single_observation_set(
                observable_type,
                link_definition.link_ends,
                observation_values,
                observation_times,
                reference_link_end,
            )
            single_set_dataset = create_observation_dataset_from_single_observation_set(
                observation_set
            )

        return self.add_observation_set_from_dataset(single_set_dataset, 0)


ObservationDataset.add_observation_set = (
    _object_deprecation_add_observation_set_python_compatibility
)


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


def _observation_collection_property(name):
    return _dataset_property_deprecation(
        name,
        "ObservationDataset",
        lambda dataset, collection_name=name: getattr(_legacy_collection(dataset), collection_name),
    )


def _observation_collection_method(name):
    return _dataset_object_deprecation(
        name,
        "ObservationDataset",
        lambda dataset, *args, collection_name=name, **kwargs: getattr(
            _legacy_collection(dataset), collection_name
        )(*args, **kwargs),
    )


for _name in dir(ObservationCollection):
    if _name.startswith("_") or hasattr(ObservationDataset, _name):
        continue

    _object = getattr(ObservationCollection, _name)
    if isinstance(_object, property):
        setattr(ObservationDataset, _name, _observation_collection_property(_name))
    elif callable(_object):
        setattr(ObservationDataset, _name, _observation_collection_method(_name))

del _name, _object
