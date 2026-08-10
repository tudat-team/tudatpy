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

_ADD_OBSERVATION_SET_SIGNATURES = (
    "Accepted signatures are "
    "add_observation_set(single_observation_set) and "
    "add_observation_set(observable_type, link_definition, observations, times, reference_link_end, "
    "dependent_variables=[], dependent_variable_bookkeeping=None, ancillary_settings=None, "
    "weights=[], residuals=[], sort_observations=False, erase_duplicate_observations=False)."
)

_ADD_OBSERVATION_SET_ARGUMENTS = (
    "observable_type",
    "link_definition",
    "observations",
    "times",
    "reference_link_end",
    "dependent_variables",
    "dependent_variable_bookkeeping",
    "ancillary_settings",
    "weights",
    "residuals",
    "sort_observations",
    "erase_duplicate_observations",
)

_ADD_OBSERVATION_SET_ALIASES = {
    "observations": ("observations", "observation_values"),
    "times": ("times", "observation_times"),
}

_MISSING_ADD_OBSERVATION_SET_ARGUMENT = object()


def _add_observation_set_type_error(detail):
    return TypeError(
        f"Invalid ObservationDataset.add_observation_set call: {detail}. "
        f"{_ADD_OBSERVATION_SET_SIGNATURES}"
    )


def _pop_add_observation_set_keyword(kwargs, canonical_name):
    aliases = _ADD_OBSERVATION_SET_ALIASES.get(canonical_name, (canonical_name,))
    present_aliases = [alias for alias in aliases if alias in kwargs]
    if len(present_aliases) > 1:
        raise _add_observation_set_type_error(
            f"received duplicate keyword aliases for '{canonical_name}': {', '.join(present_aliases)}"
        )
    if present_aliases:
        return kwargs.pop(present_aliases[0])
    return _MISSING_ADD_OBSERVATION_SET_ARGUMENT


def _normalized_component_add_observation_set_arguments(args, kwargs):
    if len(args) > len(_ADD_OBSERVATION_SET_ARGUMENTS):
        raise _add_observation_set_type_error(f"received {len(args)} positional arguments")

    keyword_arguments = dict(kwargs)
    values = {}
    provided = set()

    for argument_name, argument_value in zip(_ADD_OBSERVATION_SET_ARGUMENTS, args):
        values[argument_name] = argument_value
        provided.add(argument_name)

    for argument_name in _ADD_OBSERVATION_SET_ARGUMENTS:
        keyword_value = _pop_add_observation_set_keyword(keyword_arguments, argument_name)
        if keyword_value is _MISSING_ADD_OBSERVATION_SET_ARGUMENT:
            continue
        if argument_name in provided:
            raise _add_observation_set_type_error(f"received multiple values for '{argument_name}'")
        values[argument_name] = keyword_value
        provided.add(argument_name)

    if keyword_arguments:
        raise _add_observation_set_type_error(
            f"received unexpected keyword argument '{next(iter(keyword_arguments))}'"
        )

    missing_arguments = [
        argument_name
        for argument_name in _ADD_OBSERVATION_SET_ARGUMENTS[:5]
        if argument_name not in provided
    ]
    if missing_arguments:
        raise _add_observation_set_type_error(f"missing required argument '{missing_arguments[0]}'")

    return values, provided


def _add_single_observation_set_to_dataset(dataset, observation_set):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        single_set_dataset = create_observation_dataset_from_single_observation_set(observation_set)
    return dataset.add_observation_set_from_dataset(single_set_dataset, 0)


def _add_legacy_component_observation_set_to_dataset(dataset, values):
    # Only the 5 legacy required args plus ancillary_settings reach this path;
    # extended options are routed to the dataset-native component overload.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        observation_set = create_single_observation_set(
            values["observable_type"],
            values["link_definition"].link_ends,
            values["observations"],
            values["times"],
            values["reference_link_end"],
            values.get("ancillary_settings"),
        )
    return _add_single_observation_set_to_dataset(dataset, observation_set)


def _native_component_add_observation_set(dataset, values):
    defaulted_values = {
        "dependent_variables": [],
        "dependent_variable_bookkeeping": None,
        "ancillary_settings": None,
        "weights": [],
        "residuals": [],
        "sort_observations": False,
        "erase_duplicate_observations": False,
    }
    defaulted_values.update(values)
    return _native_add_observation_set_object_deprecation(
        dataset,
        defaulted_values["observable_type"],
        defaulted_values["link_definition"],
        defaulted_values["observations"],
        defaulted_values["times"],
        defaulted_values["reference_link_end"],
        defaulted_values["dependent_variables"],
        defaulted_values["dependent_variable_bookkeeping"],
        defaulted_values["ancillary_settings"],
        defaulted_values["weights"],
        defaulted_values["residuals"],
        defaulted_values["sort_observations"],
        defaulted_values["erase_duplicate_observations"],
    )


def _object_deprecation_add_observation_set_python_compatibility(self, *args, **kwargs):
    if len(args) == 1 and not kwargs:
        if not isinstance(args[0], SingleObservationSet):
            raise _add_observation_set_type_error(
                f"single-argument form expects SingleObservationSet, got {type(args[0]).__name__}"
            )
        return _add_single_observation_set_to_dataset(self, args[0])

    if not args and set(kwargs) <= {"single_observation_set", "observation_set"} and kwargs:
        if len(kwargs) != 1:
            raise _add_observation_set_type_error(
                "received both 'single_observation_set' and 'observation_set'"
            )
        observation_set = next(iter(kwargs.values()))
        if not isinstance(observation_set, SingleObservationSet):
            raise _add_observation_set_type_error(
                f"single-observation-set form expects SingleObservationSet, got {type(observation_set).__name__}"
            )
        return _add_single_observation_set_to_dataset(self, observation_set)

    component_values, component_provided = _normalized_component_add_observation_set_arguments(
        args, kwargs
    )
    if component_provided <= set(_ADD_OBSERVATION_SET_ARGUMENTS[:5]) | {"ancillary_settings"}:
        return _add_legacy_component_observation_set_to_dataset(self, component_values)
    return _native_component_add_observation_set(self, component_values)


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
