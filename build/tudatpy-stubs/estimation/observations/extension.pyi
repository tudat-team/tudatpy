import numpy
import pybind11_stubgen.typing_ext
from ...astro import time_representation
from ...dynamics import environment
from ...estimation.observable_models import observables_simulation
from ...estimation.observable_models_setup import links
from ...estimation.observable_models_setup import model_settings
from ...estimation.observations_setup import ancillary_settings
from ...estimation.observations_setup import observations_dependent_variables
import typing
from . import observations_geometry
from . import observations_processing
__all__ = ['ObservationCollection', 'SingleObservationSet', 'compute_residuals_and_dependent_variables', 'create_filtered_observation_collection', 'create_new_observation_collection', 'create_single_observation_set', 'filter_observations', 'merge_observation_collections', 'observations_geometry', 'observations_processing', 'single_observation_set', 'split_observation_collection', 'split_observation_set']

class ObservationCollection:
    """Class collecting all observations and associated data for use in an estimation.
    
    Class containing the full set of observations and associated data, typically for input into the estimation. When using simulated data,
    this class is instantiated via a call to the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. More information is provided
    on the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_"""

    def __init__(self, observation_sets: list[SingleObservationSet]) -> None:
        ...

    def add_dependent_variable(self, dependent_variable_settings: observations_dependent_variables.ObservationDependentVariableSettings, bodies: environment.SystemOfBodies, observation_parser: observations_processing.ObservationCollectionParser=...) -> observations_processing.ObservationCollectionParser:
        """
        No documentation found.
        """

    def append(self, observation_collection_to_append: ObservationCollection) -> None:
        ...

    def compatible_dependent_variable_settings(self, dependent_variable_settings: observations_dependent_variables.ObservationDependentVariableSettings, observation_parser: observations_processing.ObservationCollectionParser=...) -> tuple[list[list[observations_dependent_variables.ObservationDependentVariableSettings]], observations_processing.ObservationCollectionParser]:
        """
        No documentation found.
        """

    def compatible_dependent_variables_list(self, dependent_variable_settings: observations_dependent_variables.ObservationDependentVariableSettings, observation_parser: observations_processing.ObservationCollectionParser=...) -> tuple[list[list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]], observations_processing.ObservationCollectionParser]:
        """
        No documentation found.
        """

    def concatenated_dependent_variable(self, dependent_variable_settings: observations_dependent_variables.ObservationDependentVariableSettings, first_compatible_settings: bool=False, observation_parser: observations_processing.ObservationCollectionParser=...) -> tuple[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], observations_processing.ObservationCollectionParser]:
        """
        No documentation found.
        """

    def dependent_variable(self, dependent_variable_settings: observations_dependent_variables.ObservationDependentVariableSettings, first_compatible_settings: bool=False, observation_parser: observations_processing.ObservationCollectionParser=...) -> tuple[list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]], observations_processing.ObservationCollectionParser]:
        """
        No documentation found.
        """

    def dependent_variable_history(self, dependent_variable_settings: observations_dependent_variables.ObservationDependentVariableSettings, first_compatible_settings: bool=False, observation_parser: observations_processing.ObservationCollectionParser=...) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
        No documentation found.
        """

    def dependent_variable_history_per_set(self, dependent_variable_settings: observations_dependent_variables.ObservationDependentVariableSettings, first_compatible_settings: bool=False, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]]:
        """
        No documentation found.
        """

    def dependent_variable_history_per_set_time_object(self, dependent_variable_settings: observations_dependent_variables.ObservationDependentVariableSettings, first_compatible_settings: bool=False, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]]:
        """
        No documentation found.
        """

    def dependent_variable_history_time_object(self, dependent_variable_settings: observations_dependent_variables.ObservationDependentVariableSettings, first_compatible_settings: bool=False, observation_parser: observations_processing.ObservationCollectionParser=...) -> dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
        No documentation found.
        """

    @typing.overload
    def filter_observations(self, observation_filters: dict[observations_processing.ObservationCollectionParser, observations_processing.ObservationFilterBase], save_filtered_observations: bool=True) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def filter_observations(self, observation_filters: observations_processing.ObservationFilterBase, observation_parser: observations_processing.ObservationCollectionParser=..., save_filtered_observations: bool=True) -> None:
        """
        No documentation found.
        """

    def get_bodies_in_link_ends(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[str]:
        """
        No documentation found.
        """

    def get_computed_observations(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
        No documentation found.
        """

    def get_concatenated_computed_observations(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
        No documentation found.
        """

    def get_concatenated_link_definition_ids(self, observation_parser: observations_processing.ObservationCollectionParser) -> list[int]:
        """
        No documentation found.
        """

    def get_concatenated_observation_times(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[float]:
        """
        No documentation found.
        """

    def get_concatenated_observation_times_objects(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[time_representation.Time]:
        """
        No documentation found.
        """

    def get_concatenated_observations(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
        No documentation found.
        """

    def get_concatenated_observations_and_times(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> tuple[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)], list[float]]:
        """
        No documentation found.
        """

    def get_concatenated_observations_and_times_objects(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> tuple[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)], list[time_representation.Time]]:
        """
        No documentation found.
        """

    def get_concatenated_residuals(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
        No documentation found.
        """

    def get_concatenated_weights(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
        No documentation found.
        """

    def get_link_definitions_for_observables(self, observable_type: model_settings.ObservableType) -> list[links.LinkDefinition]:
        """
        No documentation found.
        """

    def get_mean_residuals(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
        No documentation found.
        """

    def get_observable_types(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[model_settings.ObservableType]:
        """
        No documentation found.
        """

    def get_observation_times(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[list[float]]:
        """
        No documentation found.
        """

    def get_observation_times_objects(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[list[time_representation.Time]]:
        """
        No documentation found.
        """

    def get_observations(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
        No documentation found.
        """

    def get_observations_and_times(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> tuple[list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], list[list[float]]]:
        """
        No documentation found.
        """

    def get_observations_and_times_objects(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> tuple[list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], list[list[time_representation.Time]]]:
        """
        No documentation found.
        """

    def get_reference_points_in_link_ends(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[str]:
        """
        No documentation found.
        """

    def get_residuals(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
        No documentation found.
        """

    def get_rms_residuals(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
        No documentation found.
        """

    def get_single_link_and_type_observations(self, observable_type: model_settings.ObservableType, link_definition: links.LinkDefinition) -> list[SingleObservationSet]:
        """
                 Function to get all observation sets for a given observable type and link definition.
        
        
                 Parameters
                 ----------
                 observable_type : :class:`ObservableType`
                     Observable type of which observations are to be simulated.
                 link_ends : LinkDefinition
                     Link ends for which observations are to be simulated.
                 Returns
                 -------
                 list[ SingleObservationSet ]
                     List of observation sets for given observable type and link definition.
        """

    def get_single_observation_sets(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[SingleObservationSet]:
        """
        No documentation found.
        """

    def get_time_bounds_list(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[tuple[float, float]]:
        """
        No documentation found.
        """

    def get_time_bounds_list_time_object(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[tuple[time_representation.Time, time_representation.Time]]:
        """
        No documentation found.
        """

    def get_time_bounds_per_set(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[tuple[float, float]]:
        """
        No documentation found.
        """

    def get_time_bounds_per_set_time_object(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[tuple[time_representation.Time, time_representation.Time]]:
        """
        No documentation found.
        """

    def get_weights(self, observation_parser: observations_processing.ObservationCollectionParser=...) -> list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
        No documentation found.
        """

    def print_observation_sets_start_and_size(self) -> None:
        """
        No documentation found.
        """

    def remove_empty_observation_sets(self) -> None:
        """
        No documentation found.
        """

    def remove_single_observation_sets(self, observation_parser: observations_processing.ObservationCollectionParser) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_constant_weight(self, weight: float, observation_parser: observations_processing.ObservationCollectionParser=...) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_constant_weight(self, weight: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)], observation_parser: observations_processing.ObservationCollectionParser=...) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_constant_weight_per_observation_parser(self, weights_per_observation_parser: dict[observations_processing.ObservationCollectionParser, float]) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_constant_weight_per_observation_parser(self, weights_per_observation_parser: dict[observations_processing.ObservationCollectionParser, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_observations(self, observations: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]) -> None:
        """
        Function to reset the full list of observable values. The order of the observations must be the same as for :attr:`~ObservationCollection.concatenated_observations`
        
        Parameters
        ----------
        observations : np.ndarray
            New list of observable values
        """

    @typing.overload
    def set_observations(self, observations: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)], observation_parser: observations_processing.ObservationCollectionParser) -> None:
        """
        Function to reset a subset of all observable values, with this subset defined by the ``observation_parser`` input.
        The order and size of the new observation vector must be the same as when calling :attr:`~ObservationCollection.concatenated_observations` on
        an ``ObservationCollection`` containing only the parsed observation.
        
        Parameters
        ----------
        observations : np.ndarray
            New list of observable values
        observation_parser : ObservationCollectionParser
            Observation parser with which to select the subset of observations that is to be reset
        """

    @typing.overload
    def set_observations(self, observations_per_parser: dict[observations_processing.ObservationCollectionParser, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]) -> None:
        """
        Function to reset a subset of all observable values, with this subset defined by a list of observation parsers input.
        Each observation parser is associated with a new set of observable values.
        The order and size of the new observation vector for each parser must be the same as when calling :attr:`~ObservationCollection.concatenated_observations` on
        an ``ObservationCollection`` containing only the parsed observation (from a single parser). NOTE: since the multiple parsers
        are handled in order (iterating over the keys of ``observations_per_parser``) some observations may be reset several times,
        in case.
        
        Parameters
        ----------
        observations : np.ndarray
            New list of observable values
        observation_parser : ObservationCollectionParser
            Observation parser with which to select the subset of observations that is to be reset
        """

    @typing.overload
    def set_reference_point(self, bodies: environment.SystemOfBodies, antenna_position: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)], antenna_name: str, spacecraft_name: str, link_end_type: links.LinkEndType, observation_parser: observations_processing.ObservationCollectionParser=...) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_reference_point(self, bodies: environment.SystemOfBodies, antenna_body_fixed_ephemeris: environment.Ephemeris, antenna_name: str, spacecraft_name: str, link_end_type: links.LinkEndType, observation_parser: observations_processing.ObservationCollectionParser=...) -> None:
        """
        No documentation found.
        """

    def set_reference_points(self, bodies: environment.SystemOfBodies, antenna_switch_history: dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]], spacecraft_name: str, link_end_type: links.LinkEndType, observation_parser: observations_processing.ObservationCollectionParser=...) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_residuals(self, residuals: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_residuals(self, residuals: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)], observation_parser: observations_processing.ObservationCollectionParser) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_residuals(self, residuals_per_parser: dict[observations_processing.ObservationCollectionParser, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_tabulated_weights(self, tabulated_weights: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)], observation_parser: observations_processing.ObservationCollectionParser=...) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_tabulated_weights(self, tabulated_weights: dict[observations_processing.ObservationCollectionParser, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]) -> None:
        """
        No documentation found.
        """

    def set_transponder_delay(self, spacecraft_name: str, transponder_delay: float, observation_parser: observations_processing.ObservationCollectionParser=...) -> None:
        """
        No documentation found.
        """

    def split_observation_sets(self, observation_set_splitter: observations_processing.ObservationSetSplitterBase, observation_parser: observations_processing.ObservationCollectionParser=...) -> None:
        """
        No documentation found.
        """

    @property
    def concatenated_link_definition_ids(self) -> list[int]:
        """
                 **read-only**
        
                 Vector containing concatenated indices identifying the link ends. Each set of link ends is assigned a unique integer identifier (for a given instance of this class). The definition of a given integer identifier with the link ends is given by this class' :func:`link_definition_ids` function. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order of the present vector.
        
                 :type: numpy.ndarray[ int ]
        """

    @property
    def concatenated_observations(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
                 **read-only**
        
                 Vector containing concatenated observable values. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order
        
                 :type: numpy.ndarray[numpy.float64[m, 1]]
        """

    @property
    def concatenated_times(self) -> list[float]:
        """
                 **read-only**
        
                 Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order
        
                 :type: numpy.ndarray[numpy.float64[m, 1]]
        """

    @property
    def concatenated_times_objects(self) -> list[time_representation.Time]:
        """
                 **read-only**
        
                 Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order
        
                 :type: numpy.ndarray[numpy.float64[m, 1]]
        """

    @property
    def concatenated_weights(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
        No documentation found.
        """

    @property
    def link_definition_ids(self) -> dict[int, dict[links.LinkEndType, links.LinkEndId]]:
        """
                 **read-only**
        
                 Dictionaty mapping a link end integer identifier to the specific link ends
        
                 :type: dict[ int, dict[ LinkEndType, LinkEndId ] ]
        """

    @property
    def link_definitions_per_observable(self) -> dict[model_settings.ObservableType, list[links.LinkDefinition]]:
        """
        No documentation found.
        """

    @property
    def link_ends_per_observable_type(self) -> dict[model_settings.ObservableType, list[dict[links.LinkEndType, links.LinkEndId]]]:
        """
        No documentation found.
        """

    @property
    def observable_type_start_index_and_size(self) -> dict[model_settings.ObservableType, tuple[int, int]]:
        """
                 **read-only**
        
                 Dictionary defining per obervable type (dict key), the index in the full observation vector (:func:`concatenated_observations`) where the given observable type starts, and the number of subsequent entries in this vector containing a value of an observable of this type
        
                 :type: dict[ ObservableType, [ int, int ] ]
        """

    @property
    def observation_set_start_index_and_size(self) -> dict[model_settings.ObservableType, dict[int, list[tuple[int, int]]]]:
        """
                 **read-only**
        
                 The nested dictionary/list returned by this property mirrors the structure of the :func:`sorted_observation_sets` property of this class. The present function provides the start index and size of the observables in the full observation vector that come from the correspoding `SingleObservationSet` in the :func:`sorted_observation_sets` Consequently, the present property returns a nested dictionary defining per obervable type, link end identifier, and `SingleObservationSet` index (for the given observable type and link end identifier), where the observables in the given `SingleObservationSet` starts, and the number of subsequent entries in this vector containing data from it.
        
                 :type: dict[ ObservableType, dict[ int, list[ int, int ] ] ]
        """

    @property
    def observation_vector_size(self) -> int:
        """
                 **read-only**
        
                 Length of the total vector of observations
        
                 :type: int
        """

    @property
    def sorted_observation_sets(self) -> dict[model_settings.ObservableType, dict[int, list[SingleObservationSet]]]:
        """
                 **read-only**
        
                 The nested dictionary/list contains the list of `SingleObservationSet` objects, in the same method as they are stored internally in the present class. Specifics on the storage order are given in the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_
        
                 :type: dict[ ObservableType, dict[ int, list[ SingleObservationSet ] ] ]
        """

    @property
    def sorted_per_set_time_bounds(self) -> dict[model_settings.ObservableType, dict[int, list[tuple[float, float]]]]:
        """
        No documentation found.
        """

    @property
    def time_bounds(self) -> tuple[float, float]:
        """
        No documentation found.
        """

    @property
    def time_bounds_time_object(self) -> tuple[time_representation.Time, time_representation.Time]:
        """
        No documentation found.
        """

class SingleObservationSet:
    """Class collecting a single set of observations and associated data, of a given observable type, link ends, and ancilliary data."""

    def compatible_dependent_variable_settings(self, arg0: observations_dependent_variables.ObservationDependentVariableSettings) -> list[observations_dependent_variables.ObservationDependentVariableSettings]:
        """
        No documentation found.
        """

    def compatible_dependent_variables_list(self, arg0: observations_dependent_variables.ObservationDependentVariableSettings) -> list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]:
        """
        No documentation found.
        """

    def filter_observations(self, filter: observations_processing.ObservationFilterBase, save_filtered_obs: bool=True) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_constant_weight(self, weight: float) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_constant_weight(self, weight: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_observations(self, observations: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_observations(self, observations: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_residuals(self, residuals: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]) -> None:
        """
        No documentation found.
        """

    @typing.overload
    def set_residuals(self, residuals: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]) -> None:
        """
        No documentation found.
        """

    def set_tabulated_weights(self, weights: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]) -> None:
        """
        No documentation found.
        """

    def single_dependent_variable(self, dependent_variable_settings: observations_dependent_variables.ObservationDependentVariableSettings, return_first_compatible_settings: bool=False) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
        No documentation found.
        """

    def single_dependent_variable_history(self, arg0: observations_dependent_variables.ObservationDependentVariableSettings, arg1: bool) -> dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
        No documentation found.
        """

    @property
    def ancilliary_settings(self) -> ancillary_settings.ObservationAncilliarySimulationSettings:
        """
                 **read-only**
        
                 Ancilliary settings all stored observations
        
                 :type: ObservationAncilliarySimulationSettings
        """

    @property
    def computed_observations(self) -> list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
        No documentation found.
        """

    @property
    def concatenad_weights(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
        No documentation found.
        """

    @property
    def concatenated_computed_observations(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
        No documentation found.
        """

    @property
    def concatenated_observations(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
                 **read-only**
        
                 Concatenated vector of all stored observations
        
                 :type: numpy.ndarray[numpy.float64[m, 1]]
        """

    @property
    def concatenated_residuals(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
        No documentation found.
        """

    @property
    def dependent_variables(self) -> list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
        No documentation found.
        """

    @dependent_variables.setter
    def dependent_variables(self, arg1: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]) -> None:
        ...

    @property
    def dependent_variables_history(self) -> dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
        No documentation found.
        """

    @property
    def dependent_variables_matrix(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
        No documentation found.
        """

    @property
    def filtered_observation_set(self) -> SingleObservationSet:
        """
        No documentation found.
        """

    @property
    def link_definition(self) -> links.LinkDefinition:
        """
                 **read-only**
        
                 Definition of the link ends for which the object stores observations
        
                 :type: LinkDefinition
        """

    @link_definition.setter
    def link_definition(self, arg1: links.LinkDefinition) -> None:
        ...

    @property
    def list_of_observations(self) -> list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 List of separate stored observations. Each entry of this list is a vector containing a single observation. In cases where the observation is single-valued (range, Doppler), the vector is size 1, but for multi-valued observations such as angular position, each vector in the list will have size >1
        
                 :type: list[ numpy.ndarray[numpy.float64[m, 1]] ]
        """

    @property
    def mean_residuals(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
        No documentation found.
        """

    @property
    def number_filtered_observations(self) -> int:
        """
        No documentation found.
        """

    @property
    def number_of_observables(self) -> int:
        """
        No documentation found.
        """

    @property
    def observable_type(self) -> model_settings.ObservableType:
        """
                 **read-only**
        
                 Type of observable for which the object stores observations
        
                 :type: ObservableType
        """

    @property
    def observation_times(self) -> list[time_representation.Time]:
        """
                 **read-only**
        
                 Reference time for each of the observations in ``list_of_observations``
        
                 :type: list[ float]
        """

    @property
    def observations_history(self) -> dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 Dictionary of observations sorted by time. Created by making a dictionary with ``observation_times`` as keys and  ``list_of_observations`` as values
        
                 :type: dict[ float, numpy.ndarray[numpy.float64[m, 1]] ]
        """

    @property
    def reference_link_end(self) -> links.LinkEndType:
        """
                 **read-only**
        
                 Reference link end for all stored observations
        
                 :type: LinkEndType
        """

    @property
    def residuals(self) -> list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
        No documentation found.
        """

    @property
    def rms_residuals(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
        No documentation found.
        """

    @property
    def single_observable_size(self) -> int:
        """
        No documentation found.
        """

    @property
    def time_bounds(self) -> tuple[time_representation.Time, time_representation.Time]:
        """
        No documentation found.
        """

    @property
    def total_observation_set_size(self) -> int:
        """
        No documentation found.
        """

    @property
    def weights(self) -> list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
        No documentation found.
        """

    @property
    def weights_vector(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
        No documentation found.
        """

    @weights_vector.setter
    def weights_vector(self, arg1: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]) -> None:
        ...

def compute_residuals_and_dependent_variables(observation_collection: ObservationCollection, observation_simulators: list[observables_simulation.ObservationSimulator], bodies: environment.SystemOfBodies) -> None:
    """No documentation found."""

@typing.overload
def create_filtered_observation_collection(original_observation_collection: ObservationCollection, observation_filters_map: dict[observations_processing.ObservationCollectionParser, observations_processing.ObservationFilterBase]) -> ObservationCollection:
    """No documentation found."""

@typing.overload
def create_filtered_observation_collection(original_observation_collection: ObservationCollection, observation_filter: observations_processing.ObservationFilterBase, observation_parser: observations_processing.ObservationCollectionParser=...) -> ObservationCollection:
    """No documentation found."""

def create_new_observation_collection(original_observation_collection: ObservationCollection, observation_parser: observations_processing.ObservationCollectionParser=...) -> ObservationCollection:
    """No documentation found."""

def create_single_observation_set(observable_type: model_settings.ObservableType, link_ends: dict[links.LinkEndType, links.LinkEndId], observations: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], observation_times: list[time_representation.Time], reference_link_end: links.LinkEndType, ancillary_settings: ancillary_settings.ObservationAncilliarySimulationSettings) -> SingleObservationSet:
    """No documentation found."""

def filter_observations(original_observation_set: SingleObservationSet, observation_filter: observations_processing.ObservationFilterBase, save_filtered_observations: bool=False) -> SingleObservationSet:
    """No documentation found."""

def merge_observation_collections(observation_collection_list: list[ObservationCollection]) -> ObservationCollection:
    ...

def single_observation_set(observable_type: model_settings.ObservableType, link_definition: links.LinkDefinition, observations: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], observation_times: list[time_representation.Time], reference_link_end: links.LinkEndType, ancilliary_settings: ancillary_settings.ObservationAncilliarySimulationSettings=None) -> SingleObservationSet:
    """No documentation found."""

def split_observation_collection(original_observation_collection: ObservationCollection, observation_set_splitter: observations_processing.ObservationSetSplitterBase, observation_parser: observations_processing.ObservationCollectionParser=...) -> ObservationCollection:
    """No documentation found."""

def split_observation_set(original_observation_set: SingleObservationSet, observation_splitter: observations_processing.ObservationSetSplitterBase, print_warning: bool=True) -> list[SingleObservationSet]:
    """No documentation found."""