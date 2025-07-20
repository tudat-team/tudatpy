import numpy
import pybind11_stubgen.typing_ext
from ....estimation.observable_models_setup import links
from ....estimation.observable_models_setup import model_settings
from ....estimation.observations_setup import observations_simulation_settings
import typing
__all__ = ['add_gaussian_noise_to_all', 'add_gaussian_noise_to_observable', 'add_gaussian_noise_to_observable_for_link_ends', 'add_noise_function_to_all', 'add_noise_function_to_observable', 'add_noise_function_to_observable_for_link_ends']

def add_gaussian_noise_to_all(observation_simulation_settings_list: list[observations_simulation_settings.ObservationSimulationSettings], noise_amplitude: float) -> None:
    """Function for adding gaussian noise function to all existing observation simulation settings.
    
    Function for including simple time-independent and time-uncorrelated Gaussian noise function to the simulation settings of one or more observable(s).
    The noise settings are added to all :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` object(s) in the `observation_simulation_settings`
    list.
    
    Note: the :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects are modified in-place by this function,
    and thus the function does not return anything.
    
    
    Parameters
    ----------
    observation_simulation_settings : List[ :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` ]
        Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects.
    noise_amplitude : float
        Standard deviation defining the un-biased Gaussian distribution for the noise.
    Returns
    -------
    None
        The :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` object(s) are changed in-place."""

def add_gaussian_noise_to_observable(observation_simulation_settings_list: list[observations_simulation_settings.ObservationSimulationSettings], noise_amplitude: float, observable_type: model_settings.ObservableType) -> None:
    """Function for adding gaussian noise function to existing observation simulation settings of a given observable type.
    
    As :func:`~tudatpy.estimation.observations_setup.random_noise.add_gaussian_noise_to_all`, except that the function only adds noise to entries of the
    `observation_simulation_settings` list that matches the specified `observable_type`.
    
    
    Parameters
    ----------
    observation_simulation_settings : List[ :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` ]
        Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects.
    noise_amplitude : float
        Standard deviation defining the un-biased Gaussian distribution for the noise.
    observable_type : :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
        Identifies the observable type in the observation simulation settings to which the noise is to be added.
    
    Returns
    -------
    None
        The :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` object(s) are changed in-place."""

def add_gaussian_noise_to_observable_for_link_ends(observation_simulation_settings_list: list[observations_simulation_settings.ObservationSimulationSettings], noise_amplitude: float, observable_type: model_settings.ObservableType, link_definition: links.LinkDefinition) -> None:
    """Function for adding gaussian noise function to existing observation simulation settings of a given observable type and link definition.
    
    As :func:`~tudatpy.estimation.observations_setup.random_noise.add_gaussian_noise_to_all`, except that the function only adds noise to entries of the
    `observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`.
    
    
    Parameters
    ----------
    observation_simulation_settings : List[ :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` ]
        Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects.
    noise_amplitude : float
        Standard deviation defining the un-biased Gaussian distribution for the noise.
    observable_type : :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
        Identifies the observable type in the observation simulation settings to which the noise is to be added.
    
    link_definition : :class:`~tudatpy.estimation.observable_models_setup.links.LinkDefinition`
        Identifies the link definition in the observation simulation settings for which the noise is to be added.
    
    Returns
    -------
    None
        The :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` object(s) are changed in-place."""

def add_noise_function_to_all(observation_simulation_settings_list: list[observations_simulation_settings.ObservationSimulationSettings], noise_amplitude: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]) -> None:
    """Function for adding a custom noise function to all existing observation simulation settings.
    
    Function for including a custom noise function to the simulation settings of all observables.
    The noise settings are added to all :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` object(s) in the `observation_simulation_settings`
    list.
    
    Note: the :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects are modified in-place by this function,
    and thus the function does not return anything.
    
    
    Parameters
    ----------
    observation_simulation_settings_list : List[ :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` ]
        Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects.
    
    noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ]
        Function providing the observation noise factors as a function of observation time.
    
    Returns
    -------
    None
        The :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` object(s) are changed in-place."""

def add_noise_function_to_observable(observation_simulation_settings_list: list[observations_simulation_settings.ObservationSimulationSettings], noise_amplitude: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], observable_type: model_settings.ObservableType) -> None:
    """Function for adding a custom noise function to selected existing observation simulation settings of a given observable type.
    
    As :func:`~tudatpy.estimation.observations_setup.random_noise.add_noise_function_to_all`, except that the function only adds noise to entries of the
    `observation_simulation_settings` list that matches the specified `observable_type`.
    
    
    Parameters
    ----------
    observation_simulation_settings_list : List[ :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` ]
        Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects.
    
    noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ]
        Function providing the observation noise factors as a function of observation time.
    
    observable_type : :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
        Identifies the observable type in the observation simulation settings to which the noise is to be added.
    
    Returns
    -------
    None
        The :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` object(s) are changed in-place."""

def add_noise_function_to_observable_for_link_ends(observation_simulation_settings_list: list[observations_simulation_settings.ObservationSimulationSettings], noise_amplitude: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], observable_type: model_settings.ObservableType, link_ends: links.LinkDefinition) -> None:
    """Function for adding a custom noise function to existing observation simulation settings of a given observable type and link definition.
    
    As :func:`~tudatpy.estimation.observations_setup.random_noise.add_noise_function_to_all`, except that the function only adds noise to entries of the
    `observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`.
    
    
    Parameters
    ----------
    observation_simulation_settings : List[ :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` ]
        Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects.
    
    noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ]
        Function providing the observation noise factors as a function of observation time.
    
    observable_type : :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
        Identifies the observable type in the observation simulation settings to which the noise is to be added.
    
    link_definition : :class:`~tudatpy.estimation.observable_models_setup.links.LinkDefinition`
        Identifies the link definition in the observation simulation settings for which the noise is to be added.
    
    Returns
    -------
    None
        The :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` object(s) are changed in-place."""