import numpy
import pybind11_stubgen.typing_ext
from ....astro import time_representation
from ....dynamics import environment
from ....estimation.observable_models import observables_simulation
from ....estimation.observable_models_setup import links
from ....estimation.observable_models_setup import model_settings
from ....estimation.observations_setup import ancillary_settings
import typing
__all__ = ['ObservationSimulationSettings', 'TabulatedObservationSimulationSettings', 'change_simulation_settings_observable_types', 'continuous_arc_simulation_settings', 'continuous_arc_simulation_settings_list', 'create_observation_simulators', 'observation_settings_from_collection', 'tabulated_simulation_settings', 'tabulated_simulation_settings_list']

class ObservationSimulationSettings:
    """Base class for defining settings for simulated observations."""

    @property
    def noise_function(self) -> typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ], default = None -
                 Function providing the observation noise as a function of observation time (can be constant or time-dependent), default is None.
        """

    @noise_function.setter
    def noise_function(self, arg1: typing.Callable[[float], float]) -> None:
        ...

    @property
    def viability_settings_list(self) -> list[...]:
        """
                 viability_settings_list : List[ ObservationViabilitySettings ], default = [ ]) -
                 Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
        """

    @viability_settings_list.setter
    def viability_settings_list(self, arg1: list[...]) -> None:
        ...

class TabulatedObservationSimulationSettings(ObservationSimulationSettings):
    """Class for defining settings for simulating observations at a predefined set of times.
    
    Class for defining settings for simulating observations at a predefined set of times.
    This type defines predefined time epochs at which applicable observations are to be simulated, stored in a rigid, "tabulated" form.
    Some observation times may be discarded due to the use of viability settings.
    Instances of this class are typically created via the :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.tabulated_simulation_settings`
    and :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.tabulated_simulation_settings_list` functions.
    
    Associated base class: :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings`
    
    Examples
    --------
    .. code-block:: python
    
        # Code snippet to show the creation of a TabulatedObservationSimulationSettings object
        import numpy as np
        from tudatpy.astro.time_representation import DateTime
        from tudatpy.estimation.observable_models_setup import links, model_settings
        from tudatpy.estimation.observations_setup import observations_simulation_settings
    
        # Set simulation start and end epochs
        simulation_start_epoch = DateTime(2000, 1, 1).epoch()
        simulation_end_epoch   = DateTime(2000, 1, 4).epoch()
    
        # Define the uplink link ends for one-way observable
        link_ends = dict()
        link_ends[links.transmitter] = links.body_origin_link_end_id("Earth")
        link_ends[links.receiver] = links.body_origin_link_end_id("Delfi-C3")
    
        # Create LinkDefinition Object and set observation settings for each link/observable
        link_definition = links.LinkDefinition(link_ends)
        observation_settings_list = [model_settings.one_way_doppler_instantaneous(link_definition)]
    
        # Define observation simulation times (separated by steps of 1 minute)
        observation_times = np.arange(simulation_start_epoch, simulation_end_epoch, 60.0)
    
        # Create TabulatedObservationSimulationSettings object
        tabulated_observation_simulation_settings = observations_simulation_settings.tabulated_simulation_settings(
            model_settings.one_way_instantaneous_doppler_type,
            link_definition,
            observation_times
        )
    
        # Show that this is indeed a TabulatedObservationSimulationSettings object
        print(tabulated_observation_simulation_settings)"""

def change_simulation_settings_observable_types(observation_simulation_settings: list[ObservationSimulationSettings], replacement_observable_types: dict[model_settings.ObservableType, model_settings.ObservableType]=...) -> None:
    """No documentation found."""

def continuous_arc_simulation_settings(observable_type: model_settings.ObservableType, link_ends: links.LinkDefinition, start_time: time_representation.Time, end_time: time_representation.Time, interval_between_observations: time_representation.Time, arc_limiting_constraints: ..., minimum_arc_duration: time_representation.Time, maximum_arc_duration: time_representation.Time, minimum_time_between_arcs: time_representation.Time, reference_link_end_type: links.LinkEndType=..., additional_viability_settings: list[...]=[], noise_function: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]=None) -> ObservationSimulationSettings:
    """Function for creating settings object for observation simulation, using observation times according to a requirement for a continuous tracking arc.
    
    Function for creating settings object for observation simulation. Unlike the :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.tabulated_simulation_settings`
    function, the resulting settings do not define the observation times explicitly. Instead, this settings object determines the observation times adaptively during the
    simulation of the observation, with the requirement that observations should be simulated over a set of contiguous arcs (if possible). The exact algorithm meets the following conditions:
    
    * Observations are only simulated within the time span of ``start_time`` and ``end_time``
    * A contiguous tracking arc has simulated observations separated by ``interval_between_observations``
    * Starting from ``start_time``, an observation is simulated each ``interval_between_observations``. Once an observation is unviable, as defined by
      the ``arc_limiting_constraints`` input, it is checked whether the arc up until that point
      is longer in duration than ``minimum_arc_duration``. If it is, the arc is added to the simulated observations. If not, the arc is discarded. In either case, a new arc is started once a
      viable is observation is encountered
    * If the current arc reaches a duration greater than ``maximum_arc_duration``, the arc is added to the existing observations, and a new arc is started
    * If defined (e.g. if not NaN), the current observation time is incremented by ``minimum_time_between_arcs`` when an arc has been added to the observations.
    
    Nominally, this algorithm ensures that any arc of observations has a minimum and maximum duration. In addition, it ensures that (if desired) there is a minimum time interval
    between two tracking arcs. This behaviour can be modified by adding ``additional_viability_settings``, which are *not* used when computing the tracking arcs, but which are instead only used
    to reduce the set of simulated observations afterwards.
    
    
    Parameters
    ----------
    observable_type : :class:`ObservableType`
        Observable type of which observations are to be simulated.
    link_ends : LinkDefinition
        Link ends for which observations are to be simulated.
    start_time : float
        First time at which an observation is to be simulated (and checked for viability).
    end_time : float
        Maximum time at which an observation is to be simulated (and checked for viability).
    interval_between_observations : float
        Cadence (in seconds) of subsequent observations in an arc
    arc_limiting_constraints : List[ :class:`ObservationViabilitySettings` ], default = [ ]
        List of settings for the creation of the viability criteria calculators, which are used to check if an observation is viable, and define
        whether an arc should be terminated.
    
    minimum_arc_duration : float
        Minimum permissible time for a tracking arc
    maximum_arc_duration : float
        Maximum permissible time for a tracking arc
    minimum_time_between_arc : float, default = NaN
        Minimum time between two tracking arcs. If NaN, this is effectively set to the ``interval_between_observations``
    additional_viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
        Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
        These settings are *not* used to determine whether an arc is to be terminated, but are instead applied after the arcs have been computed.
    
    noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ], default = None
        Function providing the observation noise factors as a function of observation time.
    Returns
    -------
    :class:`TabulatedObservationSimulationSettings`
        Instance of the :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` derived :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.TabulatedObservationSimulationSettings` class."""

def continuous_arc_simulation_settings_list(link_ends_per_observable: dict[model_settings.ObservableType, list[links.LinkDefinition]], start_time: time_representation.Time, end_time: time_representation.Time, interval_between_observations: time_representation.Time, arc_limiting_constraints: ..., minimum_arc_duration: time_representation.Time, maximum_arc_duration: time_representation.Time, minimum_time_between_arcs: time_representation.Time, reference_link_end_type: links.LinkEndType=..., additional_viability_settings: list[...]=[]) -> list[ObservationSimulationSettings]:
    """Function for creating a list of settings object for observation simulation, using observation times according to a requirement for a continuous tracking arc.
    
    Function for creating multiple settings objects for observation simulation in a list. This function is
    equivalent to calling the :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.continuous_arc_simulation_settings` repeatedly, with the different
    observables and link definition provided here through `link_ends_per_observable`.
    During a single call to this function, one simulation settings object is created for each combination of observable type and link geometry given by the `link_ends_per_observable` parameter.
    
    
    Parameters
    ----------
    link_ends_per_observable : Dict[:class:`ObservableType`, List[LinkDefinition]]]
        Link geometry per observable type of which observations are to be simulated.
    start_time : float
        First time at which an observation is to be simulated (and checked for viability).
    end_time : float
        Maximum time at which an observation is to be simulated (and checked for viability).
    interval_between_observations : float
        Cadence (in seconds) of subsequent observations in an arc
    arc_limiting_constraints : List[ :class:`ObservationViabilitySettings` ], default = [ ]
        List of settings for the creation of the viability criteria calculators, which are used to check if an observation is viable, and define
        whether an arc should be terminated.
    
    minimum_arc_duration : float
        Minimum permissible time for a tracking arc
    maximum_arc_duration : float
        Maximum permissible time for a tracking arc
    minimum_time_between_arc : float, default = NaN
        Minimum time between two tracking arcs. If NaN, this is effectively set to the ``interval_between_observations``
    additional_viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
        Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
        These settings are *not* used to determine whether an arc is to be terminated, but are instead applied after the arcs have been computed.
    
    Returns
    -------
    List[ :class:`TabulatedObservationSimulationSettings` ]
        List of :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` derived :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.TabulatedObservationSimulationSettings` objects."""

def create_observation_simulators(observation_settings: list[model_settings.ObservationSettings], bodies: environment.SystemOfBodies) -> list[observables_simulation.ObservationSimulator]:
    """Function for creating observation simulator objects.
    
    Function for creating observation simulator objects from observation settings.
    Note that each observation (i.e. combination of observable and link geometry) requires its own observation simulator object.
    
    
    Parameters
    ----------
    observation_settings : List[ ObservationSettings ]
        List of settings objects, each object defining the observation model settings for one combination of observable and link geometry that is to be simulated.
    
    bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
        Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.
    
    Returns
    -------
    List[ :class:`~tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator` ]
        List of :class:`~tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator` objects, each object hosting the functionality for simulating one combination of observable type and link geometry.
    
    Examples
    --------
    .. code-block:: python
    
        from tudatpy.estimation.observations_setup import observations_simulation_settings
    
        # Create bodies
        bodies = ...
        # Define parameters settings
        observation_settings = ...
        # Create observation simulators
        observation_simulators = observations_simulation_settings.create_observation_simulators(observation_settings, bodies)
    
    This code snippet closely follows what is done in: The following snippet closely follows what is done in: `Galilean Moons State Estimation Example <https://github.com/tudat-team/tudatpy-examples/blob/master/estimation/galilean_moons_state_estimation.ipynb>`_."""

def observation_settings_from_collection(*args, **kwargs) -> list[ObservationSimulationSettings]:
    """No documentation found."""

def tabulated_simulation_settings(observable_type: model_settings.ObservableType, link_ends: links.LinkDefinition, simulation_times: list[time_representation.Time], reference_link_end_type: links.LinkEndType=..., viability_settings: list[...]=[], noise_function: typing.Callable[[float], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]=None, ancilliary_settings: ancillary_settings.ObservationAncilliarySimulationSettings=None) -> ObservationSimulationSettings:
    """Function for creating settings object for observation simulation, using a predefined list of observation times.
    
    Function for creating single simulation settings object, using a predefined list of observation times.
    The list of resulting observations may be reduced compared to the ``simulation_times`` should be provided here, as
    only observations that meet the viability settings are retained during observation simulation (these may be
    provide directly here through the ``viability_settings`` input, or added later to the resulting settings object).
    
    
    Parameters
    ----------
    observable_type : :class:`ObservableType`
        Observable type of which observations are to be simulated.
    link_ends : LinkDefinition
        Link ends for which observations are to be simulated.
    simulation_times : List[float]
        List of times at which to perform the observation simulation.
    reference_link_end_type : :class:`LinkEndType`, default = :class:`LinkEndType.receiver`
        Defines the link end (via the :class:`LinkEndType`) which is used as a reference time for the observation.
    viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
        Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
    
    noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ], default = None
        Function providing the observation noise factors as a function of observation time.
    Returns
    -------
    :class:`TabulatedObservationSimulationSettings`
        Instance of the :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` derived :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.TabulatedObservationSimulationSettings` class."""

def tabulated_simulation_settings_list(link_ends_per_observable: dict[model_settings.ObservableType, list[links.LinkDefinition]], simulation_times: list[time_representation.Time], reference_link_end_type: links.LinkEndType=..., viability_settings: list[...]=[]) -> list[ObservationSimulationSettings]:
    """Function for creating a list of settings object for observation simulation, using a predefined list of observation times.
    
    Function for creating multiple tabulated observation simulation settings objects in a list. This function is
    equivalent to calling the :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.tabulated_simulation_settings` repeatedly, with the different
    observables and link definition provided here through `link_ends_per_observable`.
    During a single call to this function, one simulation settings object is created for each combination of observable type and link geometry given by the `link_ends_per_observable` parameter.
    
    
    Parameters
    ----------
    link_ends_per_observable : Dict[:class:`ObservableType`, List[LinkDefinition]]]
        Link geometry per observable type of which observations are to be simulated.
    simulation_times : List[ float ]
        List of times at which to perform the observation simulation.
    reference_link_end_type : :class:`LinkEndType`, default = :class:`LinkEndType.receiver`
        Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.
        The single link end specified here will be considered as the reference link end for all simulation settings object created in the function call.
    
    viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
        Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
        The single settings list given here will be considered as potential viability settings for all simulation settings object created in the function call.
    
    Returns
    -------
    List[ TabulatedObservationSimulationSettings ]
        List of :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` derived :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.TabulatedObservationSimulationSettings` objects."""