from ....estimation.observable_models_setup import links
from ....estimation.observable_models_setup import model_settings
from ....estimation.observations_setup import observations_simulation_settings
import typing
__all__ = ['ObservationViabilitySettings', 'ObservationViabilityType', 'add_viability_check_to_all', 'add_viability_check_to_observable', 'add_viability_check_to_observable_for_link_ends', 'body_avoidance_angle', 'body_avoidance_viability', 'body_avoidance_viability_list', 'body_occultation', 'body_occultation_viability', 'body_occultation_viability_list', 'elevation_angle_viability', 'elevation_angle_viability_list', 'minimum_elevation_angle']

class ObservationViabilitySettings:
    """Class for defining observation viability calculator settings.
    
    Class for defining the settings for observation viability calculator creation.
    Instances of this class are typically be created through various dedicated functions,such as :func:`~tudatpy.estimation.observations_setup.viability.elevation_angle_viability`, :func:`~tudatpy.estimation.observations_setup.viability.body_avoidance_viability` and :func:`~tudatpy.estimation.observations_setup.viability.body_occultation_viability`
    
    Examples
    --------
    .. code-block:: python
    
        # Code snippet to show the creation of an ObservationViabilitySettings object
        import numpy as np
        from tudatpy.estimation.observations_setup import viability
    
        # Create ObservationViabilitySettings object
        # In this case, we exclude observations for which the local elevation angle at link end is less 15 degrees.
        min_elevation = np.deg2rad(15)
        # We apply these settings to every ground station on Earth using the following link_end_id: [“Earth”, “”]
        viability_settings = viability.elevation_angle_viability(["Earth", ""], min_elevation)
    
        # Show that this is indeed an ObservationViabilitySettings object
        print(viability_settings)"""

class ObservationViabilityType:
    """Enumeration of observation viability criterion types.
    
    Examples
    --------
    .. code-block:: python
    
        # Code snippet to print all available Observation Viability Types
        from tudatpy.estimation.observations_setup import viability
    
        num_observation_viability_types = len(viability.ObservationViabilityType.__members__)
        print(f'The length of all available Tudatpy Observation Viability Types is: {num_observation_viability_types}')
    
        # Print all available Observation Viability Types using the "name" property
        for i in range(num_observation_viability_types):
            print(i, viability.ObservationViabilityType(i).name)
    
    
    
    
          
    
    Members:
    
      minimum_elevation_angle
    
      body_avoidance_angle
    
      body_occultation"""
    __members__: typing.ClassVar[dict[str, ObservationViabilityType]]
    body_avoidance_angle: typing.ClassVar[ObservationViabilityType]
    body_occultation: typing.ClassVar[ObservationViabilityType]
    minimum_elevation_angle: typing.ClassVar[ObservationViabilityType]

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...

def add_viability_check_to_all(observation_simulation_settings_list: list[observations_simulation_settings.ObservationSimulationSettings], viability_settings: list[ObservationViabilitySettings]) -> None:
    """Function for including viability checks into existing observation simulation settings.
    
    Function for adding viability checks to the observation simulation settings, such that only observations meeting certain conditions are retained.
    The noise settings are added to all :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` object(s) in the `observation_simulation_settings`
    list.
    Note: the :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` objects are modified in-place by this function,
    and thus the function does not return anything.
    
    
    Parameters
    ----------
    observation_simulation_settings : List[ :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` ]
        Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` objects.
    viability_settings : List[ :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` ]
        List of one or more :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, defining the viability checks to be included.
    
    Returns
    -------
    None
        The :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` object(s) are changed in-place."""

def add_viability_check_to_observable(observation_simulation_settings_list: list[observations_simulation_settings.ObservationSimulationSettings], viability_settings: list[ObservationViabilitySettings], observable_type: model_settings.ObservableType) -> None:
    """Function for including viability checks into existing observation simulation settings.
    
    As :func:`~tudatpy.estimation.observations_setup.viability.add_viability_check_to_all`, except that the function only adds viabilitt settings to entries of the
    `observation_simulation_settings` list that matches the specified `observable_type`.
    
    
    Parameters
    ----------
    observation_simulation_settings : List[ :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` ]
        Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` objects.
    viability_settings : List[ :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` ]
        List of one or more :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, defining the viability checks to be included.
    
    observable_type : :class:`ObservableType`
        Identifies the observable type in the observation simulation settings for which the viability checks are to be considered.
    
    Returns
    -------
    None
        The :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` object(s) are changed in-place."""

def add_viability_check_to_observable_for_link_ends(observation_simulation_settings_list: list[observations_simulation_settings.ObservationSimulationSettings], viability_settings: list[ObservationViabilitySettings], observable_type: model_settings.ObservableType, link_ends: links.LinkDefinition) -> None:
    """Function for including viability checks into existing observation simulation settings.
    
    As :func:`~tudatpy.estimation.observations_setup.viability.add_viability_check_to_all`, except that the function only adds noise to entries of the
    `observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`.
    
    
    Parameters
    ----------
    observation_simulation_settings : List[ :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` ]
        Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` objects.
    viability_settings : List[ :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` ]
        List of one or more :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, defining the viability checks to be included.
    
    observable_type : :class:`ObservableType`
        Identifies the observable type in the observation simulation settings for which the viability checks are to be considered.
    
    link_definition : :class:`~tudatpy.estimation.observable_models_setup.links.LinkDefinition`
        Identifies the link definition in the observation simulation settings for which the viability checks are to be considered.
    
    Returns
    -------
    None"""

def body_avoidance_viability(link_end_id: tuple[str, str], body_to_avoid: str, avoidance_angle: float) -> ObservationViabilitySettings:
    """Function for defining body avoidance observation viability settings.
    
    Function for defining body avoidance observation viability settings for single link ends.
    When simulating observations, this settings ensures that any applicable observations, for which the signal path passes 'too close' to a body, will be omitted.
    The definition of 'too close' is computed as the angle between:
    
    * The line-of-sight vector from a link end to a given third body
    * The line-of-sight between two link ends
    
    This constraint is typically used to prevent the Sun from being too close to the field-of-view of the telescope(s), as defined by
    a so-called 'SPE' (Sun-Probe-Earth) angle constraint. The present viability setting generalizes this constraint.
    
    
    Parameters
    ----------
    link_end_id : Tuple[str,str]
        Link end (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId` ), for which the viability settings are to be created.
        To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
        For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
        link end passes too close to the specified body.
    
    body_to_avoid : str
        Name of the body which the signal path should not pass 'too close' to.
    
    avoidance_angle : float
        Limit angle (generalization of SPE angle), below which no observations are produced when using the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. Note: this
        value must be in radians.
    
    Returns
    -------
    :class:`ObservationViabilitySettings`
        Instance of the :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings`, defining the settings for observation viability."""

def body_avoidance_viability_list(link_end_ids: list[tuple[str, str]], body_to_avoid: str, avoidance_angle: float) -> list[ObservationViabilitySettings]:
    """Function for defining list of body avoidance viability settings.
    
    Function for defining body avoidance viability settings for multiple link ends.
    Each entry in the returned list contains the observation viability settings for one link end.
    When simulating observations, these settings ensure that any applicable observations, for which the signal path passes 'too close' to a body, will be omitted.
    The definition of 'too close' is computed as the angle between:
    
    * The line-of-sight vector from a link end to a given third body
    * The line-of-sight between two link ends
    
    This constraint is typically used to prevent the Sun from being too close to the field-of-view of the telescope(s), as defined by
    a so-called 'SPE' (Sun-Probe-Earth) angle constraint. The present viability setting generalizes this constraint.
    
    
    Parameters
    ----------
    link_end_ids : List[ Tuple[str,str] ]
        List of individual link ends (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the elevation angle viability setting is to be created.
        To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].
    
    body_to_avoid : str
        Name of the body which the signal path should not pass 'too close' to.
    
    avoidance_angle : float
        Limit angle (generalization of SPE angle), below which no observations are produced when using the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. Note: this
        value must be in radians.
    
    Returns
    -------
    :class:`ObservationViabilitySettings`
        List of :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end."""

def body_occultation_viability(link_end_id: tuple[str, str], occulting_body: str) -> ObservationViabilitySettings:
    """Function for defining body occultation viability settings.
    
    Function for defining body occultation viability settings for single link ends.
    When simulating observations, this setting ensures that any applicable observations, for which the signal path is occulted by a given body, will be omitted.
    The occultation is computed using the shape model of the specified body.
    
    
    Parameters
    ----------
    link_end_id : Tuple[str,str]
        Link end (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the viability settings are to be created.
        To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
    
    body_to_avoid : str
        Name of the body which the signal path should not be occulted by.
    
    Returns
    -------
    :class:`ObservationViabilitySettings`
        Instance of the :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings`, defining the settings for observation viability."""

def body_occultation_viability_list(link_end_ids: list[tuple[str, str]], occulting_body: str) -> list[ObservationViabilitySettings]:
    """Function for defining body occultation viability settings.
    
    Function for defining body occultation viability settings for multiple link ends.
    Each entry in the returned list contains the observation viability settings for one link end.
    When simulating observations, these settings ensure that any applicable observations, for which the signal path is occulted by a given body, will be omitted.
    The occultation is computed using the shape model of the specified body.
    
    
    Parameters
    ----------
    link_end_ids : List[ Tuple[str,str] ]
        List of individual link ends (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the viability settings are to be created.
        To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
        For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
        link end is occulted by the specified body.
    
    body_to_avoid : str
        Name of the body which the signal path should not be occulted by.
    
    Returns
    -------
    :class:`ObservationViabilitySettings`
        List of :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end."""

def elevation_angle_viability(link_end_id: tuple[str, str], elevation_angle: float) -> ObservationViabilitySettings:
    """Function for defining single elevation angle viability setting.
    
    Function for defining elevation angle viability settings for single link end.
    When simulating observations, this setting ensures that any applicable observations, for which the local elevation angle at link end is less than some limit value, will be omitted.
    
    
    Parameters
    ----------
    link_end_id : Tuple[str,str]
        Link end (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the elevation angle viability setting is to be created.
        To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].
    
    elevation_angle : float
        Limit elevation angle, below which no observations are produced when using the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. Note: this
        value must be in radians.
    
    Returns
    -------
    :class:`ObservationViabilitySettings`
        Instance of the :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` class, defining the settings for observation viability"""

def elevation_angle_viability_list(link_end_ids: list[tuple[str, str]], elevation_angle: float) -> list[ObservationViabilitySettings]:
    """Function for defining list of elevation angle viability settings.
    
    Function for defining elevation angle viability settings for multiple link ends.
    Each entry in the returned list contains the observation viability settings for one link end.
    When simulating observations, these settings ensure that any applicable observations, for which the local elevation angle at a link end is less than some limit value, will be omitted.
    
    
    Parameters
    ----------
    link_end_ids : List[ Tuple[str,str] ]
        List of individual link ends (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the elevation angle viability setting is to be created.
        To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].
        For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
        link end violates the minimum elevation angle constraint.
    
    elevation_angle : float
        Limit elevation angle, below which no observations are produced when using the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. Note: this
        value must be in radians.
    
    Returns
    -------
    :class:`ObservationViabilitySettings`
        List of :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end."""
body_avoidance_angle: ObservationViabilityType
body_occultation: ObservationViabilityType
minimum_elevation_angle: ObservationViabilityType