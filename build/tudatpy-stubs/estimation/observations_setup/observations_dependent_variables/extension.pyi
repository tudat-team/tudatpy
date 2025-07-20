from ....dynamics import environment
from ....estimation.observable_models_setup import links
from ....estimation.observable_models_setup import model_settings
import typing
__all__ = ['IntegratedObservationPropertyHandling', 'ObservationDependentVariableSettings', 'add_dependent_variables_to_all', 'add_dependent_variables_to_observable', 'add_dependent_variables_to_observable_for_link_ends', 'angle_wrt_orbital_plane_dependent_variable', 'avoidance_angle_dependent_variable', 'azimuth_angle_dependent_variable', 'body_center_distance_dependent_variable', 'body_limb_distance_dependent_variable', 'elevation_angle_dependent_variable', 'integration_time_dependent_variable', 'interval_end', 'interval_start', 'interval_undefined', 'retransmission_delays_dependent_variable', 'target_range_between_link_ends_dependent_variable']

class IntegratedObservationPropertyHandling:
    """No documentation found.
    
    Members:
    
      interval_start
    
      interval_end
    
      interval_undefined"""
    __members__: typing.ClassVar[dict[str, IntegratedObservationPropertyHandling]]
    interval_end: typing.ClassVar[IntegratedObservationPropertyHandling]
    interval_start: typing.ClassVar[IntegratedObservationPropertyHandling]
    interval_undefined: typing.ClassVar[IntegratedObservationPropertyHandling]

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

class ObservationDependentVariableSettings:
    """Base class for setting observation dependent variables as part of the observation output.
    
    Base class for setting observation dependent variables as part of the observation output.
    The user can create instances of this class via the :func:`~tudatpy.estimation.observations_setup.observations_dependent_variables.elevation_angle_dependent_variable` function.
    Note: The associated functionality is not yet mature enough for the end user. Class is exposed for development purposes only.
    
    Examples
    --------
    .. code-block:: python
    
        # Code snippet to show the creation of an ObservationDependentVariableSettings object
        from tudatpy.estimation.observations_setup import observations_dependent_variables
        from tudatpy.estimation.observable_models_setup import links
    
        # Create ObservationDependentVariableSettings object
        elevation_angle_settings = observations_dependent_variables.elevation_angle_dependent_variable(links.receiver)
    
        # Show that this is indeed an ObservationDependentVariableSettings object
        print(elevation_angle_settings)"""

def add_dependent_variables_to_all(observation_simulation_settings: list[...], dependent_variable_settings: list[ObservationDependentVariableSettings], bodies: environment.SystemOfBodies) -> None:
    """Function for including dependent variables into all existing observation simulation settings.
    
    Function for including the computation and reporting of dependent variables into the observation simulation settings of all observables.
    Note: The associated functionality is not yet mature enough for the end user. Function is exposed for development purposes only.
    
    Modifications are applied to all given :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` object(s),
    matching each :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` object with the corresponding :class:`ObservationDependentVariableSettings` entry in the `dependent_variable_settings` parameter.
    Note that the :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects are modified in-place and thus the function does not return anything.
    
    
    Parameters
    ----------
    observation_simulation_settings : List[ :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` ]
        Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects.
    
    dependent_variable_settings : List[ :class:`~tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings` ]
        List of one or more :class:`~tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings` objects, defining the dependent variables to be considered.
    
    bodies : :class:`~tudatpy.dynamics.environment_setup.SystemOfBodies`
        Object consolidating all bodies and environment models that constitute the physical environment."""

def add_dependent_variables_to_observable(observation_simulation_settings: list[...], dependent_variable_settings: list[ObservationDependentVariableSettings], bodies: environment.SystemOfBodies, observable_type: model_settings.ObservableType) -> None:
    """Function for including dependent variables into selected existing observation simulation settings.
    
    As :func:`~tudatpy.estimation.observations_setup.observations_dependent_variables.add_dependent_variables_to_all`, except that the function only adds includes the
    computation and reporting of dependent variables to entries of the `observation_simulation_settings` list that matches the specified `observable_type`.
    
    
    Parameters
    ----------
    observation_simulation_settings : List[ :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` ]
        Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects.
    
    dependent_variable_settings : List[ :class:`~tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings` ]
        List of one or more :class:`~tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings` objects, defining the dependent variables to be considered.
    
    bodies : :class:`~tudatpy.dynamics.environment_setup.SystemOfBodies`
        Object consolidating all bodies and environment models that constitute the physical environment.
    
    observable_type : :class:`ObservableType`
        Identifies the observable type in the observation simulation settings for which the dependent variables are to be included."""

def add_dependent_variables_to_observable_for_link_ends(observation_simulation_settings: list[...], dependent_variable_settings: list[ObservationDependentVariableSettings], bodies: environment.SystemOfBodies, observable_type: model_settings.ObservableType, link_ends: links.LinkDefinition) -> None:
    """No documentation found."""

def angle_wrt_orbital_plane_dependent_variable(body_name: str, start_link_end_type: links.LinkEndType=..., end_link_end_type: links.LinkEndType=..., start_link_end_id: links.LinkEndId=..., end_link_end_id: links.LinkEndId=..., integrated_observation_handling: IntegratedObservationPropertyHandling=...) -> ObservationDependentVariableSettings:
    """No documentation found."""

def avoidance_angle_dependent_variable(body_name: str, start_link_end_type: links.LinkEndType=..., end_link_end_type: links.LinkEndType=..., start_link_end_id: links.LinkEndId=..., end_link_end_id: links.LinkEndId=..., integrated_observation_handling: IntegratedObservationPropertyHandling=...) -> ObservationDependentVariableSettings:
    """No documentation found."""

def azimuth_angle_dependent_variable(link_end_type: links.LinkEndType=..., link_end_id: links.LinkEndId=..., originating_link_end_type: links.LinkEndType=..., originating_link_end_id: links.LinkEndId=..., integrated_observation_handling: IntegratedObservationPropertyHandling=...) -> ObservationDependentVariableSettings:
    """No documentation found."""

def body_center_distance_dependent_variable(body_name: str, start_link_end_type: links.LinkEndType=..., end_link_end_type: links.LinkEndType=..., start_link_end_id: links.LinkEndId=..., end_link_end_id: links.LinkEndId=..., integrated_observation_handling: IntegratedObservationPropertyHandling=...) -> ObservationDependentVariableSettings:
    """No documentation found."""

def body_limb_distance_dependent_variable(body_name: str, start_link_end_type: links.LinkEndType=..., end_link_end_type: links.LinkEndType=..., start_link_end_id: links.LinkEndId=..., end_link_end_id: links.LinkEndId=..., integrated_observation_handling: IntegratedObservationPropertyHandling=...) -> ObservationDependentVariableSettings:
    """No documentation found."""

def elevation_angle_dependent_variable(link_end_type: links.LinkEndType=..., link_end_id: links.LinkEndId=..., originating_link_end_type: links.LinkEndType=..., originating_link_end_id: links.LinkEndId=..., integrated_observation_handling: IntegratedObservationPropertyHandling=...) -> ObservationDependentVariableSettings:
    """No documentation found."""

def integration_time_dependent_variable(observable_type: model_settings.ObservableType=...) -> ObservationDependentVariableSettings:
    """No documentation found."""

def retransmission_delays_dependent_variable(observable_type: model_settings.ObservableType=...) -> ObservationDependentVariableSettings:
    """No documentation found."""

def target_range_between_link_ends_dependent_variable(start_link_end_type: links.LinkEndType=..., end_link_end_type: links.LinkEndType=..., start_link_end_id: links.LinkEndId=..., end_link_end_id: links.LinkEndId=..., integrated_observation_handling: IntegratedObservationPropertyHandling=...) -> ObservationDependentVariableSettings:
    """No documentation found."""
interval_end: IntegratedObservationPropertyHandling
interval_start: IntegratedObservationPropertyHandling
interval_undefined: IntegratedObservationPropertyHandling