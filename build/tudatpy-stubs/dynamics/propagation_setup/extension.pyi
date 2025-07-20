import typing
from ...dynamics import environment
from . import acceleration
from . import dependent_variable
from . import integrator
from . import mass_rate
from . import propagator
from . import thrust
from . import torque
__all__ = ['acceleration', 'create_acceleration_models', 'create_mass_rate_models', 'create_torque_models', 'dependent_variable', 'integrator', 'mass_rate', 'propagator', 'thrust', 'torque']

def create_acceleration_models(body_system: environment.SystemOfBodies, selected_acceleration_per_body: dict[str, dict[str, list[acceleration.AccelerationSettings]]], bodies_to_propagate: list[str], central_bodies: list[str]) -> dict[str, dict[str, list[..., 3, 1, 0, 3, ...]]]:
    """Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
    
    Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
    bodies and central bodies are provided as two separate lists with the same order.
    
    
    Parameters
    ----------
    body_system : SystemOfBodies
        System of bodies to be used in the propagation.
    selected_acceleration_per_body : Dict[str, Dict[str, List[AccelerationSettings]]]
        Key-value container, with key denoting the body undergoing the acceleration, and the value containing an additional key-value container, with the body exerting acceleration, and list of acceleration settings exerted by this body.
    bodies_to_propagate : list
        List of bodies to propagate.
    central_bodies : list
        List of central bodies, each referred to each propagated body in the same order.
    Returns
    -------
    AccelerationMap : dict[str, list[AccelerationModel]]
       Set of accelerations acting on the bodies to propagate, provided as dual key-value container (dictionary), similar to the acceleration settings input, but now with ``AccelerationModel`` lists as inner value
    
    
    
    
    
    Examples
    --------
    In this example, the acceleration model includes a spherical harmonic (degree and order 5) gravitational acceleration
    and aerodynamic acceleration exerted by the Earth, as well as point-mass gravity exerted by the Sun and the Moon.
    The variable ``accelerations_settings_vehicle`` denotes the list of bodies exerting accelerations and the types of
    accelerations, while the variable ``acceleration_settings`` associates this list with the body undergoing the
    acceleration (``"Vehicle"``).
    
    .. code-block:: python
    
       # Define bodies that are propagated.
       bodies_to_propagate = ["Vehicle"]
    
       # Define central bodies.
       central_bodies = ["Earth"]
    
       # Define accelerations acting on Vehicle
       accelerations_settings_vehicle = dict(
           Sun=[propagation_setup.acceleration.point_mass_gravity()],
           Moon=[propagation_setup.acceleration.point_mass_gravity()],
           Earth=[propagation_setup.acceleration.spherical_harmonic_gravity(5, 5),
                  propagation_setup.acceleration.aerodynamic()]
       )
    
       # Create global accelerations settings dictionary
       acceleration_settings = {"Vehicle": accelerations_settings_vehicle}
    
       # Create acceleration models
       acceleration_models = propagation_setup.create_acceleration_models(
           bodies, acceleration_settings,  bodies_to_propagate, central_bodies)"""

def create_mass_rate_models(body_system: environment.SystemOfBodies, selected_mass_rates_per_body: dict[str, list[mass_rate.MassRateModelSettings]], acceleration_models: dict[str, dict[str, list[..., 3, 1, 0, 3, ...]]]=None) -> dict[str, list[...]]:
    """Function to create a set of mass-rate models from associated settings.
    
    Function to create a set of mass-rate models from a map of bodies and mass-rate model types.
    If the mass-rate depends on any acceleration models (e.g. thrust), the acceleration
    models must be provided as an input.
    
    
    Parameters
    ----------
    body_system : SystemOfBodies
        System of bodies to be used in the propagation.
    selected_mass_rates_per_body : Dict[str, List[MassRateModelSettings]]
        Key-value container, with key denoting the body with changing mass, and the value containing a list of mass rate settings (in most cases, this list will have only a single entry)
    acceleration_models : dict[str, list[AccelerationModel]]
        Sorted list of acceleration models, as created by :func:`create_acceleration_models`
    Returns
    -------
    MassRateModelMap
        Set of mass-rate models, as key-value container, same as the settings input, with the difference that the rate settings objects have been processed into the associated objects calculating the actual mass-rate changes."""

def create_torque_models(body_system: environment.SystemOfBodies, selected_torque_per_body: dict[str, dict[str, list[torque.TorqueSettings]]], bodies_to_propagate: list[str]) -> dict[str, dict[str, list[...]]]:
    """Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
    
    Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
    bodies is provided as a list.
    
    
    Parameters
    ----------
    body_system : SystemOfBodies
        System of bodies to be used in the propagation.
    selected_torque_per_body : Dict[str, Dict[str, List[TorqueSettings]]]
        Key-value container, with key denoting the body undergoing the torque, and the value containing an additional key-value container, with the body exerting torque, and list of torque settings exerted by this body.
    bodies_to_propagate : list
        List of bodies to propagate.
    Returns
    -------
    TorqueModelMap
        Set of torques acting on the bodies to propagate, provided as dual key-value container, similar to the torque settings input, but now with ``TorqueModel`` lists as inner value
    
    
    
    
    
    Examples
    --------
    
    In this example, the following torques are exerted on the vehicle: spherical-harmonic gravitational torque
    (up to order 4 and degree 4) and aerodynamic torque exerted by the Earth, second-degree gravitational torque
    exerted by the Sun and the Moon.
    The variable ``torques_settings_vehicle`` denotes the list of bodies exerting torques and the types of
    torques, while the variable ``torque_settings`` associates this list with the body undergoing the
    torque.
    
    .. code-block:: python
    
      # Define bodies that are propagated
      bodies_to_propagate = ["Vehicle"]
    
      # Define torques per each exerting body
      torque_settings_vehicle = dict(
          Sun=[propagation_setup.torque.second_degree_gravitational()],
          Moon=[propagation_setup.torque.second_degree_gravitational()],
          Earth=[propagation_setup.torque.spherical_harmonic_gravitational(4, 4),
                 propagation_setup.torque.aerodynamic()]
      )
    
      # Create global torque settings dictionary
      torque_settings = {"Vehicle": torque_settings_vehicle}
    
      # Create torque models
      torque_models = propagation_setup.create_torque_models(
          bodies, torque_settings,  bodies_to_propagate )"""