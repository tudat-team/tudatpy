import typing
__all__ = ['create_acceleration_models', 'create_mass_rate_models', 'create_torque_models']

def create_acceleration_models(body_system: ..., selected_acceleration_per_body: dict[str, dict[str, list[...]]], bodies_to_propagate: list[str], central_bodies: list[str]) -> dict[str, dict[str, list[..., 3, 1, 0, 3, ...]]]:
    """Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
	
	Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
	bodies and central bodies are provided as two separate lists with the same order.
	
	
	:param body_system:
			System of bodies to be used in the propagation.
	:param selected_acceleration_per_body:
			Key-value container, with key denoting the body undergoing the acceleration, and the value containing an additional key-value container, with the body exerting acceleration, and list of acceleration settings exerted by this body.
	:param bodies_to_propagate:
			List of bodies to propagate.
	:param central_bodies:
			List of central bodies, each referred to each propagated body in the same order.
	:return:
			Set of accelerations acting on the bodies to propagate, provided as dual key-value container, similar to the acceleration settings input, but now with ``AccelerationModel`` lists as inner value
	"""

def create_mass_rate_models(body_system: ..., selected_mass_rates_per_body: dict[str, list[...]], acceleration_models: dict[str, dict[str, list[..., 3, 1, 0, 3, ...]]]=None) -> dict[str, list[...]]:
    """Function to create a set of mass-rate models from associated settings.
	
	Function to create a set of mass-rate models from a map of bodies and mass-rate model types.
	If the mass-rate depends on any acceleration models (e.g. thrust), the acceleration
	models must be provided as an input.
	
	
	:param body_system:
			System of bodies to be used in the propagation.
	:param selected_mass_rates_per_body:
			Key-value container, with key denoting the body with changing mass, and the value containing a list of mass rate settings (in most cases, this list will have only a single entry)
	:param acceleration_models:
			Sorted list of acceleration models, as created by :func:`create_acceleration_models`
	:return:
			Set of mass-rate models, as key-value container, same as the settings input, with the difference that the rate settings objects have been processed into the associated objects calculating the actual mass-rate changes.
	"""

def create_torque_models(body_system: ..., selected_torque_per_body: dict[str, dict[str, list[...]]], bodies_to_propagate: list[str]) -> dict[str, dict[str, list[...]]]:
    """Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.
	
	Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
	bodies is provided as a list.
	
	
	:param body_system:
			System of bodies to be used in the propagation.
	:param selected_torque_per_body:
			Key-value container, with key denoting the body undergoing the torque, and the value containing an additional key-value container, with the body exerting torque, and list of torque settings exerted by this body.
	:param bodies_to_propagate:
			List of bodies to propagate.
	:return:
			Set of torques acting on the bodies to propagate, provided as dual key-value container, similar to the torque settings input, but now with ``TorqueModel`` lists as inner value
	"""