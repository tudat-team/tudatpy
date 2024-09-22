import typing
import numpy
__all__ = ['create_observation_simulators', 'create_parameter_set', 'print_parameter_names', 'single_type_observation_collection']

def create_observation_simulators(observation_settings: list[...], bodies: ...) -> list[..., ...]:
    """Function for creating observation simulator objects.
	
	Factory function for creating observation simulator objects from observation settings.
	Note that each observation (i.e. combination of observable and link geometry) requires its own observation simulator object.
	
	
	:param observation_settings:
			List of settings objects, each object defining the observation model settings for one combination of observable and link geometry that is to be simulated.
	
	:param bodies:
			Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.
	
	:return:
			List of :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` objects, each object hosting the functionality for simulating one combination of observable type and link geometry.
	"""

def create_parameter_set(parameter_settings: list[...], bodies: ..., propagator_settings: ...=None, consider_parameters_names: list[...]=[]) -> ...:
    """Function for creating a consolidated set of estimatable parameters.
	
	Function for creating a consolidated parameter from the given estimatable parameter settings.
	The function checks for consistency between the parameter settings and the models contained in the simulation setup (given by the `bodies` & and `propagator_settings` parameters).
	
	
	:param parameter_settings:
			List of objects that define the settings for the parameters that are to be created. Each entry in this list is typically created by a call to a factory function in the :ref:`\\`\\`parameter\\`\\`` module
	
	:param bodies:
			Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.
	
	:param propagator_settings:
			Object containing the consolidated propagation settings of the simulation.
	
	:return:
			Instance of :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet` class, consolidating all estimatable parameters and simulation models.
	"""

def print_parameter_names(parameter_set: ...) -> None:
    ...

def single_type_observation_collection(observable_type: ..., link_ends: ..., observations_list: list[numpy.ndarray], times_list: list[float], reference_link_end: ..., ancilliary_settings: ...=None) -> ...:
    ...