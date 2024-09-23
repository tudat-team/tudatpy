import numpy
import typing
__all__ = ['AccelerationModel', 'AerodynamicGuidance', 'ConstantThrustMagnitudeWrapper', 'CustomThrustMagnitudeWrapper', 'DependentVariablesInterface', 'HybridArcSimulationResults', 'HybridArcVariationalSimulationResults', 'MassRateModel', 'MultiArcSimulationResults', 'MultiArcVariationalSimulationResults', 'PropagationTerminationDetails', 'PropagationTerminationDetailsFromHybridCondition', 'PropagationTerminationReason', 'RotationalProperModeDampingResults', 'SimulationResults', 'SingleArcSimulationResults', 'SingleArcVariationalSimulationResults', 'ThrustMagnitudeWrapper', 'TorqueModel', 'combine_initial_states', 'get_damped_proper_mode_initial_rotational_state', 'get_generalized_acceleration_size', 'get_initial_rotational_state_of_body', 'get_initial_state_of_bodies', 'get_initial_state_of_body', 'get_single_integration_differential_equation_order', 'get_single_integration_size', 'get_state_of_bodies', 'nan_or_inf_detected_in_state', 'propagation_never_run', 'runtime_error_caught_in_propagation', 'termination_condition_reached', 'unknown_reason']

class AccelerationModel:

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class AerodynamicGuidance:
    angle_of_attack: float
    bank_angle: float
    sideslip_angle: float

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self) -> None:
        ...

    def updateGuidance(self, current_time: float) -> None:
        ...

class ConstantThrustMagnitudeWrapper(ThrustMagnitudeWrapper):
    constant_thrust_magnitude: float

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class CustomThrustMagnitudeWrapper(ThrustMagnitudeWrapper):

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class DependentVariablesInterface:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class HybridArcSimulationResults(SimulationResults):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def multi_arc_results(self) -> MultiArcSimulationResults:
        ...

    @property
    def single_arc_results(self) -> SingleArcSimulationResults:
        ...

class HybridArcVariationalSimulationResults(SimulationResults):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def multi_arc_results(self) -> MultiArcVariationalSimulationResults:
        ...

    @property
    def single_arc_results(self) -> SingleArcVariationalSimulationResults:
        ...

class MassRateModel:

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class MultiArcSimulationResults(SimulationResults):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def arc_end_times(self) -> list[float]:
        ...

    @property
    def arc_start_times(self) -> list[float]:
        ...

    @property
    def propagation_is_performed(self) -> bool:
        ...

    @property
    def single_arc_results(self) -> list[SingleArcSimulationResults]:
        ...

    @property
    def solution_is_cleared(self) -> bool:
        ...

class MultiArcVariationalSimulationResults(SimulationResults):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def arc_end_times(self) -> list[float]:
        ...

    @property
    def arc_start_times(self) -> list[float]:
        ...

    @property
    def propagation_is_performed(self) -> bool:
        ...

    @property
    def single_arc_results(self) -> list[SingleArcVariationalSimulationResults]:
        ...

    @property
    def solution_is_cleared(self) -> bool:
        ...

class PropagationTerminationDetails:
    """Object that provides information on the reason for the
	termination of the propagation.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def terminated_on_exact_condition(self) -> bool:
        """
        Boolean defining whether the propagation was terminated on an *exact* final condition,
        or once the propagation went *past* the determined final condition. The choice of behaviour is
        defined by the termination settings provided as input to the Simulator object. This variable only
        has a meaningful definition if the ``termination_reason`` has value ``termination_condition_reached``
        """

    @property
    def termination_reason(self) -> PropagationTerminationReason:
        """
        Enum defining the reason the propagation was terminated
        """

class PropagationTerminationDetailsFromHybridCondition(PropagationTerminationDetails):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def was_condition_met_when_stopping(self) -> list[bool]:
        ...

class PropagationTerminationReason:
    """Enumeration of types of termination of propagation.
	
	
	:member propagation_never_run:
	:member unknown_reason:
	:member termination_condition_reached:
	:member runtime_error_caught_in_propagation:
	:member nan_or_inf_detected_in_state:
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, PropagationTerminationReason]]
    nan_or_inf_detected_in_state: typing.ClassVar[PropagationTerminationReason]
    propagation_never_run: typing.ClassVar[PropagationTerminationReason]
    runtime_error_caught_in_propagation: typing.ClassVar[PropagationTerminationReason]
    termination_condition_reached: typing.ClassVar[PropagationTerminationReason]
    unknown_reason: typing.ClassVar[PropagationTerminationReason]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

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

class RotationalProperModeDampingResults:
    """
		"""
    damped_initial_state: numpy.ndarray
    forward_backward_dependent_variables: list[tuple[dict[float, numpy.ndarray], dict[float, numpy.ndarray]]]
    forward_backward_states: list[tuple[dict[float, numpy.ndarray], dict[float, numpy.ndarray]]]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class SimulationResults:
    """Base class for objects that store all results of a numerical propagation. Derived class are implemented for single-, multi- and hybrid-arc propagation of botj dynamics and variational equations
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def dependent_variable_interface(self) -> DependentVariablesInterface:
        ...

class SingleArcSimulationResults(SimulationResults):
    """Class that stores all the results (including logging data) of a single-arc propagation
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def cumulative_computation_time_history(self) -> dict[float, float]:
        """
        History of cumulative computation time in seconds needed during the propagation as key-value
        pairs. At each epoch (key) the computation time (value) in seconds is the total computation time
        used up to and including that time step. This includes the total time up to and including the current time step,
        since the beginning of the (single-arc) propagation.
        """

    @property
    def cumulative_number_of_function_evaluations_history(self) -> dict[float, int]:
        ...

    @property
    def dependent_variable_history(self) -> dict[float, numpy.ndarray]:
        """
        Dependent variables computed during the propagation as key-value pairs.
        They are returned in the order with the same order of the DependentVariableSaveSettings object as values,
        as value, with the epoch as key.
        """

    @property
    def dependent_variable_ids(self) -> dict[tuple[int, int], str]:
        """
        Key-value container with the starting entry of the dependent variables saved (key), along with associated ID (value).
        """

    @property
    def initial_and_final_times(self) -> tuple[float, float]:
        ...

    @property
    def integration_completed_successfully(self) -> bool:
        """
        Boolean defining whether the last propagation was finished
        successfully, as defined by the termination conditions, or if
        it was terminated prematurely (for instance due to an
        exception, or an Inf/NaN state entry being detected).
        """

    @property
    def ordered_dependent_variable_settings(self) -> dict[tuple[int, int], ...]:
        ...

    @property
    def processed_state_ids(self) -> dict[tuple[int, int], str]:
        """
        Key-value container with the starting entry of the states (key), along with associated ID (value).
        """

    @property
    def propagated_state_ids(self) -> dict[tuple[int, int], str]:
        """
        Key-value container with the starting entry of the states (key), along with associated ID (value).
        """

    @property
    def propagated_state_vector_length(self) -> int:
        ...

    @property
    def propagation_is_performed(self) -> bool:
        ...

    @property
    def solution_is_cleared(self) -> bool:
        ...

    @property
    def state_history(self) -> dict[float, numpy.ndarray]:
        """
        Numerical solution of the equations of motion as key-value pairs. The key denotes the epoch. The value contains the
        numerically calculated state at this epoch. For this function, the states are always converted to so-called
        'conventional' formulations (e.g. Cartesian states for translational dynamics), see `here <https://tudat-space.readthedocs.io/en/latest/_src_api/propagation_setup/settings/conventional_vs_propagated_coordinates.html>`_
        for details. For the history of the states that were actually propagated, use the ``unprocessed_state_history``.
        
        .. note:: The propagated state at each epoch contains the state types in the following order: Translational ( **C** ), Rotational ( **R** ), Mass ( **M** ), and Custom ( **C** ).
                  When propagating two bodies, an example of what the output state would look like is for instance:
                  [ **T** Body 1, **T** Body 2, **R** Body 1, **R** Body 2, **M** Body 1, **M** Body 2 ] The specifics can be retrieved using the :attr:`state_ids` attribute of this class
        
        .. note:: For propagation of translational dynamics using cowell
                  propagator, the conventional and propagated
                  coordinates are identical.
        """

    @property
    def termination_details(self) -> PropagationTerminationDetails:
        """
        Object describing the details of the event that triggered the termination of the last propagation.
        """

    @property
    def total_computation_time(self) -> float:
        ...

    @property
    def total_number_of_function_evaluations(self) -> float:
        ...

    @property
    def unordered_dependent_variable_settings(self) -> list[...]:
        ...

    @property
    def unprocessed_state_history(self) -> dict[float, numpy.ndarray]:
        """
        Numerical solution of the equations of motion as key-value pairs, without any processing applied. The key denotes the epoch. The value contains the
        numerically calculated state at this epoch. This attribute contains the states of the propagated bodies expressed in the
        "raw" form in which the propagation took place. For instance, when using a Gauss-Kepler propagation scheme, this
        attribute will contain the numerically propagated Keplerian elements at each time epoch
        """

class SingleArcVariationalSimulationResults(SimulationResults):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def dynamics_results(self) -> SingleArcSimulationResults:
        ...

    @property
    def sensitivity_matrix_history(self) -> dict[float, numpy.ndarray]:
        ...

    @property
    def state_transition_matrix_history(self) -> dict[float, numpy.ndarray]:
        ...

class ThrustMagnitudeWrapper:

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class TorqueModel:

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

def combine_initial_states(propagator_settings_per_type: dict[..., list[..., ...]]) -> numpy.ndarray:
    """Function to retrieve the initial state for a list of propagator settings.
	
	Function to retrieve the initial state for a list of propagator settings. This way, the initial state for
	different quantities to be propagated (e.g., translational state, rotational state, mass) are retrieved and
	organized in a single container.
	
	
	:param propagator_settings_per_type:
			Propagator settings where the type of propagation is reported as key and the respective list of propagator settings as value.
	:return:
			Vector of initial states, sorted in order of IntegratedStateType, and then in the order of the vector of SingleArcPropagatorSettings of given type.
	"""

def get_damped_proper_mode_initial_rotational_state(*args, **kwargs) -> RotationalProperModeDampingResults:
    ...

def get_generalized_acceleration_size(state_type: ...) -> int:
    ...

def get_initial_rotational_state_of_body(body_to_propagate: str, base_orientation: str, bodies: ..., initial_time: float) -> numpy.ndarray:
    ...

def get_initial_state_of_bodies(bodies_to_propagate: list[str], central_bodies: list[str], body_system: ..., initial_time: float) -> numpy.ndarray:
    ...

def get_initial_state_of_body(body_to_propagate: str, central_body: str, bodies: ..., initial_time: float) -> numpy.ndarray:
    ...

def get_single_integration_differential_equation_order(state_type: ...) -> int:
    ...

def get_single_integration_size(state_type: ...) -> int:
    ...

def get_state_of_bodies(bodies_to_propagate: list[str], central_bodies: list[str], body_system: ..., initial_time: float) -> numpy.ndarray:
    """Function to get the states of a set of bodies, with respect to some set of central bodies, at the requested time.
	
	Function to get the states of a set of bodies, with respect to some set of central bodies, at the requested time.
	
	
	:param bodies_to_propagate:
			List of bodies to be propagated.
	:param central_bodies:
			List of central bodies, each referred to a body being propagated (in the same order).
	:param bodies_to_propagate:
			System of bodies used in the propagation.
	:param initial_time:
			Initial time of the propagation.
	:return:
			Time at which the states should be retrieved.
	"""
nan_or_inf_detected_in_state: PropagationTerminationReason
propagation_never_run: PropagationTerminationReason
runtime_error_caught_in_propagation: PropagationTerminationReason
termination_condition_reached: PropagationTerminationReason
unknown_reason: PropagationTerminationReason