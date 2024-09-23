import numpy
import typing
__all__ = ['DynamicsSimulator', 'Estimator', 'HybridArcDynamicsSimulator', 'IntegratorSettings', 'MultiArcDynamicsSimulator', 'SingleArcSimulator', 'SingleArcVariationalSimulator', 'Time', 'create_dynamics_simulator', 'create_variational_equations_solver', 'get_integrated_type_and_body_list', 'get_single_integration_size']

class DynamicsSimulator:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class Estimator:
    """Class for consolidating all estimation functionality.
	
	Class for consolidating all functionality required to perform an estimation.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @staticmethod
    def compute_covariance(*args, **kwargs) -> ...:
        """
        Function to perform a covariance analysis for the given observations and parameters
        
        
        	Function to perform a covariance analysis for the given observations and parameters. The observations are provided through the
        	``covariance_analysis_input`` input, as are the weights :math:`\\mathbf{W}` and inverse a priori covariance :math:`(\\mathbf{P}_{0})^{-1}`.
        	Calling this function uses the environment and propagator settings provided to the constructor of this `Estimator` class to simulate
        	the dynamics of any relevant bodies for the observations (and associated variational equations). The observations are then
        	computed using the observation models created by the settings provided to the constructor of this `Estimator` class, as is the
        	associated design matrix :math:`\\mathbf{H}`. This function then produces the covariance :math:`\\mathbf{P}` (omitting the normalization used
        	internally for numerical stability)
        
        	.. math::
        	   \\mathbf{P}=\\left(\\mathbf{H}^{T}\\mathbf{W}\\mathbf{H}+(\\mathbf{P}_{0})^{-1}\right)^{-1}
        
        	Note that, although the actual observations are formally not required for a covariance analysis, all additional data (e.g. observation time, type, link ends, etc.)
        	are. And, as such, the ``covariance_analysis_input`` does require the full set of observations and associated information, for consistency purposes (e.g., same input as
        	``perform_estimation`` function) .
        
        
        	:param covariance_analysis_input:
        		Object consolidating all relevant settings for the covariance analysis
        		This includes foremost the simulated observations, as well as a priori information about the estimatable parameters
        
        	:return:
        		Object containing all outputs from the estimation process.
        """

    @staticmethod
    def perform_estimation(*args, **kwargs) -> ...:
        """
        Function to trigger the parameter estimation.
        
        
        	Function to trigger the parameter estimation. Much of the process and requirements are similar to those described in the
        	:func:`~tudatpy.numerical_simulation.Estimator.compute_covariance` function. This function uses an iterative least-squares
        	estimate process to fit the data (inside ``estimation_input``) to the model defined by the inputs to the ``Estimator`` constructor.s
        
        
        	:param estimation_input:
        		Object consolidating all relevant settings for the estimation
        		This includes foremost the simulated observations, as well as a priori information about the estimatable parameters and convergence criteria for the least squares estimation.
        
        	:return:
        		Object containing all outputs from the estimation process.
        """

    def __init__(self, bodies: ..., estimated_parameters: ..., observation_settings: list[...], propagator_settings: ..., integrate_on_creation: bool=True) -> None:
        """
        Class constructor.
        
        	Constructor through which the user can create instances of this class.
        	Defines environment, propagation and integrations models, as well as a number of settings related
        	to the estimatable parameters and observation settings.
        
        	.. note:: When using default settings, creating an object of
        	          this type automatically triggers the propagation
        
        
        	:param bodies:
        		Object defining the physical environment, with all
        		properties of artificial and natural bodies.
        
        	:param estimated_parameters:
        		Object defining a consolidated set of estimatable parameters,
        		linked to the environment and acceleration settings of the simulation.
        
        	:param observation_settings:
        		List of settings objects, each object defining the observation model settings for one
        		combination of observable and link geometry that is to be simulated.
        
        	:param integrator_settings:
        		Settings to create the numerical integrator that is to be
        		used for the integration of the equations of motion
        
        	:param propagator_settings:
        		Settings to create the propagator that is to be
        		used for the propagation of dynamics
        
        	:param integrate_on_creation:
        		Boolean defining whether the propagation should be
        		performed immediately (default), or at a later time
        		(when calling the :func:`perform_estimation` member function.
        """

    @property
    def observation_managers(self) -> dict[..., ..., ...]:
        """
        Observation managers contained in the Estimator object. A single observation manager can simulate observations and
        calculate observation partials for all link ends involved in the given observable type.
        """

    @property
    def observation_simulators(self) -> list[..., ...]:
        """
        Observation simulators contained in the Estimator object. A single observation simulator hosts
        the functionality for simulating a given observable over the defined link geometry.
        """

    @property
    def state_transition_interface(self) -> ...:
        """
        State transition and sensitivity matrix interface, setting the variational equations/dynamics in the
        Estimator object.
        """

    @property
    def variational_solver(self) -> ...:
        """
        Variational equations solver, which is used to manage and execute the numerical integration of
        equations of motion and variational equations/dynamics in the Estimator object.
        """

class HybridArcDynamicsSimulator(DynamicsSimulator):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def propagation_results(self) -> ...:
        ...

class IntegratorSettings:
    """Functional base class to define settings for integrators.
	
	Class to define settings for numerical integrators, for instance for use in numerical integration of equations of motion/variational equations. This class can be used for simple integrators such as fixed step RK and Euler. Integrators that require more settings to define have their own derived class.
	"""
    initial_time: float

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class MultiArcDynamicsSimulator(DynamicsSimulator):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def propagation_results(self) -> ...:
        ...

class SingleArcSimulator(DynamicsSimulator):
    """Class for consolidating single arc dynamics simulation functionality.
	
	Class for consolidating all functionality required to perform single arc dynamics simulations.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, bodies: ..., integrator_settings: IntegratorSettings, propagator_settings: ..., are_equations_of_motion_to_be_integrated: bool=True, clear_numerical_solutions: bool=False, set_integrated_result: bool=False, print_number_of_function_evaluations: bool=False, print_dependent_variable_data: bool=True, print_state_data: bool=True) -> None:
        ...

    @property
    def cumulative_computation_time_history(self) -> dict[float, float]:
        ...

    @property
    def cumulative_number_of_function_evaluations(self) -> dict[float, int]:
        ...

    @property
    def dependent_variable_history(self) -> dict[float, numpy.ndarray]:
        ...

    @property
    def environment_updater(self) -> ...:
        """
        # Object used in the propagation to update the environment, it uses the current time and numerically calculated state
        to update the translational state, rotational state, flight conditions, etc. of all bodies in the simulation to be
        consistent with this time and state.  Typically, this class is NOT used directly by users.
        """

    @property
    def integration_completed_successfully(self) -> bool:
        ...

    @property
    def integrator_settings(self) -> IntegratorSettings:
        """
        Settings to create the numerical integrator that is to be used
        for the integration of the equations of motion
        """

    @property
    def propagation_results(self) -> ...:
        """
        This function retrieves all the results of the numerical propagation, stored
        in a single wrapper object
        """

    @property
    def propagation_termination_details(self) -> ...:
        ...

    @property
    def state_derivative_function(self) -> typing.Callable[[float, numpy.ndarray], numpy.ndarray]:
        """
        Function that performs a single state derivative function evaluation. This function takes the numerically propagated
        state, and current independent variable (time) as input, and returns the derivative of the state that is then used
        by the numerical integration routine. Typically, this function is NOT used directly by users.
        """

    @property
    def state_history(self) -> dict[float, numpy.ndarray]:
        ...

    @property
    def unprocessed_state_history(self) -> dict[float, numpy.ndarray]:
        ...

class SingleArcVariationalSimulator:
    """Class for consolidating single arc variational dynamics functionality.
	
	Class for consolidating all functionality required to perform single arc variational dynamics simulations.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, bodies: ..., integrator_settings: IntegratorSettings, propagator_settings: ..., estimated_parameters: ..., integrate_equations_concurrently: bool=True, variational_only_integrator_settings: IntegratorSettings=None, clear_numerical_solutions: bool=False, integrate_on_creation: bool=True, set_integrated_result: bool=False, print_dependent_variable_data: bool=True) -> None:
        """
        Class constructor.
        
        	Constructor through which the user can create instances of this class.
        	Defines environment, propagation and integrations models, as well as a number of settings related
        	to the (estimatable) parameters, w.r.t. which the variational equations are defined.
        
        	.. note:: When using default settings, creating an object of
        	          this type automatically triggers the propagation
        
        
        	:param bodies:
        		Object defining the physical environment, with all
        		properties of artificial and natural bodies.
        
        	:param integrator_settings:
        		Settings to create the numerical integrator that is to be used for the integration of the equations of motion.
        
        	:param propagator_settings:
        		Settings to create the propagator that is to be used for the propagation of the dynamics.
        
        	:param estimated_parameters:
        		Object defining a consolidated set of (estimatable) parameters (w.r.t. variational equations are defined),
        		linked to the environment and acceleration settings of the simulation.
        
        	:param integrate_equations_concurrently:
        		Boolean defining whether equations of motion and variational equations are to be propagated concurrently
        		(if true) or sequentially (of false).
        
        	:param variational_only_integrator_settings:
        		Settings to create the numerical integrator that is to be used for integration the variational equations.
        		If none is given (default), the numerical integration settings are taken to be the same as the ones applied
        		in the integration of the equations of motions (specified by the `integrator_settings` parameter).
        
        	:param clear_numerical_solutions:
        		Boolean to determine whether to clear the raw numerical solution member variables
        		and to reset the state transition interface after propagation.
        
        	:param integrate_on_creation:
        		Boolean defining whether the propagation should be performed immediately (default), or at a later time
        		(when calling the :func:`integrate_full_equations` or :func:`integrate_equations_of_motion_only` member function).
        
        	:param set_integrated_result:
        		Boolean to determine whether to automatically use the integrated results to set ephemerides for the
        		propagated bodies.
        """

    def integrate_equations_of_motion_only(self, initial_states: numpy.ndarray) -> None:
        """
        Function to trigger the integration of the (regular) equations of motion.
        
        
        	Function to trigger the integration only of the (regular) equations of motion, resulting in a `state_history`.
        	This step does not yet use variational dynamics. In order to also solve the variational equations,
        	use the :func:`integrate_full_equations` member function.
        
        	:return:
        		Creates / modifies the `state_history` property of the :class:`~tudatpy.numerical_simulation.SingleArcVariationalSolver` object.
        """

    def integrate_full_equations(self, initial_states: numpy.ndarray, integrate_equations_concurrently: bool=True) -> None:
        """
        Function to trigger the integration of variational and dynamical equations (equations of motion).
        
        
        	Function to trigger the integration of the (regular) equations of motion as well as the variational equations,
        	solving for `state_history` and `variational_equations_history`
        	(in its two components `state_transition_matrix_history` & `sensitivity_matrix_history`).
        
        
        	:param initial_states:
        		Initial state to be used for the parameters in the equations of motion.
        
        	:param integrate_equations_concurrently:
        		Boolean defining whether equations of motion and variational equations are to be propagated concurrently
        		(if true) or sequentially (of false).
        
        	:return:
        		Creates / modifies the properties of the VariationalSolver object.
        """

    @property
    def dynamics_simulator(self) -> SingleArcSimulator:
        """
        Simulator object containing all functionality for solving of the (regular) equations of motion.
        """

    @property
    def parameter_vector(self) -> ...:
        """
        Consolidated set of (estimatable) parameters
        w.r.t. the variational dynamics in the Variational Simulator are defined.
        """

    @parameter_vector.setter
    def parameter_vector(self, arg1: numpy.ndarray, arg2: bool) -> None:
        ...

    @property
    def sensitivity_matrix_history(self) -> dict[float, numpy.ndarray]:
        """
        Sensitivity matrix history, given as epoch with propagation epochs as keys.
        This is (alongside the `state_transition_matrix_history`) the solution of the variational equations.
        """

    @property
    def state_history(self) -> dict[float, numpy.ndarray]:
        """
        State history, given as epoch with propagation epochs as keys.
        This is the solution of the (propagated) equations of motion, describing the states along which
        the variational dynamics are solved.
        """

    @property
    def state_transition_matrix_history(self) -> dict[float, numpy.ndarray]:
        """
        State transition matrix history, given as epoch with propagation epochs as keys.
        This is (alongside the `sensitivity_matrix_history`) the solution of the variational equations.
        """

    @property
    def variational_equations_history(self) -> list[dict[float, numpy.ndarray]]:
        """
        List containing the solution of the variational equations, i.e. the
        state transition matrix history (first entry) and sensitivity matrix history (second vector entry).
        """

class Time:
    """
		"""
    __hash__: typing.ClassVar[None] = None

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @typing.overload
    def __add__(self, arg0: Time) -> Time:
        ...

    @typing.overload
    def __add__(self, arg0: float) -> Time:
        ...

    @typing.overload
    def __eq__(self, arg0: Time) -> bool:
        ...

    @typing.overload
    def __eq__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __eq__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __ge__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __ge__(self, arg0: Time) -> bool:
        ...

    @typing.overload
    def __ge__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __gt__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __gt__(self, arg0: Time) -> bool:
        ...

    @typing.overload
    def __gt__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __iadd__(self, arg0: Time) -> None:
        ...

    @typing.overload
    def __iadd__(self, arg0: float) -> None:
        ...

    def __imul__(self, arg0: float) -> None:
        ...

    def __init__(self, full_periods: int, seconds_into_full_period: float) -> None:
        ...

    @typing.overload
    def __isub__(self, arg0: Time) -> None:
        ...

    @typing.overload
    def __isub__(self, arg0: float) -> None:
        ...

    def __itruediv__(self, arg0: float) -> None:
        ...

    @typing.overload
    def __le__(self, arg0: Time) -> bool:
        ...

    @typing.overload
    def __le__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __le__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __lt__(self, arg0: Time) -> bool:
        ...

    @typing.overload
    def __lt__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __lt__(self, arg0: float) -> bool:
        ...

    def __mul__(self, arg0: float) -> Time:
        ...

    @typing.overload
    def __ne__(self, arg0: Time) -> bool:
        ...

    @typing.overload
    def __ne__(self, arg0: float) -> bool:
        ...

    @typing.overload
    def __ne__(self, arg0: float) -> bool:
        ...

    def __radd__(self, arg0: float) -> Time:
        ...

    def __rmul__(self, arg0: float) -> Time:
        ...

    def __rsub__(self, arg0: float) -> Time:
        ...

    @typing.overload
    def __sub__(self, arg0: Time) -> Time:
        ...

    @typing.overload
    def __sub__(self, arg0: float) -> Time:
        ...

    def __truediv__(self, arg0: float) -> Time:
        ...

def create_dynamics_simulator(bodies: ..., propagator_settings: ..., simulate_dynamics_on_creation: bool=True) -> ...:
    """Function to create object that propagates the dynamics.
	
	Function to create object that propagates the dynamics, as specified by propagator settings, and the physical environment.
	Depending on the specific input type (e.g. which factory function from the :ref:`\\`\\`propagator\\`\\`` module was used),
	a single-, multi- or hybrid-arc simulator is created. The environment is typically created by the :func:`~tudatpy.numerical_simulation.environment_setup.create_system_of_bodies`
	function. When using default settings, calling this function will automatically propagate the dynamics.
	
	
	:param bodies:
			Object defining the physical environment, with all
			properties of artificial and natural bodies.
	
	:param propagator_settings:
			Settings to be used for the numerical propagation (dynamics type, termination conditions, integrator, etc.)
	
	:param simulate_dynamics_on_creation:
			Boolean defining whether to propagate the dynamics upon creation of the Simulator. If false, the dynamics c
			can be propagated at a later time by calling the :class:`~tudatpy.numerical_simulation.Simulator.integrate_equations_of_motion` function
	
	:return:
			Object that propagates the dynamics, and processes the results.
	"""

def create_variational_equations_solver(bodies: ..., propagator_settings: ..., parameters_to_estimate: ..., simulate_dynamics_on_creation: bool=True) -> ...:
    """Function to create object that propagates the dynamics.
	
	Function to create object that propagates the dynamics, as specified by propagator settings, and the physical environment.
	Depending on the specific input type (e.g. which factory function from the :ref:`\\`\\`propagator\\`\\`` module was used),
	a single-, multi- or hybrid-arc simulator is created. The environment is typically created by the :func:`~tudatpy.numerical_simulation.environment_setup.create_system_of_bodies`
	function. When using default settings, calling this function will automatically propagate the dynamics.
	
	
	:param bodies:
			Object defining the physical environment, with all
			properties of artificial and natural bodies.
	
	:param propagator_settings:
			Settings to be used for the numerical propagation (dynamics type, termination conditions, integrator, etc.)
	
	:param simulate_dynamics_on_creation:
			Boolean defining whether to propagate the dynamics upon creation of the Simulator. If false, the dynamics c
			can be propagated at a later time by calling the :class:`~tudatpy.numerical_simulation.Simulator.integrate_equations_of_motion` function
	
	:return:
			Object that propagates the dynamics, and processes the results.
	"""

def get_integrated_type_and_body_list(*args, **kwargs) -> dict[..., list[tuple[str, str, ...]]]:
    ...

def get_single_integration_size(state_type: ...) -> int:
    ...