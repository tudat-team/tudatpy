import numpy
import pybind11_stubgen.typing_ext
from ...astro import time_representation
from ...dynamics import environment
from ...dynamics import propagation
from ...dynamics.propagation_setup import integrator
from ...dynamics.propagation_setup import propagator
import typing
__all__ = ['CombinedStateTransitionAndSensitivityMatrixInterface', 'DynamicsSimulator', 'HybridArcSimulator', 'MultiArcSimulator', 'SingleArcSimulator', 'SingleArcVariationalSimulator', 'VariationalSimulator', 'create_dynamics_simulator', 'create_variational_equations_solver']

class CombinedStateTransitionAndSensitivityMatrixInterface:
    """Class establishing an interface with the simulation's State Transition and Sensitivity Matrices.
    
    Class establishing an interface to the State Transition and Sensitivity Matrices.
    Instances of this class are instantiated automatically upon creation of :class:`~tudatpy.estimation.estimation_analysis.Estimator` objects,
    using the simulation information in the observation, propagation and integration settings that the :class:`~tudatpy.estimation.estimation_analysis.Estimator` instance is linked to."""

    def full_state_transition_sensitivity_at_epoch(self, time: float, add_central_body_dependency: bool=True, arc_defining_bodies: list[str]=[]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 Parameters
                 ----------
                 time : float
                     Time at which full concatenated state transition and sensitivity matrix are to be retrieved.
                 Returns
                 -------
                 numpy.ndarray[numpy.float64[m, n]]
                     Full concatenated state transition and sensitivity matrix at a given time.
        """

    def state_transition_sensitivity_at_epoch(self, time: float, add_central_body_dependency: bool=True, arc_defining_bodies: list[str]=[]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 Function to get the concatenated state transition and sensitivity matrix at a given time.
        
                 Function to get the concatenated state transition and sensitivity matrix at a given time.
                 Entries corresponding to parameters which are not active at the current arc are omitted.
        
        
                 Parameters
                 ----------
                 time : float
                     Time at which concatenated state transition and sensitivity matrix are to be retrieved.
                 Returns
                 -------
                 numpy.ndarray[numpy.float64[m, n]]
                     Concatenated state transition and sensitivity matrix at a given time.
        """

    @property
    def full_parameter_size(self) -> int:
        """
                 **read-only**
        
                 Full amount of parameters w.r.t. which partials have been set up via State Transition and Sensitivity Matrices.
        
                 :type: int
        """

    @property
    def sensitivity_size(self) -> int:
        """
                 **read-only**
        
                 Number of columns in the sensitivity matrix.
        
                 :type: int
        """

    @property
    def state_transition_size(self) -> int:
        """
                 **read-only**
        
                 Size of the (square) state transition matrix.
        
                 :type: int
        """

class DynamicsSimulator:
    """Base class for propagation of dynamics (with derived classes implementing single-, multi- or hybrid-arc."""

class HybridArcSimulator(DynamicsSimulator):
    """Class for propagation of hybrid-arc dynamics.
    
    Class for propagation of hybrid-arc dynamics, from propagation settings and environment models, typically
    instantiated using :func:`~create_dynamics_simulator` function. See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagating_dynamics.html>`_ for more details."""

    def integrate_equations_of_motion(self, initial_state: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> None:
        """
         Function to reintegrate the equations of motion with a new initial state.
        
         Parameters
         ----------
         initial_state : numpy.ndarray
                New initial state from which the dynamics is to be propagated, consisting of concatenated initial state single-arc portion of dynamics, followed by concatenated initial states of the constituent arcs of the multi-arc portion
        """

    @property
    def multi_arc_simulator(self) -> MultiArcSimulator:
        """
                 **read-only**
        
                 Object used to propagate the multi-arc portion of the hybrid-arc dynamics
        
                 :type: MultiArcSimulator
        """

    @property
    def propagation_results(self) -> propagation.HybridArcSimulationResults:
        """
                 **read-only**
        
                 Object containing all the results of the numerical propagation, stored
                 in a single wrapper object
        
        
                 :type: HybridArcSimulationResults
        """

    @property
    def single_arc_simulator(self) -> SingleArcSimulator:
        """
                 **read-only**
        
                 Object used to propagate the single-arc portion of the hybrid-arc dynamics
        
        
                 :type: SingleArcSimulator
        """

class MultiArcSimulator(DynamicsSimulator):
    """Class for propagation of multi-arc dynamics.
    
    Class for propagation of multi-arc dynamics, from propagation settings and environment models, typically
    instantiated using :func:`~create_dynamics_simulator` function. See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagating_dynamics.html>`_ for more details."""

    def integrate_equations_of_motion(self, initial_state: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> None:
        """
         Function to reintegrate the equations of motion with a new initial state.
        
         Parameters
         ----------
         initial_state : numpy.ndarray
                New initial state from which the dynamics is to be propagated, consisting of concatenated initial states of each arc
        """

    @property
    def propagation_results(self) -> propagation.MultiArcSimulationResults:
        """
                 **read-only**
        
                 Object containing all the results of the numerical propagation, stored
                 in a single wrapper object
        
        
                 :type: MultiArcSimulationResults
        """

    @property
    def single_arc_simulators(self) -> list[SingleArcSimulator]:
        """
                 **read-only**
        
                 List of single-arc simulators, that were used to propagate the dynamics on the constituent single arcs
        
        
                 :type: list[SingleArcSimulator]
        """

class SingleArcSimulator(DynamicsSimulator):
    """Class for propagation of single arc dynamics.
    
    Class for propagation of single arc dynamics, from propagation settings and environment models, typically
    instantiated using :func:`~create_dynamics_simulator` function. See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagating_dynamics.html>`_ for more details."""

    def integrate_equations_of_motion(self, initial_state: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> None:
        """
         Function to reintegrate the equations of motion with a new initial state.
        
         Parameters
         ----------
         initial_state : numpy.ndarray
                New initial state from which the dynamics is to be propagated
        """

    @property
    def bodies(self) -> environment.SystemOfBodies:
        """
                 **read-only**
        
                 Object storing the set of bodies that comprise the physical environment.
        
                 :type: SystemOfBodies
        """

    @property
    def cumulative_computation_time_history(self) -> dict[float, float]:
        """
                 **read-only**
        
                 Shorthand for propagation_results.cumulative_computation_time_history
        
                 :type: dict[float, float]
        """

    @property
    def cumulative_number_of_function_evaluations(self) -> dict[float, int]:
        """
                 **read-only**
        
                 Shorthand for propagation_results.cumulative_number_of_function_evaluations
        
                 :type: dict[float, int]
        """

    @property
    def dependent_variable_history(self) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 Shorthand for propagation_results.dependent_variable_history
        
                 :type: dict[float, numpy.ndarray]
        """

    @property
    def dependent_variable_history_time_object(self) -> dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 Shorthand for propagation_results.dependent_variable_history_time_object
        
                 :type: dict[float, numpy.ndarray]
        """

    @property
    def environment_updater(self) -> ...:
        """
                 **read-only**
        
                 Object used in the propagation to update the environment, it uses the current time and numerically calculated state
                 to update the translational state, rotational state, flight conditions, etc. of all bodies in the simulation to be
                 consistent with this time and state.  Typically, this class is NOT used directly by users, but can be useful in specific situations.
        
        
                 :type: EnvironmentUpdater
        """

    @property
    def integration_completed_successfully(self) -> bool:
        """
                 **read-only**
        
                 Shorthand for propagation_results.integration_completed_successfully
        
                 :type: bool
        """

    @property
    def propagation_results(self) -> propagation.SimulationResults:
        """
                 **read-only**
        
                 Object containing all the results of the numerical propagation, stored
                 in a single wrapper object.
        
        
                 :type: SingleArcSimulationResults
        """

    @property
    def propagation_termination_details(self) -> propagation.PropagationTerminationDetails:
        """
                 **read-only**
        
                 Shorthand for propagation_results.termination_details
        
                 :type: PropagationTerminationDetails
        """

    @property
    def propagator_settings(self) -> propagation_setup.propagator.SingleArcPropagatorSettings:
        """
                 **read-only**
        
                 Settings for the propagation used to initialize the simulator
        
                 :type: SingleArcPropagatorProcessingSettings
        """

    @property
    def state_derivative_function(self) -> typing.Callable[[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]:
        """
                 **read-only**
        
                 Function that performs a single state derivative function evaluation. This function takes the numerically propagated
                 state, and current independent variable (time) as input, and returns the derivative of the state that is then used
                 by the numerical integration routine. Typically, this function is NOT used directly by users.
        
        
                 :type: Callable[[float, numpy.ndarray], numpy.ndarray]
        """

    @property
    def state_history(self) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 Shorthand for propagation_results.state_history
        
                 :type: dict[float, numpy.ndarray]
        """

    @property
    def state_history_time_object(self) -> dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 Shorthand for propagation_results.state_history_time_object
        
                 :type: dict[float, numpy.ndarray]
        """

    @property
    def unprocessed_state_history(self) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 Shorthand for propagation_results.unprocessed_state_history
        
                 :type: dict[float, numpy.ndarray]
        """

    @property
    def unprocessed_state_history_time_object(self) -> dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 Shorthand for propagation_results.unprocessed_state_history_time_object
        
                 :type: dict[float, numpy.ndarray]
        """

class SingleArcVariationalSimulator:
    """Class for single arc variational equations propagation."""

    def __init__(self, bodies: environment.SystemOfBodies, integrator_settings: propagation_setup.integrator.IntegratorSettings, propagator_settings: propagation_setup.propagator.PropagatorSettings, estimated_parameters: ..., integrate_equations_concurrently: bool=True, variational_only_integrator_settings: ...=None, clear_numerical_solutions: bool=False, integrate_on_creation: bool=True, set_integrated_result: bool=False, print_dependent_variable_data: bool=True) -> None:
        """
                 Class constructor.
        
                 Constructor through which the user can create instances of this class.
                 Defines environment, propagation and integrations models, as well as a number of settings related
                 to the (estimatable) parameters, w.r.t. which the variational equations are defined.
        
                 .. note:: When using default settings, creating an object of
                           this type automatically triggers the propagation
        
        
                 Parameters
                 ----------
                 bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
                     Object defining the physical environment, with all
                     properties of artificial and natural bodies.
        
                 integrator_settings : :class:`~tudatpy.dynamics.propagation_setup.integrator.IntegratorSettings`
                     Settings to create the numerical integrator that is to be used for the integration of the equations of motion.
        
                 propagator_settings : :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`
                     Settings to create the propagator that is to be used for the propagation of the dynamics.
        
                 estimated_parameters : :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`
                     Object defining a consolidated set of (estimatable) parameters (w.r.t. variational equations are defined),
                     linked to the environment and acceleration settings of the simulation.
        
                 integrate_equations_concurrently : bool, default = True
                     Boolean defining whether equations of motion and variational equations are to be propagated concurrently
                     (if true) or sequentially (of false).
        
                 variational_only_integrator_settings : :class:`~tudatpy.dynamics.propagation_setup.integrator.IntegratorSettings`, default = []
                     Settings to create the numerical integrator that is to be used for integration the variational equations.
                     If none is given (default), the numerical integration settings are taken to be the same as the ones applied
                     in the integration of the equations of motions (specified by the `integrator_settings` parameter).
        
                 clear_numerical_solutions : bool, default = False
                     Boolean to determine whether to clear the raw numerical solution member variables
                     and to reset the state transition interface after propagation.
        
                 integrate_on_creation : bool, default = True
                     Boolean defining whether the propagation should be performed immediately (default), or at a later time
                     (when calling the :func:`integrate_full_equations` or :func:`integrate_equations_of_motion_only` member function).
        
                 set_integrated_result : bool, default = True
                     Boolean to determine whether to automatically use the integrated results to set ephemerides for the
                     propagated bodies.
        """

    def integrate_equations_of_motion_only(self, initial_states: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]) -> None:
        """
                 Function to trigger the integration of the (regular) equations of motion.
        
        
                 Function to trigger the integration only of the (regular) equations of motion, resulting in a `state_history`.
                 This step does not yet use variational dynamics. In order to also solve the variational equations,
                 use the :func:`integrate_full_equations` member function.
        
                 Returns
                 -------
                 None
                     Creates / modifies the `state_history` property of the :class:`~tudatpy.dynamics.simulator.SingleArcVariationalSimulator` object.
        """

    def integrate_full_equations(self, initial_states: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)], integrate_equations_concurrently: bool=True) -> None:
        """
                 Function to trigger the integration of variational and dynamical equations (equations of motion).
        
        
                 Function to trigger the integration of the (regular) equations of motion as well as the variational equations,
                 solving for `state_history` and `variational_equations_history`
                 (in its two components `state_transition_matrix_history` & `sensitivity_matrix_history`).
        
        
                 Parameters
                 ----------
                 initial_states : numpy.ndarray([m, 1])
                     Initial state to be used for the parameters in the equations of motion.
        
                 integrate_equations_concurrently : bool, default = True
                     Boolean defining whether equations of motion and variational equations are to be propagated concurrently
                     (if true) or sequentially (of false).
        
                 Returns
                 -------
                 None
                     Creates / modifies the properties of the VariationalSolver object.
        """

    @property
    def dynamics_simulator(self) -> SingleArcSimulator:
        """
                 **read-only**
        
                 Simulator object containing all functionality for solving of the (regular) equations of motion.
        
        
                 :type: :class:`~tudatpy.dynamics.simulator.SingleArcSimulator`
        """

    @property
    def parameter_vector(self) -> ...:
        """
                 Consolidated set of (estimatable) parameters
                 w.r.t. the variational dynamics in the Variational Simulator are defined.
        
        
                 :type: :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`
        """

    @parameter_vector.setter
    def parameter_vector(self, arg1: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)], arg2: bool) -> None:
        ...

    @property
    def sensitivity_matrix_history(self) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]:
        """
                 **read-only**
        
                 Sensitivity matrix history, given as epoch with propagation epochs as keys.
                 This is (alongside the `state_transition_matrix_history`) the solution of the variational equations.
        
        
                 :type: dict[float, numpy.ndarray]
        """

    @property
    def state_history(self) -> dict[time_representation.Time, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 **read-only**
        
                 State history, given as epoch with propagation epochs as keys.
                 This is the solution of the (propagated) equations of motion, describing the states along which
                 the variational dynamics are solved.
        
        
                 :type: dict[float, numpy.ndarray]
        """

    @property
    def state_transition_matrix_history(self) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]:
        """
                 **read-only**
        
                 State transition matrix history, given as epoch with propagation epochs as keys.
                 This is (alongside the `sensitivity_matrix_history`) the solution of the variational equations.
        
        
                 :type: dict[float, numpy.ndarray]
        """

    @property
    def variational_equations_history(self) -> list[dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]]:
        """
                 **read-only**
        
                 List containing the solution of the variational equations, i.e. the
                 state transition matrix history (first entry) and sensitivity matrix history (second vector entry).
        
        
                 :type: list[ dict[float, numpy.ndarray] ]
        """

class VariationalSimulator:
    '''Base class for variational equations propagation.
    
    Base class for variational equations propagation.
    Derived classes :class:`~SingleArcVariationalSimulator`, :class:`~MultiArcVariationalSimulator` and
    :class:`~HybridArcVariationalSimulator` implement single-, multi- and hybrid-arc functionality, respectively."'''

def create_dynamics_simulator(bodies: environment.SystemOfBodies, propagator_settings: propagation_setup.propagator.PropagatorSettings, simulate_dynamics_on_creation: bool=True) -> DynamicsSimulator:
    """Function to create object that propagates the dynamics.
    
    Function to create object that propagates the dynamics, as specified by propagator settings, and the physical environment.
    Depending on the specific input type (e.g. which function from the :ref:`propagator` module was used),
    a single-, multi- or hybrid-arc simulator is created. The environment is typically created by the :func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies`
    function.
    
    .. note::
    
        When using default settings, calling this function will automatically propagate the dynamics.
    
    
    Parameters
    ----------
    bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
        Object defining the physical environment, with all
        properties of artificial and natural bodies.
    
    propagator_settings : :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`
        Settings to be used for the numerical propagation (dynamics type, termination conditions, integrator, etc.)
    
    simulate_dynamics_on_creation : bool, default=True
        Boolean defining whether to propagate the dynamics upon creation of the Simulator. If false, the dynamics c
        can be propagated at a later time by calling the :func:`~tudatpy.dynamics.simulator.DynamicsSimulator.integrate_equations_of_motion` function
    
    Returns
    -------
    :class:`~tudatpy.dynamics.simulator.DynamicsSimulator`
        Object that propagates the dynamics, and processes the results.
        Depending on the ``propagator_settings``, this object can be a single-, multi- or hybrid-arc simulator."""

def create_variational_equations_solver(bodies: environment.SystemOfBodies, propagator_settings: propagation_setup.propagator.PropagatorSettings, parameters_to_estimate: ..., simulate_dynamics_on_creation: bool=True) -> VariationalSimulator:
    """Function to create object that propagates the dynamics and variational equations.
    
    Function to create object that propagates the dynamics and variational equations, as specified by propagator settings, the physical environment, and a set of parameters for which to compute the partials.
    Depending on the specific input type (e.g. which function from the :ref:`propagator` module was used to define the propagator settings),
    a single-, multi- or hybrid-arc variational solver is created. The environment is typically created by the :func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies`
    function. When using default settings, calling this function will automatically propagate the dynamics.
    
    Parameters
    ----------
    bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
     Object defining the physical environment, with all
     properties of artificial and natural bodies.
    
    propagator_settings : :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`
     Settings to be used for the numerical propagation (dynamics type, termination conditions, integrator, etc.)
    
    parameters_to_estimate : :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`
     Object defining a consolidated set of (estimatable) parameters (w.r.t. variational equations are defined),
     linked to the environment and acceleration settings of the simulation.
    
    simulate_dynamics_on_creation : bool, default=True
     Boolean defining whether to propagate the dynamics upon creation of the Simulator. If false, the dynamics
     can be propagated at a later time by calling the :func:`~tudatpy.dynamics.DynamicsSimulator.integrate_equations_of_motion` function
    
    Returns
    -------
    :class:`~tudatpy.dynamics.simulator.VariationalSimulator`
     Object that propagates the dynamics, and processes the results."""