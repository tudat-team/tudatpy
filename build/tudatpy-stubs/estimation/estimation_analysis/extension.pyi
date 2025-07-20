import numpy
import pybind11_stubgen.typing_ext
from ...astro import time_representation
from ...dynamics import environment
from ...dynamics import parameters
from ...dynamics import parameters_setup
from ...dynamics import propagation
from ...dynamics.propagation_setup import integrator
from ...dynamics.propagation_setup import propagator
from ...dynamics import simulator
from ...estimation.observable_models import observables_simulation
from ...estimation.observable_models_setup import links
from ...estimation.observable_models_setup import model_settings
from ...estimation import observations
import typing
__all__ = ['CovarianceAnalysisInput', 'CovarianceAnalysisOutput', 'EstimationConvergenceChecker', 'EstimationInput', 'EstimationOutput', 'Estimator', 'PodInput', 'PodOutput', 'create_best_fit_to_ephemeris', 'estimation_convergence_checker', 'propagate_covariance', 'propagate_covariance_rsw_split_output', 'propagate_covariance_split_output', 'propagate_formal_errors', 'propagate_formal_errors_rsw_split_output', 'propagate_formal_errors_split_output']

class CovarianceAnalysisInput:
    """Class for defining all specific inputs to a covariance analysis."""

    def __init__(self, observations_and_times: observations.ObservationCollection, inverse_apriori_covariance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]=..., consider_covariance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]=...) -> None:
        """
                 Class constructor.
        
                 Constructor through which the user can create instances of this class. Note that the weight are all initiated as 1.0, and the default settings of ``define_covariance_settings`` are used.
        
        
                 Parameters
                 ----------
                 observations_and_times : ObservationCollection
                     Total data structure of observations and associated times/link ends/type/etc.
                 inverse_apriori_covariance : numpy.ndarray[numpy.float64[m, n]], default = [ ]
                     A priori covariance matrix (unnormalized) of estimated parameters. This should be either a size 0x0 matrix (no a priori information), or a square matrix with the same size as the number of parameters that are considered
                 Returns
                 -------
                 :class:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisInput`
                     Instance of the :class:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisInput` class, defining the data and other settings to be used for the covariance analysis.
        """

    def define_covariance_settings(self, reintegrate_equations_on_first_iteration: bool=True, reintegrate_variational_equations: bool=True, save_design_matrix: bool=True, print_output_to_terminal: bool=True, limit_condition_number_for_warning: float=100000000.0) -> None:
        """
                 Function to define specific settings for covariance analysis process
        
                 Function to define specific settings for covariance analysis process
        
        
                 Parameters
                 ----------
                 reintegrate_equations : bool, default = True
                     Boolean denoting whether the dynamics and variational equations are to be reintegrated
                     or if existing values are to be used to perform first iteration.
        
                 reintegrate_variational_equations : bool, default = True
                     Boolean denoting whether the variational equations are to be reintegrated during estimation
                     (if this is set to False, and ``reintegrate_equations`` to true, only the dynamics are re-integrated)
        
                 save_design_matrix : bool, default = True
                     Boolean denoting whether to save the partials matrix (also called design matrix) :math:`\\mathbf{H}` in the output. Setting this to false makes the
                     :math:`\\mathbf{H}` matrix unavailable to the user, with the advantage of lower RAM usage.
        
                 print_output_to_terminal : bool, default = True
                     Boolean denoting whether to print covariance-analysis-specific output to the terminal when running the estimation.
        
                 Returns
                 -------
                 None
                     Function modifies the object in-place.
        """

    def set_constant_single_observable_and_link_end_vector_weight(self, observable_type: model_settings.ObservableType, link_ends: dict[links.LinkEndType, links.LinkEndId], weight: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]) -> None:
        """
        No documentation found.
        """

    def set_constant_single_observable_and_link_end_weight(self, observable_type: model_settings.ObservableType, link_ends: dict[links.LinkEndType, links.LinkEndId], weight: float) -> None:
        """
        No documentation found.
        """

    def set_constant_single_observable_vector_weight(self, observable_type: model_settings.ObservableType, weight: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]) -> None:
        """
        No documentation found.
        """

    def set_constant_single_observable_weight(self, observable_type: model_settings.ObservableType, weight: float) -> None:
        """
        No documentation found.
        """

    def set_constant_vector_weight_per_observable(self, weight_per_observable: dict[model_settings.ObservableType, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]) -> None:
        """
        No documentation found.
        """

    def set_constant_weight(self, weight: float) -> None:
        """
                 Function to set a constant weight matrix for all observables.
        
                 Function to set a constant weight matrix for all observables.
                 The weights are applied to all observations managed by the given PodInput object.
        
        
                 Parameters
                 ----------
                 constant_weight : float
                     Constant weight factor that is to be applied to all observations.
                 Returns
                 -------
                 None
                     Function modifies the object in-place.
        """

    def set_constant_weight_per_observable(self, weight_per_observable: dict[model_settings.ObservableType, float]) -> None:
        """
                 Function to set a constant weight matrix for a given type of observable.
        
                 Function to set a constant weight matrix for a given type of observable.
                 The weights are applied to all observations of the observable type specified by the `weight_per_observable` parameter.
        
        
                 Parameters
                 ----------
                 constant_weight : Dict[ :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservableType`, float ]
                     Constant weight factor that is to be applied to all observations.
                 Returns
                 -------
                 None
                     Function modifies the object in-place.
        """

    def set_total_single_observable_and_link_end_vector_weight(self, observable_type: model_settings.ObservableType, link_ends: dict[links.LinkEndType, links.LinkEndId], weight_vector: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]) -> None:
        """
        No documentation found.
        """

    def set_weights_from_observation_collection(self) -> None:
        """
        No documentation found.
        """

    @property
    def weight_matrix_diagonal(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
                 **read-only**
        
                 Complete diagonal of the weights matrix that is to be used
        
                 :type: numpy.ndarray[numpy.float64[n, 1]]
        """

    @weight_matrix_diagonal.setter
    def weight_matrix_diagonal(self, arg1: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]) -> None:
        ...

class CovarianceAnalysisOutput:
    """Class collecting all outputs from the covariance analysis process."""

    @property
    def consider_covariance_contribution(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
        No documentation found.
        """

    @property
    def consider_normalization_factors(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
        No documentation found.
        """

    @property
    def correlations(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 **read-only**
        
                 Correlation matrix of the estimation result. Entry :math:`i,j` is equal to :math:`P_{i,j}/(\\sigma_{i}\\sigma_{j})`
        
                 :type: numpy.ndarray[numpy.float64[m, m]]
        """

    @property
    def covariance(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 **read-only**
        
                 (Unnormalized) estimation covariance matrix :math:`\\mathbf{P}`.
        
                 :type: numpy.ndarray[numpy.float64[m, m]]
        """

    @property
    def design_matrix(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 **read-only**
        
                 Matrix of unnormalized partial derivatives :math:`\\mathbf{H}=\\frac{\\partial\\mathbf{h}}{\\partial\\mathbf{p}}`.
        
                 :type: numpy.ndarray[numpy.float64[m, n]]
        """

    @property
    def formal_errors(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
                 **read-only**
        
                 Formal error vector :math:`\\boldsymbol{\\sigma}` of the estimation result (e.g. square root of diagonal entries of covariance)s
        
                 :type: numpy.ndarray[numpy.float64[m, 1]]s
        """

    @property
    def inverse_covariance(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 **read-only**
        
                 (Unnormalized) inverse estimation covariance matrix :math:`\\mathbf{P}^{-1}`.
        
                 :type: numpy.ndarray[numpy.float64[m, m]]
        """

    @property
    def inverse_normalized_covariance(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 **read-only**
        
                 Normalized inverse estimation covariance matrix :math:`\\mathbf{\\tilde{P}}^{-1}`.
        
                 :type: numpy.ndarray[numpy.float64[m, m]]
        """

    @property
    def normalization_terms(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
                 **read-only**
        
                 Vector of normalization terms used for covariance and design matrix
        
                 :type: numpy.ndarray[numpy.float64[m, 1]]
        """

    @property
    def normalized_covariance(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 **read-only**
        
                 Normalized estimation covariance matrix :math:`\\mathbf{\\tilde{P}}`.
        
                 :type: numpy.ndarray[numpy.float64[m, m]]
        """

    @property
    def normalized_covariance_with_consider_parameters(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
        No documentation found.
        """

    @property
    def normalized_design_matrix(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 **read-only**
        
                 Matrix of normalized partial derivatives :math:`\\tilde{\\mathbf{H}}`.
        
                 :type: numpy.ndarray[numpy.float64[m, n]]
        """

    @property
    def normalized_design_matrix_consider_parameters(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
        No documentation found.
        """

    @property
    def unnormalized_covariance_with_consider_parameters(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
        No documentation found.
        """

    @property
    def weighted_design_matrix(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 **read-only**
        
                 Matrix of weighted partial derivatives, equal to :math:`\\mathbf{W}^{1/2}{\\mathbf{H}}`
        
                 :type: numpy.ndarray[numpy.float64[m, n]]
        """

    @property
    def weighted_normalized_design_matrix(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 **read-only**
        
                 Matrix of weighted, normalized partial derivatives, equal to :math:`\\mathbf{W}^{1/2}\\tilde{\\mathbf{H}}`
        
                 :type: numpy.ndarray[numpy.float64[m, n]]
        """

class EstimationConvergenceChecker:
    """Class defining the convergence criteria for an estimation.
    
    Class defining the convergence criteria for an estimation.
    The user typically creates instances of this class via the :func:`~tudatpy.estimation.estimation_analysis.estimation_convergence_checker` function."""

class EstimationInput(CovarianceAnalysisInput):
    """Class for defining all inputs to the estimation."""

    def __init__(self, observations_and_times: observations.ObservationCollection, inverse_apriori_covariance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]=..., convergence_checker: EstimationConvergenceChecker=..., consider_covariance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]=..., consider_parameters_deviations: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]=..., apply_final_parameter_correction: bool=True) -> None:
        """
                 Class constructor.
        
                 Constructor through which the user can create instances of this class.
        
        
                 Parameters
                 ----------
                 observations_and_times : ObservationCollection
                     Total data structure of observations and associated times/link ends/type/etc.
                 inverse_apriori_covariance : numpy.ndarray[numpy.float64[m, n]], default = [ ]
                     A priori covariance matrix (unnormalized) of estimated parameters. This should be either a size 0x0 matrix (no a priori information), or a square matrix with the same size as the number of parameters that are considered
                 convergence_checker : :class:`~tudatpy.estimation.estimation_analysis.EstimationConvergenceChecker`, default = :func:`~tudatpy.estimation.estimation_analysis.estimation_convergence_checker`
                     Object defining when the estimation is converged.
                 Returns
                 -------
                 :class:`~tudatpy.estimation.estimation_analysis.EstimationInput`
                     Instance of the :class:`~tudatpy.estimation.estimation_analysis.EstimationInput` class, defining the data and other settings to be used for the estimation.
        """

    def define_estimation_settings(self, reintegrate_equations_on_first_iteration: bool=True, reintegrate_variational_equations: bool=True, save_design_matrix: bool=True, print_output_to_terminal: bool=True, save_residuals_and_parameters_per_iteration: bool=True, save_state_history_per_iteration: bool=False, limit_condition_number_for_warning: float=100000000.0, condition_number_warning_each_iteration: bool=True) -> None:
        """
                 Function to define specific settings for the estimation process
        
                 Function to define specific settings for covariance analysis process
        
        
                 Parameters
                 ----------
                 reintegrate_equations_on_first_iteration : bool, default = True
                     Boolean denoting whether the dynamics and variational equations are to be reintegrated
                     or if existing values are to be used to perform first iteration.
        
                 reintegrate_variational_equations : bool, default = True
                     Boolean denoting whether the variational equations are to be reintegrated during estimation
                     (if this is set to False, and ``reintegrate_equations_on_first_iteration`` to true, only the dynamics are re-integrated)
        
                 save_design_matrix : bool, default = True
                     Boolean denoting whether to save the partials matrix (also called design matrix) :math:`\\mathbf{H}` in the output. Setting this to false makes the
                     :math:`\\mathbf{H}` matrix unavailable to the user, with the advantage of lower RAM usage.
        
                 print_output_to_terminal : bool, default = True
                     Boolean denoting whether to print covariance-analysis-specific output to the terminal when running the estimation.
        
                 save_residuals_and_parameters_per_iteration : bool, default = True
                     Boolean denoting whether the residuals and parameters from the each iteration are to be saved.
        
                 save_state_history_per_iteration : bool, default = False
                     Boolean denoting whether the state history and dependent variables are to be saved on each iteration.
        
                 Returns
                 -------
                 None
                     Function modifies the object in-place.
        """

class EstimationOutput(CovarianceAnalysisOutput):
    """Class collecting all outputs from the iterative estimation process."""

    @property
    def best_iteration(self) -> int:
        """
        No documentation found.
        """

    @property
    def final_parameters(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
        No documentation found.
        """

    @property
    def final_residuals(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
                 **read-only**
        
                 Vector of post-fit observation residuals, for the iteration with the lowest rms residuals.
        
                 :type: numpy.ndarray[numpy.float64[m, 1]]
        """

    @property
    def parameter_history(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 **read-only**
        
                 Parameter vectors, concatenated per iteration into a matrix. Column 0 contains pre-estimation values. The :math:`(i+1)^{th}` column has the residuals from the :math:`i^{th}` iteration.
        
                 :type: numpy.ndarray[numpy.float64[m, n]]
        """

    @property
    def residual_history(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 **read-only**
        
                 Residual vectors, concatenated per iteration into a matrix; the :math:`i^{th}` column has the residuals from the :math:`i^{th}` iteration.
        
                 :type: numpy.ndarray[numpy.float64[m, n]]
        """

    @property
    def simulation_results_per_iteration(self) -> list[propagation.SimulationResults]:
        """
                 **read-only**
        
                 List of complete numerical propagation results, with the :math:`i^{th}` entry of thee list thee results of the :math:`i^{th}` propagation
        
                 :type: list[SimulationResults]
        """

class Estimator:
    """Class for consolidating all estimation functionality.
    
    Class for consolidating all functionality required to perform an estimation."""

    def __init__(self, bodies: environment.SystemOfBodies, estimated_parameters: parameters.EstimatableParameterSet, observation_settings: list[model_settings.ObservationSettings], propagator_settings: propagation_setup.propagator.PropagatorSettings, integrate_on_creation: bool=True) -> None:
        """
                 Class constructor.
        
                 Constructor through which the user can create instances of this class.
                 Defines environment, propagation and integrations models, as well as a number of settings related
                 to the estimatable parameters and observation settings.
        
                 .. note:: When using default settings, creating an object of
                           this type automatically triggers the propagation
        
        
                 Parameters
                 ----------
                 bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
                     Object defining the physical environment, with all
                     properties of artificial and natural bodies.
        
                 estimated_parameters : :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`
                     Object defining a consolidated set of estimatable parameters,
                     linked to the environment and acceleration settings of the simulation.
        
                 observation_settings : :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings`
                     List of settings objects, each object defining the observation model settings for one
                     combination of observable and link geometry that is to be simulated.
        
                 integrator_settings : :class:`~tudatpy.dynamics.propagation_setup.integrator.IntegratorSettings`
                     Settings to create the numerical integrator that is to be
                     used for the integration of the equations of motion
        
                 propagator_settings : :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`
                     Settings to create the propagator that is to be
                     used for the propagation of dynamics
        
                 integrate_on_creation : bool, default = True
                     Boolean defining whether the propagation should be
                     performed immediately (default), or at a later time
                     (when calling the :func:`perform_estimation` member function.
        """

    def compute_covariance(self, covariance_analysis_input: ..., tudat: ...) -> ...:
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
                    \\mathbf{P}=\\left(\\mathbf{H}^{T}\\mathbf{W}\\mathbf{H}+(\\mathbf{P}_{0})^{-1}\\right)^{-1}
        
                 Note that, although the actual observations are formally not required for a covariance analysis, all additional data (e.g. observation time, type, link ends, etc.)
                 are. And, as such, the ``covariance_analysis_input`` does require the full set of observations and associated information, for consistency purposes (e.g., same input as
                 ``perform_estimation`` function) .
        
        
                 Parameters
                 ----------
                 covariance_analysis_input : :class:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisInput`
                     Object consolidating all relevant settings for the covariance analysis
                     This includes foremost the simulated observations, as well as a priori information about the estimatable parameters
        
                 Returns
                 -------
                 :class:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput`
                     Object containing all outputs from the estimation process.
        """

    def perform_estimation(self, estimation_input: ..., tudat: ...) -> ...:
        """
                 Function to trigger the parameter estimation.
        
        
                 Function to trigger the parameter estimation. Much of the process and requirements are similar to those described in the
                 :func:`~tudatpy.estimation.estimation_analysis.Estimator.compute_covariance` function. This function uses an iterative least-squares
                 estimate process to fit the data (inside ``estimation_input``) to the model defined by the inputs to the ``Estimator`` constructor.s
        
        
                 Parameters
                 ----------
                 estimation_input : :class:`~tudatpy.estimation.estimation_analysis.EstimationInput`
                     Object consolidating all relevant settings for the estimation
                     This includes foremost the simulated observations, as well as a priori information about the estimatable parameters and convergence criteria for the least squares estimation.
        
                 Returns
                 -------
                 :class:`~tudatpy.estimation.estimation_analysis.EstimationOutput`
                     Object containing all outputs from the estimation process.
        """

    @property
    def observation_managers(self) -> dict[model_settings.ObservableType, ..., ...]:
        """
                 **read-only**
        
                 Observation managers contained in the Estimator object. A single observation manager can simulate observations and
                 calculate observation partials for all link ends involved in the given observable type.
        
        
                 :type: dict[ :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservableType`, :class:`~tudatpy.estimation.observations.ObservationManager` ]
        """

    @property
    def observation_simulators(self) -> list[observables_simulation.ObservationSimulator]:
        """
                 **read-only**
        
                 Observation simulators contained in the Estimator object. A single observation simulator hosts
                 the functionality for simulating a given observable over the defined link geometry.
        
        
                 :type: list[ :class:`~tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator` ]
        """

    @property
    def state_transition_interface(self) -> simulator.CombinedStateTransitionAndSensitivityMatrixInterface:
        """
                 **read-only**
        
                 State transition and sensitivity matrix interface, setting the variational equations/dynamics in the
                 Estimator object.
        
        
                 :type: :class:`~tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface`
        """

    @property
    def variational_solver(self) -> simulator.VariationalSimulator:
        """
                 **read-only**
        
                 Variational equations solver, which is used to manage and execute the numerical integration of
                 equations of motion and variational equations/dynamics in the Estimator object.
        
                 :type: :class:`~tudatpy.dynamics.simulator.SingleArcVariationalSimulator`
        """

def create_best_fit_to_ephemeris(bodies: environment.SystemOfBodies, acceleration_models: dict[str, dict[str, list[propagation.AccelerationModel]]], observed_bodies: list[str], central_bodies: list[str], integrator_settings: propagation_setup.integrator.IntegratorSettings, initial_time: time_representation.Time, final_time: time_representation.Time, data_point_interval: time_representation.Time, additional_parameter_names: list[parameters_setup.EstimatableParameterSettings]=[], number_of_iterations: int=3, reintegrate_variational_equations: bool=True, results_print_frequency: float=0.0) -> EstimationOutput:
    """No documentation found."""

def estimation_convergence_checker(maximum_iterations: int=5, minimum_residual_change: float=0.0, minimum_residual: float=0.0, number_of_iterations_without_improvement: int=2) -> EstimationConvergenceChecker:
    """Function for creating an :class:`~tudatpy.estimation.estimation_analysis.EstimationConvergenceChecker` object.
    
    Function for creating an :class:`~tudatpy.estimation.estimation_analysis.EstimationConvergenceChecker` object, which is required for defining the convergence criteria of an estimation.
    
    
    Parameters
    ----------
    maximum_iterations : int, default = 5
        Maximum number of allowed iterations for estimation.
    minimum_residual_change : float, default = 0.0
        Minimum required change in residual between two iterations.
    minimum_residual : float, default = 0.0
        Minimum value of observation residual below which estimation is converged.
    number_of_iterations_without_improvement : int, default = 2
        Number of iterations without reduction of residual.
    Returns
    -------
    :class:`~tudatpy.estimation.estimation_analysis.EstimationConvergenceChecker`
        Instance of the :class:`~tudatpy.estimation.estimation_analysis.EstimationConvergenceChecker` class, defining the convergence criteria for an estimation."""

def propagate_covariance(initial_covariance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], state_transition_interface: simulator.CombinedStateTransitionAndSensitivityMatrixInterface, output_times: list[float]) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]:
    """Function to propagate system covariance through time.
    
    Function to propagate the covariance of a given system through time.
    The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.
    
    
    Parameters
    ----------
    initial_covariance : numpy.ndarray[numpy.float64[m, n]]
        System covariance matrix (symmetric and positive semi-definite) at initial time.
        Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)
    
    state_transition_interface : :class:`~tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface`
        Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.
    
    output_times : List[ float ]
        Times at which the propagated covariance matrix shall be reported.
        Note that this argument has no impact on the integration time-steps of the covariance propagation,
        which always adheres to the integrator settings that the `state_transition_interface` links to.
        Output times which do not coincide with integration time steps are calculated via interpolation.
    
    Returns
    -------
    Dict[ float, numpy.ndarray[numpy.float64[m, n]] ]
        Dictionary reporting the propagated covariances at each output time."""

def propagate_covariance_rsw_split_output(covariance_output: CovarianceAnalysisOutput, estimator: Estimator, output_times: list[float]) -> tuple[list[float], list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]]:
    """No documentation found."""

def propagate_covariance_split_output(initial_covariance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], state_transition_interface: simulator.CombinedStateTransitionAndSensitivityMatrixInterface, output_times: list[float]) -> tuple[list[float], list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]]:
    """No documentation found."""

def propagate_formal_errors(initial_covariance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], state_transition_interface: simulator.CombinedStateTransitionAndSensitivityMatrixInterface, output_times: list[float]) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
    """Function to propagate system formal errors through time.
    
    Function to propagate the formal errors of a given system through time.
    Note that in practice the entire covariance matrix is propagated, but only the formal errors (variances) are reported at the output times.
    The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.
    
    
    Parameters
    ----------
    initial_covariance : numpy.ndarray[numpy.float64[m, n]]
        System covariance matrix (symmetric and positive semi-definite) at initial time.
        Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)
    
    state_transition_interface : :class:`~tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface`
        Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.
    
    output_times : List[ float ]
        Times at which the propagated covariance matrix shall be reported.
        Note that this argument has no impact on the integration time-steps of the covariance propagation,
        which always adheres to the integrator settings that the `state_transition_interface` links to.
        Output times which do not coincide with integration time steps are calculated via interpolation.
    
    Returns
    -------
    Dict[ float, numpy.ndarray[numpy.float64[m, 1]] ]
        Dictionary reporting the propagated formal errors at each output time."""

def propagate_formal_errors_rsw_split_output(covariance_output: CovarianceAnalysisOutput, estimator: Estimator, output_times: list[float]) -> tuple[list[float], list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]]:
    """No documentation found."""

def propagate_formal_errors_split_output(initial_covariance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], state_transition_interface: simulator.CombinedStateTransitionAndSensitivityMatrixInterface, output_times: list[float]) -> tuple[list[float], list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]]:
    """Function to propagate system formal errors through time.
    
    Function to propagate the formal errors of a given system through time.
    Note that in practice the entire covariance matrix is propagated, but only the formal errors (variances) are reported at the output times.
    The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.
    
    
    Parameters
    ----------
    initial_covariance : numpy.ndarray[numpy.float64[m, n]]
        System covariance matrix (symmetric and positive semi-definite) at initial time.
        Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)
    
    state_transition_interface : :class:`~tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface`
        Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.
    
    output_times : List[ float ]
        Times at which the propagated covariance matrix shall be reported.
        Note that this argument has no impact on the integration time-steps of the covariance propagation,
        which always adheres to the integrator settings that the `state_transition_interface` links to.
        Output times which do not coincide with integration time steps are calculated via interpolation.
    
    Returns
    -------
    Dict[ float, numpy.ndarray[numpy.float64[m, 1]] ]
        Dictionary reporting the propagated formal errors at each output time."""
PodInput = EstimationInput
PodOutput = EstimationOutput