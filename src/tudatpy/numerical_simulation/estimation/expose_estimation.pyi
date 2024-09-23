import typing
import numpy
__all__ = ['CombinedStateTransitionAndSensitivityMatrixInterface', 'CovarianceAnalysisInput', 'CovarianceAnalysisOutput', 'EstimatableParameterSet', 'EstimationConvergenceChecker', 'EstimationInput', 'EstimationOutput', 'ObservationCollection', 'ObservationSimulator', 'ObservationSimulator_1', 'ObservationSimulator_2', 'ObservationSimulator_3', 'ObservationSimulator_6', 'ObservationViabilityCalculator', 'PodInput', 'PodOutput', 'SingleObservationSet', 'compute_target_angles_and_range', 'compute_target_angles_and_range_vectors', 'create_pseudo_observations_and_models', 'estimation_convergence_checker', 'propagate_covariance', 'propagate_covariance_rsw_split_output', 'propagate_covariance_split_output', 'propagate_formal_errors', 'propagate_formal_errors_rsw_split_output', 'propagate_formal_errors_split_output', 'set_existing_observations', 'simulate_observations', 'single_observation_set']

class CombinedStateTransitionAndSensitivityMatrixInterface:
    """Class establishing an interface with the simulation's State Transition and Sensitivity Matrices.
	
	Class establishing an interface to the State Transition and Sensitivity Matrices.
	Instances of this class are instantiated automatically upon creation of :class:`~tudatpy.numerical_simulation.Estimator` objects,
	using the simulation information in the observation, propagation and integration settings that the :class:`~tudatpy.numerical_simulation.Estimator` instance is linked to.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def full_state_transition_sensitivity_at_epoch(self, time: float, add_central_body_dependency: bool=True, arc_defining_bodies: list[str]=[]) -> numpy.ndarray:
        """
        	:param time:
        		Time at which full concatenated state transition and sensitivity matrix are to be retrieved.
        	:return:
        		Full concatenated state transition and sensitivity matrix at a given time.
        """

    def state_transition_sensitivity_at_epoch(self, time: float, add_central_body_dependency: bool=True, arc_defining_bodies: list[str]=[]) -> numpy.ndarray:
        """
        Function to get the concatenated state transition and sensitivity matrix at a given time.
        
        	Function to get the concatenated state transition and sensitivity matrix at a given time.
        	Entries corresponding to parameters which are not active at the current arc are omitted.
        
        
        	:param time:
        		Time at which concatenated state transition and sensitivity matrix are to be retrieved.
        	:return:
        		Concatenated state transition and sensitivity matrix at a given time.
        """

    @property
    def full_parameter_size(self) -> int:
        """
        Full amount of parameters w.r.t. which partials have been set up via State Transition and Sensitivity Matrices.
        	
        """

    @property
    def sensitivity_size(self) -> int:
        """
        Number of columns in the sensitivity matrix.
        	
        """

    @property
    def state_transition_size(self) -> int:
        """
        Size of the (square) state transition matrix.
        	
        """

class CovarianceAnalysisInput:
    """Class for defining all specific inputs to a covariance analysis.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, observations_and_times: ObservationCollection, inverse_apriori_covariance: numpy.ndarray=..., consider_covariance: numpy.ndarray=...) -> None:
        """
        Class constructor.
        
        	Constructor through which the user can create instances of this class. Note that the weight are all initiated as 1.0, and the default settings of ``define_covariance_settings`` are used.
        
        
        	:param observations_and_times:
        		Total data structure of observations and associated times/link ends/type/etc.
        	:param inverse_apriori_covariance:
        		A priori covariance matrix (unnormalized) of estimated parameters. This should be either a size 0x0 matrix (no a priori information), or a square matrix with the same size as the number of parameters that are considered
        	:return:
        		Instance of the :class:`~tudatpy.numerical_simulation.estimation.CovarianceAnalysisInput` class, defining the data and other settings to be used for the covariance analysis.
        """

    def define_covariance_settings(self, reintegrate_equations_on_first_iteration: bool=True, reintegrate_variational_equations: bool=True, save_design_matrix: bool=True, print_output_to_terminal: bool=True, limit_condition_number_for_warning: float=100000000.0) -> None:
        """
        Function to define specific settings for covariance analysis process
        
        	Function to define specific settings for covariance analysis process
        
        
        	:param reintegrate_equations:
        		Boolean denoting whether the dynamics and variational equations are to be reintegrated
        		or if existing values are to be used to perform first iteration.
        
        	:param reintegrate_variational_equations:
        		Boolean denoting whether the variational equations are to be reintegrated during estimation
        		(if this is set to False, and ``reintegrate_equations`` to true, only the dynamics are re-integrated)
        
        	:param save_design_matrix:
        		Boolean denoting whether to save the partials matrix (also called design matrix) :math:`\\mathbf{H}` in the output. Setting this to false makes the
        		:math:`\\mathbf{H}` matrix unavailable to the user, with the advantage of lower RAM usage.
        
        	:param print_output_to_terminal:
        		Boolean denoting whether to print covariance-analysis-specific output to the terminal when running the estimation.
        
        	:return:
        		Function modifies the object in-place.
        """

    def set_constant_single_observable_and_link_end_vector_weight(self, observable_type: ..., link_ends: dict[..., ...], weight: numpy.ndarray) -> None:
        ...

    def set_constant_single_observable_and_link_end_weight(self, observable_type: ..., link_ends: dict[..., ...], weight: float) -> None:
        ...

    def set_constant_single_observable_vector_weight(self, observable_type: ..., weight: numpy.ndarray) -> None:
        ...

    def set_constant_single_observable_weight(self, observable_type: ..., weight: float) -> None:
        ...

    def set_constant_vector_weight_per_observable(self, weight_per_observable: dict[..., numpy.ndarray]) -> None:
        ...

    def set_constant_weight(self, weight: float) -> None:
        """
        Function to set a constant weight matrix for all observables.
        
        	Function to set a constant weight matrix for all observables.
        	The weights are applied to all observations managed by the given PodInput object.
        
        
        	:param constant_weight:
        		Constant weight factor that is to be applied to all observations.
        	:return:
        		Function modifies the object in-place.
        """

    def set_constant_weight_per_observable(self, weight_per_observable: dict[..., float]) -> None:
        """
        Function to set a constant weight matrix for a given type of observable.
        
        	Function to set a constant weight matrix for a given type of observable.
        	The weights are applied to all observations of the observable type specified by the `weight_per_observable` parameter.
        
        
        	:param constant_weight:
        		Constant weight factor that is to be applied to all observations.
        	:return:
        		Function modifies the object in-place.
        """

    def set_total_single_observable_and_link_end_vector_weight(self, observable_type: ..., link_ends: dict[..., ...], weight_vector: numpy.ndarray) -> None:
        ...

    @property
    def weight_matrix_diagonal(self) -> numpy.ndarray:
        """
        Complete diagonal of the weights matrix that is to be used
        	
        """

    @weight_matrix_diagonal.setter
    def weight_matrix_diagonal(self, arg1: numpy.ndarray) -> None:
        ...

class CovarianceAnalysisOutput:
    """Class collecting all outputs from the covariance analysis process.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def consider_covariance_contribution(self) -> numpy.ndarray:
        ...

    @property
    def consider_normalization_factors(self) -> numpy.ndarray:
        ...

    @property
    def correlations(self) -> numpy.ndarray:
        """
        Correlation matrix of the estimation result. Entry :math:`i,j` is equal to :math:`P_{i,j}/(\\sigma_{i}\\sigma_{j})`
        	
        """

    @property
    def covariance(self) -> numpy.ndarray:
        """
        (Unnormalized) estimation covariance matrix :math:`\\mathbf{P}`.
        	
        """

    @property
    def design_matrix(self) -> numpy.ndarray:
        """
        Matrix of unnormalized partial derivatives :math:`\\mathbf{H}=\x0crac{\\partial\\mathbf{h}}{\\partial\\mathbf{p}}`.
        	
        """

    @property
    def formal_errors(self) -> numpy.ndarray:
        """
        Formal error vector :math:`\x08oldsymbol{\\sigma}` of the estimation result (e.g. square root of diagonal entries of covariance)s
        	
        """

    @property
    def inverse_covariance(self) -> numpy.ndarray:
        """
        (Unnormalized) inverse estimation covariance matrix :math:`\\mathbf{P}^{-1}`.
        	
        """

    @property
    def inverse_normalized_covariance(self) -> numpy.ndarray:
        """
        Normalized inverse estimation covariance matrix :math:`\\mathbf{	ilde{P}}^{-1}`.
        	
        """

    @property
    def normalization_terms(self) -> numpy.ndarray:
        """
        Vector of normalization terms used for covariance and design matrix
        	
        """

    @property
    def normalized_covariance(self) -> numpy.ndarray:
        """
        Normalized estimation covariance matrix :math:`\\mathbf{	ilde{P}}`.
        	
        """

    @property
    def normalized_covariance_with_consider_parameters(self) -> numpy.ndarray:
        ...

    @property
    def normalized_design_matrix(self) -> numpy.ndarray:
        """
        Matrix of normalized partial derivatives :math:`	ilde{\\mathbf{H}}`.
        	
        """

    @property
    def normalized_design_matrix_consider_parameters(self) -> numpy.ndarray:
        ...

    @property
    def unnormalized_covariance_with_consider_parameters(self) -> numpy.ndarray:
        ...

    @property
    def weighted_design_matrix(self) -> numpy.ndarray:
        """
        Matrix of weighted partial derivatives, equal to :math:`\\mathbf{W}^{1/2}{\\mathbf{H}}`
        	
        """

    @property
    def weighted_normalized_design_matrix(self) -> numpy.ndarray:
        """
        Matrix of weighted, normalized partial derivatives, equal to :math:`\\mathbf{W}^{1/2}	ilde{\\mathbf{H}}`
        	
        """

class EstimatableParameterSet:
    """Class containing a consolidated set of estimatable parameters.
	
	Class containing a consolidated set of estimatable parameters, linked to the environment and acceleration settings of the simulation.
	The user typically creates instances of this class via the :func:`~tudatpy.numerical_simulation.estimation_setup.create_parameters_to_estimate` factory function.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def indices_for_parameter_type(self, parameter_type: tuple[..., tuple[str, str]]) -> list[tuple[int, int]]:
        """
        Function to retrieve the indices of a given type of parameter.
        
        	Function to retrieve the index of all parameters of a given type from the parameter set.
        	This function can be very useful, since the order of parameters within the parameter set does not necessarily correspond to the order in which the elements were added to the set.
        
        
        	:param parameter_type:
        		help
        	:return:
        		help
        """

    @property
    def constraints_size(self) -> int:
        """
        Total size of linear constraint that is to be applied during estimation.
        	
        """

    @property
    def initial_multi_arc_states_size(self) -> int:
        """
        Amount of initial state parameters in the set, which are treated in a multi-arc fashion.
        	
        """

    @property
    def initial_single_arc_states_size(self) -> int:
        """
        Amount of initial state parameters in the set, which are treated in a single-arc fashion.
        	
        """

    @property
    def initial_states_size(self) -> int:
        """
        Amount of initial state parameters contained in the set.
        	
        """

    @property
    def parameter_set_size(self) -> int:
        """
        Size of the parameter set, i.e. amount of estimatable parameters contained in the set.
        	
        """

    @property
    def parameter_vector(self) -> numpy.ndarray:
        """
        Vector containing the parameter values of all parameters in the set.
        	
        """

    @parameter_vector.setter
    def parameter_vector(self, arg1: numpy.ndarray) -> None:
        ...

class EstimationConvergenceChecker:
    """Class defining the convergence criteria for an estimation.
	
	Class defining the convergence criteria for an estimation.
	The user typically creates instances of this class via the :func:`~tudatpy.numerical_simulation.estimation.estimation_convergence_checker` factory function.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class EstimationInput(CovarianceAnalysisInput):
    """Class for defining all inputs to the estimation.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, observations_and_times: ObservationCollection, inverse_apriori_covariance: numpy.ndarray=..., convergence_checker: EstimationConvergenceChecker=..., consider_covariance: numpy.ndarray=..., consider_parameters_deviations: numpy.ndarray=..., apply_final_parameter_correction: bool=True) -> None:
        """
        Class constructor.
        
        	Constructor through which the user can create instances of this class.
        
        
        	:param observations_and_times:
        		Total data structure of observations and associated times/link ends/type/etc.
        	:param inverse_apriori_covariance:
        		A priori covariance matrix (unnormalized) of estimated parameters. This should be either a size 0x0 matrix (no a priori information), or a square matrix with the same size as the number of parameters that are considered
        	:param convergence_checker:
        		Object defining when the estimation is converged.
        	:return:
        		Instance of the :class:`~tudatpy.numerical_simulation.estimation.EstimationInput` class, defining the data and other settings to be used for the estimation.
        """

    def define_estimation_settings(self, reintegrate_equations_on_first_iteration: bool=True, reintegrate_variational_equations: bool=True, save_design_matrix: bool=True, print_output_to_terminal: bool=True, save_residuals_and_parameters_per_iteration: bool=True, save_state_history_per_iteration: bool=False, limit_condition_number_for_warning: float=100000000.0, condition_number_warning_each_iteration: bool=True) -> None:
        """
        Function to define specific settings for the estimation process
        
        	Function to define specific settings for covariance analysis process
        
        
        	:param reintegrate_equations_on_first_iteration:
        		Boolean denoting whether the dynamics and variational equations are to be reintegrated
        		or if existing values are to be used to perform first iteration.
        
        	:param reintegrate_variational_equations:
        		Boolean denoting whether the variational equations are to be reintegrated during estimation
        		(if this is set to False, and ``reintegrate_equations_on_first_iteration`` to true, only the dynamics are re-integrated)
        
        	:param save_design_matrix:
        		Boolean denoting whether to save the partials matrix (also called design matrix) :math:`\\mathbf{H}` in the output. Setting this to false makes the
        		:math:`\\mathbf{H}` matrix unavailable to the user, with the advantage of lower RAM usage.
        
        	:param print_output_to_terminal:
        		Boolean denoting whether to print covariance-analysis-specific output to the terminal when running the estimation.
        
        	:param save_residuals_and_parameters_per_iteration:
        		Boolean denoting whether the residuals and parameters from the each iteration are to be saved.
        
        	:param save_state_history_per_iteration:
        		Boolean denoting whether the state history and dependent variables are to be saved on each iteration.
        
        	:return:
        		Function modifies the object in-place.
        """

class EstimationOutput(CovarianceAnalysisOutput):
    """Class collecting all outputs from the iterative estimation process.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def best_iteration(self) -> int:
        ...

    @property
    def final_parameters(self) -> numpy.ndarray:
        ...

    @property
    def final_residuals(self) -> numpy.ndarray:
        """
        Vector of post-fit observation residuals.
        	
        """

    @property
    def parameter_history(self) -> numpy.ndarray:
        """
        Parameter vectors, concatenated per iteration into a matrix. Column 0 contains pre-estimation values. The :math:`(i+1)^{th}` column has the residuals from the :math:`i^{th}` iteration.
        	
        """

    @property
    def residual_history(self) -> numpy.ndarray:
        """
        Residual vectors, concatenated per iteration into a matrix; the :math:`i^{th}` column has the residuals from the :math:`i^{th}` iteration.
        	
        """

    @property
    def simulation_results_per_iteration(self) -> list[..., ...]:
        """
        List of complete numerical propagation results, with the :math:`i^{th}` entry of thee list thee results of the :math:`i^{th}` propagation
        	
        """

class ObservationCollection:
    """Class collecting all observations and associated data for use in an estimation.
	
	Class containing the full set of observations and associated data, typically for input into the estimation. When using simulated data,
	this class is instantiated via a call to the :func:`~tudatpy.numerical_simulation.estimation.simulate_observations` function. More information is provided
	on the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, observation_sets: list[..., double, ...]) -> None:
        ...

    def get_link_definitions_for_observables(self, observable_type: ...) -> list[...]:
        ...

    def get_single_link_and_type_observations(self, observable_type: ..., link_definition: ...) -> list[..., double, ...]:
        """
        Function to get all observation sets for a given observable type and link definition.
        
        	:param observable_type:
        		Observable type of which observations are to be simulated.
        	:param link_ends:
        		Link ends for which observations are to be simulated.
        	:return:
        		List of observation sets for given observable type and link definition.
        """

    @property
    def concatenated_link_definition_ids(self) -> list[int]:
        """
        Vector containing concatenated indices identifying the link ends. Each set of link ends is assigned a unique integer identifier (for a given instance of this class). The definition of a given integer identifier with the link ends is given by this class' :func:`link_definition_ids` function. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order of the present vector.
        	
        """

    @property
    def concatenated_observations(self) -> numpy.ndarray:
        """
        Vector containing concatenated observable values. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order
        	
        """

    @property
    def concatenated_times(self) -> list[float]:
        """
        Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order
        	
        """

    @property
    def link_definition_ids(self) -> dict[int, dict[..., ...]]:
        """
        Dictionaty mapping a link end integer identifier to the specific link ends
        	
        """

    @property
    def link_definitions_per_observable(self) -> dict[..., list[...]]:
        ...

    @property
    def link_ends_per_observable_type(self) -> dict[..., list[dict[..., ...]]]:
        ...

    @property
    def observable_type_start_index_and_size(self) -> dict[..., tuple[int, int]]:
        """
        Dictionary defining per obervable type (dict key), the index in the full observation vector (:func:`concatenated_observations`) where the given observable type starts, and the number of subsequent entries in this vector containing a value of an observable of this type
        	
        """

    @property
    def observation_set_start_index_and_size(self) -> dict[..., dict[int, list[tuple[int, int]]]]:
        """
        The nested dictionary/list returned by this property mirrors the structure of the :func:`sorted_observation_sets` property of this class. The present function provides the start index and size of the observables in the full observation vector that come from the correspoding `SingleObservationSet` in the :func:`sorted_observation_sets` Consequently, the present property returns a nested dictionary defining per obervable type, link end identifier, and `SingleObservationSet` index (for the given observable type and link end identifier), where the observables in the given `SingleObservationSet` starts, and the number of subsequent entries in this vector containing data from it.
        	
        """

    @property
    def observation_vector_size(self) -> int:
        """
        Length of the total vector of observations
        	
        """

    @property
    def sorted_observation_sets(self) -> dict[..., dict[int, list[..., double, ...]]]:
        """
        The nested dictionary/list contains the list of `SingleObservationSet` objects, in the same method as they are stored internally in the present class. Specifics on the storage order are given in the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_
        	
        """

class ObservationSimulator:
    """Class hosting the functionality for simulating observations.
	
	Class hosting the functionality for simulating a given observable over a defined link geometry.
	Instances of this class are automatically created from the given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` objects upon instantiation of the :class:`~tudatpy.numerical_simulation.Estimator` class.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class ObservationSimulator_1(ObservationSimulator):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class ObservationSimulator_2(ObservationSimulator):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class ObservationSimulator_3(ObservationSimulator):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class ObservationSimulator_6(ObservationSimulator):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class ObservationViabilityCalculator:
    """Template class for observation viability calculators.
	
	Template class for classes which conducts viability calculations on simulated observations.
	Instances of the applicable ObservationViabilityCalculators are automatically created from the given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects during the simulation of observations (:func:`~tudatpy.numerical_simulation.estimation.simulate_observations`).
	The user typically does not interact directly with this class.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def is_observation_viable(self, link_end_states: list[numpy.ndarray], link_end_times: list[float]) -> bool:
        """
        Function to check whether an observation is viable.
        
        	Function to check whether an observation is viable.
        	The calculation is performed based on the given times and link end states.
        	Note, that this function is called automatically during the simulation of observations.
        	Direct calls to this function are generally not required.
        
        
        	:param link_end_states:
        		Vector of states of the link ends involved in the observation.
        	:param link_end_times:
        		Vector of times at the link ends involved in the observation.
        	:return:
        		True if observation is viable, false if not.
        """

class SingleObservationSet:
    """Class collecting a single set of observations and associated data, of a given observable type, link ends, and ancilliary data.
	"""
    weights_vector: numpy.ndarray

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def ancilliary_settings(self) -> ...:
        """
        Ancilliary settings all stored observations
        	
        """

    @property
    def concatenated_observations(self) -> numpy.ndarray:
        """
        Concatenated vector of all stored observations
        	
        """

    @property
    def link_definition(self) -> ...:
        """
        Definition of the link ends for which the object stores observations
        	
        """

    @property
    def list_of_observations(self) -> list[numpy.ndarray]:
        """
        List of separate stored observations. Each entry of this list is a vector containing a single observation. In cases where the observation is single-valued (range, Doppler), the vector is size 1, but for multi-valued observations such as angular position, each vector in the list will have size >1
        	
        """

    @property
    def observable_type(self) -> ...:
        """
        Type of observable for which the object stores observations
        	
        """

    @property
    def observation_times(self) -> list[float]:
        """
        Reference time for each of the observations in ``list_of_observations``
        	
        """

    @property
    def observations_history(self) -> dict[float, numpy.ndarray]:
        """
        Dictionary of observations sorted by time. Created by making a dictionaty with ``observation_times`` as keys and  ``list_of_observations`` as values
        	
        """

    @property
    def reference_link_end(self) -> ...:
        """
        Reference link end for all stored observations
        	
        """

def compute_target_angles_and_range(bodies: ..., station_id: tuple[str, str], target_body: str, observation_times: list[float], is_station_transmitting: bool) -> dict[float, numpy.ndarray]:
    """Function to compute the azimuth angle, elevation angle and range at a ground station.
	
	Function to compute the azimuth angle, elevation angle and range at a ground station. This functions is provided as a function of
	convenience, to prevent users having to manually define the relevant settings for this often-needed functionality. This function
	takes an observing station and a target body as input, and provides the observed angles and current range (without correction for aberrations, with correction for light time)
	as observed at that station
	
	
	:param bodies:
			System of bodies that defines the full physical environment
	
	:param station_id:
			Identifier for the observing station, as a pair of strings: the body name and the station name.
	
	:param target_body:
			Name of body which is observed by ground station
	
	:param observation_times:
			List of times at which the ground station observations are to be analyzed
	
	:param is_station_transmitting:
			Boolean defining whether the observation times define times at which the station is transmitting to, or receiving from, the ground station.
			This has an impact on the whether the light-time is computed forward or backward in time from the ground station to the target
	
	:return:
			Dictionary with the required output. Key defines the observation time, the value is an array of size three containing entry 0 - elevation angle, entry 1 - azimuth angle, entry 2 - range
	"""

def compute_target_angles_and_range_vectors(bodies: ..., station_id: tuple[str, str], target_body: str, observation_times: list[float], is_station_transmitting: bool) -> tuple[list[float], list[numpy.ndarray]]:
    """Function to compute the azimuth angle, elevation angle and range at a ground station.
	
	Function to compute the azimuth angle, elevation angle and range at a ground station. This functions is provided as a function of
	convenience, to prevent users having to manually define the relevant settings for this often-needed functionality. This function
	takes an observing station and a target body as input, and provides the observed angles and current range (without correction for aberrations, with correction for light time)
	as observed at that station
	
	
	:param bodies:
			System of bodies that defines the full physical environment
	
	:param station_id:
			Identifier for the observing station, as a pair of strings: the body name and the station name.
	
	:param target_body:
			Name of body which is observed by ground station
	
	:param observation_times:
			List of times at which the ground station observations are to be analyzed
	
	:param is_station_transmitting:
			Boolean defining whether the observation times define times at which the station is transmitting to, or receiving from, the ground station.
			This has an impact on the whether the light-time is computed forward or backward in time from the ground station to the target
	
	:return:
			Dictionary with the required output. Key defines the observation time, the value is an array of size three containing entry 0 - elevation angle, entry 1 - azimuth angle, entry 2 - range
	"""

def create_pseudo_observations_and_models(bodies: ..., observed_bodies: list[str], central_bodies: list[str], initial_time: float, final_time: float, time_step: float) -> tuple[list[...], ..., double, ...]:
    ...

def estimation_convergence_checker(maximum_iterations: int=5, minimum_residual_change: float=0.0, minimum_residual: float=0.0, number_of_iterations_without_improvement: int=2) -> EstimationConvergenceChecker:
    """Function for creating an :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` object.
	
	Function for creating an :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` object, which is required for defining the convergence criteria of an estimation.
	
	
	:param maximum_iterations:
			Maximum number of allowed iterations for estimation.
	:param minimum_residual_change:
			Minimum required change in residual between two iterations.
	:param minimum_residual:
			Minimum value of observation residual below which estimation is converged.
	:param number_of_iterations_without_improvement:
			Number of iterations without reduction of residual.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` class, defining the convergence criteria for an estimation.
	"""

def propagate_covariance(initial_covariance: numpy.ndarray, state_transition_interface: CombinedStateTransitionAndSensitivityMatrixInterface, output_times: list[float]) -> dict[float, numpy.ndarray]:
    """Function to propagate system covariance through time.
	
	Function to propagate the covariance of a given system through time.
	The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.
	
	
	:param initial_covariance:
			System covariance matrix (symmetric and positive semi-definite) at initial time.
			Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)
	
	:param state_transition_interface:
			Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.
	
	:param output_times:
			Times at which the propagated covariance matrix shall be reported.
			Note that this argument has no impact on the integration time-steps of the covariance propagation,
			which always adheres to the integrator settings that the `state_transition_interface` links to.
			Output times which do not coincide with integration time steps are calculated via interpolation.
	
	:return:
			Dictionary reporting the propagated covariances at each output time.
	"""

def propagate_covariance_rsw_split_output(*args, **kwargs) -> tuple[list[float], list[numpy.ndarray]]:
    ...

def propagate_covariance_split_output(initial_covariance: numpy.ndarray, state_transition_interface: CombinedStateTransitionAndSensitivityMatrixInterface, output_times: list[float]) -> tuple[list[float], list[numpy.ndarray]]:
    ...

def propagate_formal_errors(initial_covariance: numpy.ndarray, state_transition_interface: CombinedStateTransitionAndSensitivityMatrixInterface, output_times: list[float]) -> dict[float, numpy.ndarray]:
    """Function to propagate system formal errors through time.
	
	Function to propagate the formal errors of a given system through time.
	Note that in practice the entire covariance matrix is propagated, but only the formal errors (variances) are reported at the output times.
	The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.
	
	
	:param initial_covariance:
			System covariance matrix (symmetric and positive semi-definite) at initial time.
			Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)
	
	:param state_transition_interface:
			Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.
	
	:param output_times:
			Times at which the propagated covariance matrix shall be reported.
			Note that this argument has no impact on the integration time-steps of the covariance propagation,
			which always adheres to the integrator settings that the `state_transition_interface` links to.
			Output times which do not coincide with integration time steps are calculated via interpolation.
	
	:return:
			Dictionary reporting the propagated formal errors at each output time.
	"""

def propagate_formal_errors_rsw_split_output(*args, **kwargs) -> tuple[list[float], list[numpy.ndarray]]:
    ...

def propagate_formal_errors_split_output(initial_covariance: numpy.ndarray, state_transition_interface: CombinedStateTransitionAndSensitivityMatrixInterface, output_times: list[float]) -> tuple[list[float], list[numpy.ndarray]]:
    """Function to propagate system formal errors through time.
	
	Function to propagate the formal errors of a given system through time.
	Note that in practice the entire covariance matrix is propagated, but only the formal errors (variances) are reported at the output times.
	The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.
	
	
	:param initial_covariance:
			System covariance matrix (symmetric and positive semi-definite) at initial time.
			Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)
	
	:param state_transition_interface:
			Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.
	
	:param output_times:
			Times at which the propagated covariance matrix shall be reported.
			Note that this argument has no impact on the integration time-steps of the covariance propagation,
			which always adheres to the integrator settings that the `state_transition_interface` links to.
			Output times which do not coincide with integration time steps are calculated via interpolation.
	
	:return:
			Dictionary reporting the propagated formal errors at each output time.
	"""

def set_existing_observations(observations: dict[..., tuple[dict[..., ...], tuple[list[numpy.ndarray], list[float]]]], reference_link_end: ..., ancilliary_settings_per_observatble: dict[..., ...]={}) -> ...:
    ...

def simulate_observations(simulation_settings: list[...], observation_simulators: list[ObservationSimulator], bodies: ...) -> ...:
    """Function to simulate observations.
	
	Function to simulate observations from set observation simulators and observation simulator settings.
	Automatically iterates over all provided observation simulators, generating the full set of simulated observations.
	
	
	:param observation_to_simulate:
			List of settings objects, each object providing the observation time settings for simulating one type of observable and link end set.
	
	:param observation_simulators:
			List of :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` objects, each object hosting the functionality for simulating one type of observable and link end set.
	
	:param bodies:
			Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.
	
	:return:
			Object collecting all products of the observation simulation.
	"""

def single_observation_set(observable_type: ..., link_definition: ..., observations: list[numpy.ndarray], observation_times: list[float], reference_link_end: ..., ancilliary_settings: ...=None) -> SingleObservationSet:
    ...
PodInput = EstimationInput
PodOutput = EstimationOutput