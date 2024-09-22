import numpy
import typing
__all__ = ['CustomAccelerationPartialSettings', 'EmpiricalAccelerationComponents', 'EmpiricalAccelerationFunctionalShapes', 'EstimatableParameterSettings', 'EstimatableParameterTypes', 'absolute_observation_bias', 'across_track_empirical_acceleration_component', 'along_track_empirical_acceleration_component', 'arc_wise_constant_drag_coefficient_type', 'arc_wise_empirical_acceleration_coefficients_type', 'arc_wise_initial_body_state_type', 'arc_wise_radiation_pressure_coefficient_type', 'arc_wise_time_drift_observation_bias_type', 'arcwise_absolute_observation_bias', 'arcwise_constant_additive_observation_bias_type', 'arcwise_constant_drag_coefficient', 'arcwise_constant_empirical_acceleration_terms', 'arcwise_constant_relative_observation_bias_type', 'arcwise_empirical_accelerations', 'arcwise_radiation_pressure_coefficient', 'arcwise_relative_observation_bias', 'arcwise_time_drift_observation_bias', 'constant_additive_observation_bias_type', 'constant_drag_coefficient', 'constant_drag_coefficient_type', 'constant_empirical', 'constant_empirical_acceleration_terms', 'constant_relative_observation_bias_type', 'constant_rotation_rate', 'constant_rotation_rate_type', 'constant_time_drift_observation_bias_type', 'core_factor', 'core_factor_type', 'cosine_empirical', 'custom_analytical_partial', 'custom_numerical_partial', 'custom_parameter', 'desaturation_delta_v_values_type', 'direct_dissipation_tidal_time_lag_type', 'direct_tidal_dissipation_time_lag', 'empirical_acceleration_coefficients_type', 'empirical_accelerations', 'equivalence_principle_lpi_violation_parameter_type', 'free_core_nutation_rate', 'free_core_nutation_rate_type', 'full_degree_tidal_love_number_type', 'full_empirical_acceleration_terms', 'gravitational_parameter', 'gravitational_parameter_type', 'ground_station_position', 'ground_station_position_type', 'initial_body_state_type', 'initial_rotational_body_state_type', 'initial_states', 'inverse_tidal_quality_factor', 'inverse_tidal_quality_factor_type', 'mean_moment_of_inertia', 'mean_moment_of_inertia_type', 'monomial_full_block_gravity_field_variation_amplitudes', 'monomial_gravity_field_variation_amplitudes', 'order_invariant_k_love_number', 'order_varying_k_love_number', 'periodic_spin_variation_type', 'periodic_spin_variations', 'polar_motion_amplitude_type', 'polar_motion_amplitudes', 'polynomial_gravity_field_variation_amplitudes', 'ppn_parameter_beta', 'ppn_parameter_beta_type', 'ppn_parameter_gamma', 'ppn_parameter_gamma_type', 'quasi_impulsive_shots', 'radial_empirical_acceleration_component', 'radiation_pressure_coefficient', 'radiation_pressure_coefficient_type', 'radiation_pressure_target_direction_scaling', 'radiation_pressure_target_perpendicular_direction_scaling', 'relative_observation_bias', 'rotation_pole_position', 'rotation_pole_position_type', 'scaled_longitude_libration_amplitude', 'sine_empirical', 'single_degree_variable_tidal_love_number_type', 'spherical_harmonics_c_coefficients', 'spherical_harmonics_c_coefficients_block', 'spherical_harmonics_cosine_coefficient_block_type', 'spherical_harmonics_s_coefficients', 'spherical_harmonics_s_coefficients_block', 'spherical_harmonics_sine_coefficient_block_type', 'time_drift_observation_bias', 'yarkovsky_parameter']

class CustomAccelerationPartialSettings:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class EmpiricalAccelerationComponents:
    """Members:
	
	radial_empirical_acceleration_component
	
	along_track_empirical_acceleration_component
	
	across_track_empirical_acceleration_component
	"""
    __members__: typing.ClassVar[dict[str, EmpiricalAccelerationComponents]]
    across_track_empirical_acceleration_component: typing.ClassVar[EmpiricalAccelerationComponents]
    along_track_empirical_acceleration_component: typing.ClassVar[EmpiricalAccelerationComponents]
    radial_empirical_acceleration_component: typing.ClassVar[EmpiricalAccelerationComponents]

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

class EmpiricalAccelerationFunctionalShapes:
    """Members:
	
	constant_empirical
	
	sine_empirical
	
	cosine_empirical
	"""
    __members__: typing.ClassVar[dict[str, EmpiricalAccelerationFunctionalShapes]]
    constant_empirical: typing.ClassVar[EmpiricalAccelerationFunctionalShapes]
    cosine_empirical: typing.ClassVar[EmpiricalAccelerationFunctionalShapes]
    sine_empirical: typing.ClassVar[EmpiricalAccelerationFunctionalShapes]

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

class EstimatableParameterSettings:
    """Base class to defining settings of parameter to be estimated.
	
	Functional (base) class for settings of model parameter to be estimated.
	Settings of simple parameters types are managed via this class, more complex parameter types are handled by specialised derivates of this class.
	Instances of either base or derived class can be created via dedicated factory functions.
	"""
    custom_partial_settings: list[CustomAccelerationPartialSettings]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class EstimatableParameterTypes:
    """Enumeration of model parameters that are available for estimation.
	In order to establish a parameter estimation settings for a parameter of a certain type, use the factory function dedicated to this parameter type.
	Note that not all of the listed types might be accessible via factory functions in the python interface yet.
	
	
	
		:member arc_wise_initial_body_state_type:
		:member initial_body_state_type:
		:member initial_rotational_body_state_type:
		:member constant_drag_coefficient_type:
		:member arc_wise_constant_drag_coefficient_type:
		:member radiation_pressure_coefficient_type:
		:member arc_wise_radiation_pressure_coefficient_type:
		:member empirical_acceleration_coefficients_type:
		:member arc_wise_empirical_acceleration_coefficients_type:
		:member desaturation_delta_v_values_type:
		:member gravitational_parameter_type:
		:member spherical_harmonics_cosine_coefficient_block_type:
		:member spherical_harmonics_sine_coefficient_block_type:
		:member mean_moment_of_inertia_type:
		:member constant_rotation_rate_type:
		:member rotation_pole_position_type:
		:member polar_motion_amplitude_type:
		:member core_factor_type:
		:member free_core_nutation_rate_type:
		:member periodic_spin_variation_type:
		:member constant_additive_observation_bias_type:
		:member arcwise_constant_additive_observation_bias_type:
		:member constant_relative_observation_bias_type:
		:member arcwise_constant_relative_observation_bias_type:
		:member ground_station_position_type:
		:member full_degree_tidal_love_number_type:
		:member single_degree_variable_tidal_love_number_type:
		:member direct_dissipation_tidal_time_lag_type:
		:member ppn_parameter_gamma_type:
		:member ppn_parameter_beta_type:
		:member equivalence_principle_lpi_violation_parameter_type:
	
	
	Members:
	
	arc_wise_initial_body_state_type
	
	initial_body_state_type
	
	initial_rotational_body_state_type
	
	gravitational_parameter_type
	
	constant_drag_coefficient_type
	
	radiation_pressure_coefficient_type
	
	arc_wise_radiation_pressure_coefficient_type
	
	spherical_harmonics_cosine_coefficient_block_type
	
	spherical_harmonics_sine_coefficient_block_type
	
	constant_rotation_rate_type
	
	rotation_pole_position_type
	
	constant_additive_observation_bias_type
	
	arcwise_constant_additive_observation_bias_type
	
	constant_relative_observation_bias_type
	
	arcwise_constant_relative_observation_bias_type
	
	ppn_parameter_gamma_type
	
	ppn_parameter_beta_type
	
	ground_station_position_type
	
	equivalence_principle_lpi_violation_parameter_type
	
	empirical_acceleration_coefficients_type
	
	arc_wise_empirical_acceleration_coefficients_type
	
	full_degree_tidal_love_number_type
	
	single_degree_variable_tidal_love_number_type
	
	direct_dissipation_tidal_time_lag_type
	
	mean_moment_of_inertia_type
	
	arc_wise_constant_drag_coefficient_type
	
	periodic_spin_variation_type
	
	polar_motion_amplitude_type
	
	core_factor_type
	
	free_core_nutation_rate_type
	
	desaturation_delta_v_values_type
	
	constant_time_drift_observation_bias_type
	
	arc_wise_time_drift_observation_bias_type
	
	inverse_tidal_quality_factor_type
	"""
    __members__: typing.ClassVar[dict[str, EstimatableParameterTypes]]
    arc_wise_constant_drag_coefficient_type: typing.ClassVar[EstimatableParameterTypes]
    arc_wise_empirical_acceleration_coefficients_type: typing.ClassVar[EstimatableParameterTypes]
    arc_wise_initial_body_state_type: typing.ClassVar[EstimatableParameterTypes]
    arc_wise_radiation_pressure_coefficient_type: typing.ClassVar[EstimatableParameterTypes]
    arc_wise_time_drift_observation_bias_type: typing.ClassVar[EstimatableParameterTypes]
    arcwise_constant_additive_observation_bias_type: typing.ClassVar[EstimatableParameterTypes]
    arcwise_constant_relative_observation_bias_type: typing.ClassVar[EstimatableParameterTypes]
    constant_additive_observation_bias_type: typing.ClassVar[EstimatableParameterTypes]
    constant_drag_coefficient_type: typing.ClassVar[EstimatableParameterTypes]
    constant_relative_observation_bias_type: typing.ClassVar[EstimatableParameterTypes]
    constant_rotation_rate_type: typing.ClassVar[EstimatableParameterTypes]
    constant_time_drift_observation_bias_type: typing.ClassVar[EstimatableParameterTypes]
    core_factor_type: typing.ClassVar[EstimatableParameterTypes]
    desaturation_delta_v_values_type: typing.ClassVar[EstimatableParameterTypes]
    direct_dissipation_tidal_time_lag_type: typing.ClassVar[EstimatableParameterTypes]
    empirical_acceleration_coefficients_type: typing.ClassVar[EstimatableParameterTypes]
    equivalence_principle_lpi_violation_parameter_type: typing.ClassVar[EstimatableParameterTypes]
    free_core_nutation_rate_type: typing.ClassVar[EstimatableParameterTypes]
    full_degree_tidal_love_number_type: typing.ClassVar[EstimatableParameterTypes]
    gravitational_parameter_type: typing.ClassVar[EstimatableParameterTypes]
    ground_station_position_type: typing.ClassVar[EstimatableParameterTypes]
    initial_body_state_type: typing.ClassVar[EstimatableParameterTypes]
    initial_rotational_body_state_type: typing.ClassVar[EstimatableParameterTypes]
    inverse_tidal_quality_factor_type: typing.ClassVar[EstimatableParameterTypes]
    mean_moment_of_inertia_type: typing.ClassVar[EstimatableParameterTypes]
    periodic_spin_variation_type: typing.ClassVar[EstimatableParameterTypes]
    polar_motion_amplitude_type: typing.ClassVar[EstimatableParameterTypes]
    ppn_parameter_beta_type: typing.ClassVar[EstimatableParameterTypes]
    ppn_parameter_gamma_type: typing.ClassVar[EstimatableParameterTypes]
    radiation_pressure_coefficient_type: typing.ClassVar[EstimatableParameterTypes]
    rotation_pole_position_type: typing.ClassVar[EstimatableParameterTypes]
    single_degree_variable_tidal_love_number_type: typing.ClassVar[EstimatableParameterTypes]
    spherical_harmonics_cosine_coefficient_block_type: typing.ClassVar[EstimatableParameterTypes]
    spherical_harmonics_sine_coefficient_block_type: typing.ClassVar[EstimatableParameterTypes]

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

def absolute_observation_bias(link_ends: ..., observable_type: ...) -> EstimatableParameterSettings:
    """Function for defining parameter settings for an absolute observation bias.
	
	Factory function for creating a (linear sensitivity) parameter settings object for an observation's absolute bias parameter.
	Using the absolute observation bias as estimatable parameter requires:
	  * The observation model (corresponding to the `link_ends` and `observable_type`) to include an absolute bias (:func:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias`)
	
	
	:param link_ends:
			Set of link ends that define the geometry of the biased observations.
	
	:param observable_type:
			Observable type of the biased observations.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ConstantObservationBiasEstimatableParameterSettings`
			for the specified observation's arc-wise absolute bias.
	"""

def arcwise_absolute_observation_bias(link_ends: ..., observable_type: ..., arc_start_times: list[float], time_link_end: ...) -> EstimatableParameterSettings:
    """Function for defining parameter settings for arc-wise absolute observation bias.
	
	Factory function for creating a (linear sensitivity) parameter settings object for the arc-wise treatment of an observation's absolute bias parameter.
	Using the arc-wise absolute observation bias as estimatable parameter requires
	  * The observation model (corresponding to the `link_ends` and `observable_type`) to include an arc-wise absolute bias (:func:`~tudatpy.numerical_simulation.estimation_setup.observation.arcwise_absolute_bias`)
	
	
	:param link_ends:
			Set of link ends that define the geometry of the biased observations.
	
	:param observable_type:
			Observable type of the biased observations.
	:param arc_start_times:
			List of times at which the arcs over which the bias is to be estimated will start.
	:param time_link_end:
			The link end type (transmitter, receiver, etc.) at which the arc_start_times is evaluated.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ArcWiseConstantObservationBiasEstimatableParameterSettings`
			for the specified observation's arc-wise absolute bias.
	"""

def arcwise_constant_drag_coefficient(body: str, arc_initial_times: list[float]) -> EstimatableParameterSettings:
    """Function for defining parameter settings for arc-wise constant drag coefficients.
	
	Factory function for creating (linear sensitivity) parameter settings object for arc-wise constant drag coefficients (arc-wise version of :func:``~tudatpy.numerical_simulation.estimation_setup.parameter.constant_drag_coefficient`).
	Using the arc-wise constant drag coefficient as an estimatable parameter requires:
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.constant` aerodynamic interface to be defined for the body specified by the ``body`` parameter
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.aerodynamic` acceleration
	
	Note: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the drag coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.
	
	
	:param body:
			Name of the body, with whose drag acceleration model the estimatable parameter is associated.
	
	:param arc_initial_times:
			List of times at which the arcs over which the drag coefficient is to be estimated will start.
	
	:return:
			Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ArcWiseDragCoefficientEstimatableParameterSettings` class
			for arc-wise treatment of the specified body's constant drag coefficient.
	"""

def arcwise_constant_empirical_acceleration_terms(body: str, centralBody: str, arc_start_times: list[float]) -> EstimatableParameterSettings:
    """Function for defining parameter settings for arc-wise constant empirical acceleration terms.
	
	As :func:`~tudatpy.numerical_simulation.estimation_setup.parameter.arcwise_empirical_accelerations`, but only using the constant R, S and W components (no sine or cosine term estimation). This function is added as a function of convenience
	
	
	:param body:
			Name of the body, with whose empirical acceleration model the estimatable parameter is associated.
	
	:param centralBody:
			Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration
	
	:param arc_initial_times:
			List of times at which the arcs over which the empirical accelerations are to be estimated will start.
	
	:return:
			Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EmpiricalAccelerationEstimatableParameterSettings` class
			for the specified body's arc-wise constant empirical acceleration terms.
	"""

def arcwise_empirical_accelerations(body: str, centralBody: str, acceleration_components: dict[EmpiricalAccelerationComponents, list[EmpiricalAccelerationFunctionalShapes]], arc_start_times: list[float]) -> EstimatableParameterSettings:
    """Function for defining parameter settings for arc-wise empirical acceleration magnitudes.
	
	Factory function for creating a (linear sensitivity) parameter settings object for arc-wise empirical acceleration magnitudes (arc-wise version of :func:``~tudatpy.numerical_simulation.estimation_setup.parameter.empirical_accelerations`).
	Using the empirical acceleration terms as estimatable parameters requires
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.empirical` acceleration, which include constant (in RSW frame) terms
	
	Note: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the radiation pressure coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.
	
	
	:param body:
			Name of the body, with whose empirical acceleration model the estimatable parameter is associated.
	
	:param centralBody:
			Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration
	
	:param acceleration_components:
			Dictionary of components of the empirical acceleration which are to be estimated. There are two 'degrees of freedom' in these components: the direction of the acceleration (e.g. R, S or W direction) and the temporal signature (constant, sine of true anomaly or cosine of true anomaly). With this input, any subset may be selected. This parameter is a dictionary, with the key denoting the direction of the acceleration, and the value a list of the temporal signatures to estimate for this empirical acceleration direction.
	
	:param arc_initial_times:
			List of times at which the arcs over which the empirical accelerations are to be estimated will start.
	
	:return:
			Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EmpiricalAccelerationEstimatableParameterSettings` class
			for the specified body's arc-wise empirical acceleration terms.
	"""

def arcwise_radiation_pressure_coefficient(body: str, arc_initial_times: list[float]) -> EstimatableParameterSettings:
    """Function for defining parameter settings for arc-wise radiation pressure coefficients.
	
	Factory function for creating a (linear sensitivity) parameter settings object for arc-wise radiation pressure coefficients (arc-wise version of :func:``~tudatpy.numerical_simulation.estimation_setup.parameter.radiation_pressure_coefficient`).
	Using the radiation pressure coefficient as an estimatable parameter requires
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.cannonball` radiation pressure interface to be defined for the body specified by the ``body`` parameter
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.cannonball_radiation_pressure` acceleration
	The radiation pressure coefficient is defined according to the universal convention for a cannonball model and thus no further model information is given.
	
	Note: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the radiation pressure coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.
	
	
	:param body:
			Name of the body, with whose radiation pressure acceleration model the estimatable parameter is associated.
	
	:param arc_initial_times:
			List of times at which the arcs over which the radiation pressure coefficient is to be estimated will start.
	
	:return:
			Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ArcWiseRadiationPressureCoefficientEstimatableParameterSettings` class
			for arc-wise treatment of the specified body's radiation pressure coefficient.
	"""

def arcwise_relative_observation_bias(link_ends: ..., observable_type: ..., arc_start_times: list[float], time_link_end: ...) -> EstimatableParameterSettings:
    """Function for defining parameter settings for arc-wise absolute observation bias.
	
	Factory function for creating a (linear sensitivity) parameter settings object for the arc-wise treatment of an observation's relative bias parameter.
	Using the arc-wise relative observation bias as estimatable parameter requires
	  * The observation model (corresponding to the `link_ends` and `observable_type`) to include an arc-wise relative bias (:func:`~tudatpy.numerical_simulation.estimation_setup.observation.arcwise_relative_bias`)
	Note: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the bias may, but need not, correspond to the arcs used for a multi-arc propagation.
	
	
	:param link_ends:
			Set of link ends that define the geometry of the biased observations.
	
	:param observable_type:
			Observable type of the biased observations.
	:param arc_start_times:
			List of times at which the arcs over which the bias is to be estimated will start.
	:param time_link_end:
			The link end type (transmitter, receiver, etc.) at which the arc_start_times is evaluated.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ArcWiseConstantObservationBiasEstimatableParameterSettings`
			for the specified observation's arc-wise relative bias.
	"""

def arcwise_time_drift_observation_bias(link_ends: dict[..., ...], observable_type: ..., arc_start_times: list[float], ref_epochs: list[float], time_link_end: ...) -> EstimatableParameterSettings:
    ...

def constant_drag_coefficient(body: str) -> EstimatableParameterSettings:
    """Function for defining parameter settings for constant drag coefficients.
	
	Factory function for creating a (linear sensitivity) parameter settings object for constant drag coefficients.
	Using the constant drag coefficient as an estimatable parameter requires:
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.constant` aerodynamic interface to be defined for the body specified by the ``body`` parameter
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.aerodynamic` acceleration
	
	
	:param body:
			Name of the body, with whose drag acceleration model the estimatable parameter is associated.
	
	:return:
			:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's constant drag coefficient.
	"""

def constant_empirical_acceleration_terms(body: str, centralBody: str) -> EstimatableParameterSettings:
    """Function for defining parameter settings for constant empirical acceleration terms.
	
	As :func:`~tudatpy.numerical_simulation.estimation_setup.parameter.empirical_accelerations`, but only using the constant R, S and W components (no sine or cosine term estimation). This function is added as a function of convenience
	
	
	:param body:
			Name of the body, with whose empirical acceleration model the estimatable parameter is associated.
	
	:param centralBody:
			Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration
	
	:return:
			Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EmpiricalAccelerationEstimatableParameterSettings` class
			for the specified body's empirical acceleration terms.
	"""

def constant_rotation_rate(body: str) -> EstimatableParameterSettings:
    """Function for defining parameter settings for a body's constant rotation rate.
	
	Factory function for creating a (linear sensitivity) parameter settings object for a body's constant rotation rate parameter.
	Using the constant rotation rate as estimatable parameter requires:
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple` or :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple_from_spice` rotation model specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter
	
	
	:param body:
			Name of the body, with whose rotation model the estimatable parameter is associated.
	
	:return:
			:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's constant spin rate.
	"""

def core_factor(body: str) -> EstimatableParameterSettings:
    """Function for defining parameter settings for a body's core factor.
	
	Factory function for creating a (linear sensitivity) parameter settings object for a body's core factor.
	Using the core factor as estimatable parameter requires
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter
	
	
	:param body:
			Name of the body, with whose rotation model the estimatable parameter is associated.
	
	:return:
			:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's core factor.
	"""

def custom_analytical_partial(analytical_partial_function: typing.Callable[[float, numpy.ndarray], numpy.ndarray], body_undergoing_acceleration: str, body_exerting_acceleration: str, acceleration_type: ...) -> CustomAccelerationPartialSettings:
    ...

def custom_numerical_partial(parameter_perturbation: numpy.ndarray, body_undergoing_acceleration: str, body_exerting_acceleration: str, acceleration_type: ..., environment_updates: dict[..., list[str]]={}) -> CustomAccelerationPartialSettings:
    ...

def custom_parameter(custom_id: str, parameter_size: int, get_parameter_function: typing.Callable[[], numpy.ndarray], set_parameter_function: typing.Callable[[numpy.ndarray], None]) -> EstimatableParameterSettings:
    ...

@typing.overload
def direct_tidal_dissipation_time_lag(body: str, deforming_body: str) -> EstimatableParameterSettings:
    ...

@typing.overload
def direct_tidal_dissipation_time_lag(body: str, deforming_body: list[str]) -> EstimatableParameterSettings:
    ...

def empirical_accelerations(body: str, centralBody: str, acceleration_components: dict[EmpiricalAccelerationComponents, list[EmpiricalAccelerationFunctionalShapes]]) -> EstimatableParameterSettings:
    """Function for defining parameter settings for empirical acceleration magnitudes.
	
	Factory function for creating a (linear sensitivity) parameter settings object for empirical acceleration magnitudes.
	Using the empirical acceleration terms as estimatable parameters requires
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.empirical` acceleration, which include constant (in RSW frame) terms
	
	
	:param body:
			Name of the body, with whose empirical acceleration model the estimatable parameter is associated.
	
	:param centralBody:
			Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration
	
	:param acceleration_components:
			Dictionary of components of the empirical acceleration which are to be estimated. There are two 'degrees of freedom' in these components: the direction of the acceleration (e.g. R, S or W direction) and the temporal signature (constant, sine of true anomaly or cosine of true anomaly). With this input, any subset may be selected. This parameter is a dictionary, with the key denoting the direction of the acceleration, and the value a list of the temporal signatures to estimate for this empirical acceleration direction.
	
	:return:
			Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EmpiricalAccelerationEstimatableParameterSettings` class
			for the specified body's empirical acceleration terms.
	"""

def free_core_nutation_rate(body: str) -> EstimatableParameterSettings:
    """Function for defining parameter settings for a body's free core nutation rate.
	
	Factory function for creating a (linear sensitivity) parameter settings object for a body's free core nutation rate.
	Using the free core nutation rate as estimatable parameter requires
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter
	
	
	:param body:
			Name of the body, with whose rotation model the estimatable parameter is associated.
	
	:return:
			:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's free core nutation rate.
	"""

def full_empirical_acceleration_terms(body: str, centralBody: str) -> EstimatableParameterSettings:
    """Function for defining parameter settings for empirical acceleration magnitudes.
	
	Factory function for creating a (linear sensitivity) parameter settings object for empirical acceleration magnitudes.
	Using the empirical acceleration terms as estimatable parameters requires
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.empirical` acceleration, which include constant (in RSW frame) terms
	
	
	:param body:
			Name of the body, with whose empirical acceleration model the estimatable parameter is associated.
	
	:param centralBody:
			Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration
	
	:param acceleration_components:
			Dictionary of components of the empirical acceleration which are to be estimated. There are two 'degrees of freedom' in these components: the direction of the acceleration (e.g. R, S or W direction) and the temporal signature (constant, sine of true anomaly or cosine of true anomaly). With this input, any subset may be selected. This parameter is a dictionary, with the key denoting the direction of the acceleration, and the value a list of the temporal signatures to estimate for this empirical acceleration direction.
	
	:return:
			Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EmpiricalAccelerationEstimatableParameterSettings` class
			for the specified body's empirical acceleration terms.
	"""

def gravitational_parameter(body: str) -> EstimatableParameterSettings:
    """Function for defining parameter settings for a massive body's gravitational parameter.
	
	Factory function for creating a (linear sensitivity) parameter settings object for the gravitational parameter of massive bodies.
	Using the gravitational parameter as estimatable parameter requires
	  * The body specified by the ``body`` parameter to be endowed with a gravity field (see :ref:`\\`\\`gravity_field\\`\\`` module for options)
	  * Any dynamical or observational model to depend on the gravitational parameter of the body specified by the ``body`` parameter
	
	
	:param body:
			Name of the body, with whose gravitational model the estimatable parameter is associated.
	
	:return:
			:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's gravitational parameter.
	"""

def ground_station_position(body: str, ground_station_name: str) -> EstimatableParameterSettings:
    """Function for defining parameter settings for ground station position bias.
	
	Factory function for creating a (linear sensitivity) parameter settings object for a ground station's body-fixed Cartesian position.
	Using the ground station position bias as estimatable parameter requires:
	  * At least one observation model to rely on the specified ground station
	
	
	:param body:
			Body name identifying the body, with which the ground station is associated.
	:param ground_station_name:
			Name which identifies the position-biased ground station.
	:return:
			:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified ground station's position bias.
	"""

def initial_states(propagator_settings: ..., bodies: ..., arc_initial_times: list[float]=[]) -> list[EstimatableParameterSettings]:
    """Function for defining parameter settings for initial state parameters.
	
	Factory function for creating a (linear sensitivity) parameter settings object for initial state parameters.
	The factory function uses the propagator settings to determine which type of initial state parameter (single/multi/hybrid-arc; translational/rotational/... dynamics) is to be estimated,
	e.g. if a single-arc translational state propagator is defined, the function will automatically create the parameters for the associated initial state parameter
	
	Note: These function return lists of parameter settings objects.
	This means that which the return of this function cannot simply be added to the parameter settings objects of single parameters in a list creation statement.
	Instead, list concatenation is recommended. Please see the following example:
	
	.. code-block:: python
	
	   # define single estimatable parameters
	   single_parameter_1 = ...
	   single_parameter_2 = ...
	   ...
	
	   # bad: list creation statement --> will result in nested list, undesired!
	   list_of_all_parameters = [estimation_setup.parameter.initial_states(...), single_parameter_1, single_parameter_2, ...]
	
	   # better: list concatenation --> will result in simple list, desired!
	   list_of_all_parameters = estimation_setup.parameter.initial_states(...) + [single_parameter_1, single_parameter_2, ...]
	
	
	:param propagator_settings:
			Object containing the consolidated propagation settings of the simulation in the context of which the given model parameters are to be estimated.
	
	:param bodies:
			Object consolidating all bodies and environment models that constitute the physical environment.
	
	:param arc_initial_times:
			Initial times of arcs, only required if arc-wise propagation settings are passed via the `propagator_settings` object.
	
	:return:
			List of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` objects, one per component of each initial state in the simulation.
	"""

@typing.overload
def inverse_tidal_quality_factor(body: str, deforming_body: str) -> EstimatableParameterSettings:
    ...

@typing.overload
def inverse_tidal_quality_factor(body: str, deforming_body: list[str]) -> EstimatableParameterSettings:
    ...

def mean_moment_of_inertia(body: str) -> EstimatableParameterSettings:
    """Function for defining parameter settings for a body's mean moment of inertia.
	
	Factory function for creating a (linear sensitivity) parameter settings object for a body's mean moment of inertia.
	In most cases, the mean moment of inertia will not influence the dynamics/observation directly and sensitivity to this parameter will not be included. The dynamics/observation will be sensitive to this parameter if the rotational dynamics of a relevant body is estimated.
	Using the mean moment of inertia as estimatable parameter requires
	  * The estimation of an initial rotational state of the body specified by the ``body`` parameter
	
	
	:param body:
			Name of the body, with whose body model the estimatable parameter is associated.
	
	:return:
			:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's mean moment of inertia.
	"""

def monomial_full_block_gravity_field_variation_amplitudes(body_name: str, power: int, minimum_degree: int, minimum_order: int, maximum_degree: int, maximum_order: int) -> EstimatableParameterSettings:
    ...

def monomial_gravity_field_variation_amplitudes(body_name: str, power: int, cosine_indices: list[tuple[int, int]], sine_indices: list[tuple[int, int]]) -> EstimatableParameterSettings:
    ...

@typing.overload
def order_invariant_k_love_number(deformed_body: str, degree: int, deforming_body: str, use_complex_love_number: bool=0) -> EstimatableParameterSettings:
    ...

@typing.overload
def order_invariant_k_love_number(deformed_body: str, degree: int, deforming_bodies: list[str], use_complex_love_number: bool=0) -> EstimatableParameterSettings:
    ...

@typing.overload
def order_invariant_k_love_number(deformed_body: str, degree: int, use_complex_love_number: bool=0) -> EstimatableParameterSettings:
    ...

@typing.overload
def order_varying_k_love_number(deformed_body: str, degree: int, orders: list[int], deforming_body: str, use_complex_love_number: bool=0) -> EstimatableParameterSettings:
    ...

@typing.overload
def order_varying_k_love_number(deformed_body: str, degree: int, orders: list[int], deforming_bodies: list[str], use_complex_love_number: bool=0) -> EstimatableParameterSettings:
    ...

@typing.overload
def order_varying_k_love_number(deformed_body: str, degree: int, orders: list[int], use_complex_love_number: bool=0) -> EstimatableParameterSettings:
    ...

def periodic_spin_variations(body: str) -> EstimatableParameterSettings:
    """Function for defining parameter settings for a body's periodic spin variations.
	
	Factory function for creating a (linear sensitivity) parameter settings object for a body's periodic spin variation parameters.
	Using the mean moment of inertia as estimatable parameter requires
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter
	
	
	:param body:
			Name of the body, with whose rotation model the estimatable parameter is associated.
	
	:return:
			:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's periodic spin variations.
	"""

def polar_motion_amplitudes(body: str) -> EstimatableParameterSettings:
    """Function for defining parameter settings for a body's polar motion amplitudes.
	
	Factory function for creating a (linear sensitivity) parameter settings object for a body's polar motion amplitudes.
	Using the polar motion amplitudes as estimatable parameter requires
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter
	
	
	:param body:
			Name of the body, with whose rotation model the estimatable parameter is associated.
	
	:return:
			:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's polar motion amplitudes.
	"""

def polynomial_gravity_field_variation_amplitudes(body_name: str, cosine_indices_per_power: dict[int, list[tuple[int, int]]], sine_indices_per_power: dict[int, list[tuple[int, int]]]) -> EstimatableParameterSettings:
    ...

def ppn_parameter_beta() -> EstimatableParameterSettings:
    """Function for defining parameter settings for post-newtonian beta parameter.
	
	Factory function for creating a (linear sensitivity) parameter settings object for a global PPN :math:`\x08eta` parameter.
	Using the post-newtonian gamma parameter as estimatable parameter requires at least one of the following:
	 * An acceleration model depending on this parameter, such as :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.relativistic_correction`
	 * An observation model with a light-time correction depending on this parameter (none yet implemented)
	
	:return:
			:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for a global post-newtonian :math:`\x08eta` parameter.
	"""

def ppn_parameter_gamma() -> EstimatableParameterSettings:
    """Function for defining parameter settings for post-newtonian gamma parameter.
	
	Factory function for creating a (linear sensitivity) parameter settings object for a global PPN :math:`\\gamma` parameter.
	Using the post-newtonian gamma parameter as estimatable parameter requires at least one of the following:
	 * An acceleration model depending on this parameter, such as :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.relativistic_correction`
	 * An observation model with a light-time correction depending on this parameter, such as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.first_order_relativistic_light_time_correction`
	
	:return:
			:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for a global post-newtonian :math:`\\gamma` parameter.
	"""

def quasi_impulsive_shots(body: str) -> EstimatableParameterSettings:
    """Function for defining parameter settings for quasi-impulsive shots.
	
	Factory function for creating a (linear sensitivity) parameter settings object for so-called 'quasi-impulsive shots', such as desaturation maneuvers. With this parameter, the total :math:`\\Delta V` vector of a set of such shots/maneuvers can be estimated.
	Using the quasi-impulsive shots as an estimatable parameter requires
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.quasi_impulsive_shots_acceleration` acceleration
	Note: this parameter considers *all* shots/maneuvers used in the above acceleration model, and estimates the value of the 'delta_v_values' input of this acceleration.
	
	
	:param body:
			Name of the body, with which the quasi-impulsive shot estimatable parameter is associated.
	
	:return:
			:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's quasi-impulsive shots
	"""

def radiation_pressure_coefficient(body: str) -> EstimatableParameterSettings:
    """Function for defining parameter settings for radiation pressure coefficients.
	
	Factory function for creating a (linear sensitivity) parameter settings object for a radiation pressure coefficient.
	Using the radiation pressure coefficient as an estimatable parameter requires
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.cannonball` radiation pressure interface to be defined for the body specified by the ``body`` parameter
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.cannonball_radiation_pressure` acceleration
	
	
	:param body:
			Name of the body, with whose radiation pressure acceleration model the estimatable parameter is associated.
	
	:return:
			:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's radiation pressure coefficient.
	"""

def radiation_pressure_target_direction_scaling(target_body: str, source_body: str) -> EstimatableParameterSettings:
    ...

def radiation_pressure_target_perpendicular_direction_scaling(target_body: str, source_body: str) -> EstimatableParameterSettings:
    ...

def relative_observation_bias(link_ends: ..., observable_type: ...) -> EstimatableParameterSettings:
    """Function for defining parameter settings for an relative observation bias.
	
	Factory function for creating a (linear sensitivity) parameter settings object for an observation's relative bias parameter.
	Using the relative observation bias as estimatable parameter requires
	  * The observation model (corresponding to the `link_ends` and `observable_type`) to include a relative bias (:func:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias`)
	Note: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the bias may, but need not, correspond to the arcs used for a multi-arc propagation.
	
	
	:param link_ends:
			Set of link ends that define the geometry of the biased observations.
	
	:param observable_type:
			Observable type of the biased observations.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ConstantObservationBiasEstimatableParameterSettings`
			for the specified observation's arc-wise relative bias.
	"""

def rotation_pole_position(body: str) -> EstimatableParameterSettings:
    """Function for defining parameter settings for a body's rotation pole position.
	
	Factory function for creating a (linear sensitivity) parameter settings object for a body's rotation pole position, parameterized by the constant pole rotation angles (:math:`alpha` and :math:`\\delta`).
	Using the rotation pole position as estimatable parameter requires:
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple` or :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple_from_spice` rotation model specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter
	
	
	:param body:
			Name of the body, with whose rotation model the estimatable parameter is associated.
	
	:return:
			:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's rotation pole position.
	"""

def scaled_longitude_libration_amplitude(body_name: str) -> EstimatableParameterSettings:
    ...

def spherical_harmonics_c_coefficients(body: str, minimum_degree: int, minimum_order: int, maximum_degree: int, maximum_order: int) -> EstimatableParameterSettings:
    """Function for defining parameter settings for the cosine coefficients of body's spherical harmonics gravitational model.
	
	Factory function for creating a (linear sensitivity) parameter settings object for the spherical harmonics cosine-coefficients (:math:`\x08ar{C}_{lm}`) of a hody with a spherical harmonic gravity field. Using this function, a 'full' set of spherical harmonic coefficients between an minimum/maximum degree/order are estimated. For instance, for minimum degree/order of 2/0, and maximum degree/order 4/4, all spherical harmonic cosine coefficients of degrees 2, 3 and 4 are estimated. If the maximum degree/order is set to 4/2, only coefficients with an order of 0, 1 and 2 are included. The entries in the parameter are sorted first by degree, and then by order (both in ascending order)
	Using the spherical harmonics cosine coefficients as estimatable parameter requires:
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic` (or derived) gravity model to be defined for the body specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the estimated cosine coefficients of the body specified by the ``body`` parameter. Typically, this dependency will be a :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.spherical_harmonic` acceleration
	
	
	:param body:
			Name of the body, with whose gravitational model the estimatable parameters are associated.
	
	:param minimum_degree:
			Minimum degree of c-coefficients to be included.
	:param minimum_order:
			Minimum order of c-coefficients to be included.
	:param maximum_degree:
			Maximum degree of c-coefficients to be included.
	:param maximum_order:
			Maximum order of c-coefficients to be included.
	:return:
			Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.SphericalHarmonicEstimatableParameterSettings` class
			for the applicable spherical harmonics c-coefficients of the specified body's gravitational model.
	"""

def spherical_harmonics_c_coefficients_block(body: str, block_indices: list[tuple[int, int]]) -> EstimatableParameterSettings:
    """Function for defining parameter settings for the cosine coefficients of body's spherical harmonics gravitational model.
	
	As :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.spherical_harmonics_c_coefficients`, but with a manually defined set of coefficients.
	
	
	:param body:
			Name of the body, with whose gravitational model the estimatable parameters are associated.
	
	:param block_indices:
			List of block indices. The length of this list can be arbitrary, as long as the pairs are unique.
			For each pair, the first value is the degree and the second the order of the coefficient to be included.
	
	:return:
			Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.SphericalHarmonicEstimatableParameterSettings` class
			for the applicable spherical harmonics c-coefficients of the specified body's gravitational model.
	"""

def spherical_harmonics_s_coefficients(body: str, minimum_degree: int, minimum_order: int, maximum_degree: int, maximum_order: int) -> EstimatableParameterSettings:
    """Function for defining parameter settings for the sine coefficients of body's spherical harmonics gravitational model.
	
	Factory function for creating a (linear sensitivity) parameter settings object for the spherical harmonics sine-coefficients (:math:`\x08ar{S}_{lm}`) of a hody with a spherical harmonic gravity field. Using this function, a 'full' set of spherical harmonic coefficients between an minimum/maximum degree/order are estimated. For instance, for minimum degree/order of 2/1 (there is no order 0 sine coefficient), and maximum degree/order 4/4, all spherical harmonic sine coefficients of degrees 2, 3 and 4 are estimated. If the maximum degree/order is set to 4/2, only coefficients with an order of 1 and 2 are included. The entries in the parameter are sorted first by degree, and then by order (both in ascending order)
	Using the spherical harmonics cosine coefficients as estimatable parameter requires:
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic` (or derived) gravity model to be defined for the body specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the estimated cosine coefficients of the body specified by the ``body`` parameter. Typically, this dependency will be a :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.spherical_harmonic` acceleration
	
	
	:param body:
			Name of the body, with whose gravitational model the estimatable parameters are associated.
	
	:param minimum_degree:
			Minimum degree of s-coefficients to be included.
	:param minimum_order:
			Minimum order of s-coefficients to be included.
	:param maximum_degree:
			Maximum degree of s-coefficients to be included.
	:param maximum_order:
			Maximum order of s-coefficients to be included.
	:return:
			Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.SphericalHarmonicEstimatableParameterSettings` class
			for the applicable spherical harmonics s-coefficients of the specified body's gravitational model.
	"""

def spherical_harmonics_s_coefficients_block(body: str, block_indices: list[tuple[int, int]]) -> EstimatableParameterSettings:
    """Function for defining parameter settings for the sine coefficients of body's spherical harmonics gravitational model.
	
	As :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.spherical_harmonics_s_coefficients`, but with a manually defined set of coefficients.
	
	
	:param body:
			Name of the body, with whose gravitational model the estimatable parameters are associated.
	
	:param block_indices:
			List of block indices. The length of this list can be arbitrary, as long as the pairs are unique.
			For each pair, the first value is the degree and the second the order of the coefficient to be included.
	
	:return:
			Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.SphericalHarmonicEstimatableParameterSettings` class
			for the applicable spherical harmonics s-coefficients of the specified body's gravitational model.
	"""

def time_drift_observation_bias(link_ends: dict[..., ...], observable_type: ..., ref_epoch: float, time_link_end: ...) -> EstimatableParameterSettings:
    ...

def yarkovsky_parameter(body_name: str, central_body_name: str='Sun') -> EstimatableParameterSettings:
    ...
across_track_empirical_acceleration_component: EmpiricalAccelerationComponents
along_track_empirical_acceleration_component: EmpiricalAccelerationComponents
arc_wise_constant_drag_coefficient_type: EstimatableParameterTypes
arc_wise_empirical_acceleration_coefficients_type: EstimatableParameterTypes
arc_wise_initial_body_state_type: EstimatableParameterTypes
arc_wise_radiation_pressure_coefficient_type: EstimatableParameterTypes
arc_wise_time_drift_observation_bias_type: EstimatableParameterTypes
arcwise_constant_additive_observation_bias_type: EstimatableParameterTypes
arcwise_constant_relative_observation_bias_type: EstimatableParameterTypes
constant_additive_observation_bias_type: EstimatableParameterTypes
constant_drag_coefficient_type: EstimatableParameterTypes
constant_empirical: EmpiricalAccelerationFunctionalShapes
constant_relative_observation_bias_type: EstimatableParameterTypes
constant_rotation_rate_type: EstimatableParameterTypes
constant_time_drift_observation_bias_type: EstimatableParameterTypes
core_factor_type: EstimatableParameterTypes
cosine_empirical: EmpiricalAccelerationFunctionalShapes
desaturation_delta_v_values_type: EstimatableParameterTypes
direct_dissipation_tidal_time_lag_type: EstimatableParameterTypes
empirical_acceleration_coefficients_type: EstimatableParameterTypes
equivalence_principle_lpi_violation_parameter_type: EstimatableParameterTypes
free_core_nutation_rate_type: EstimatableParameterTypes
full_degree_tidal_love_number_type: EstimatableParameterTypes
gravitational_parameter_type: EstimatableParameterTypes
ground_station_position_type: EstimatableParameterTypes
initial_body_state_type: EstimatableParameterTypes
initial_rotational_body_state_type: EstimatableParameterTypes
inverse_tidal_quality_factor_type: EstimatableParameterTypes
mean_moment_of_inertia_type: EstimatableParameterTypes
periodic_spin_variation_type: EstimatableParameterTypes
polar_motion_amplitude_type: EstimatableParameterTypes
ppn_parameter_beta_type: EstimatableParameterTypes
ppn_parameter_gamma_type: EstimatableParameterTypes
radial_empirical_acceleration_component: EmpiricalAccelerationComponents
radiation_pressure_coefficient_type: EstimatableParameterTypes
rotation_pole_position_type: EstimatableParameterTypes
sine_empirical: EmpiricalAccelerationFunctionalShapes
single_degree_variable_tidal_love_number_type: EstimatableParameterTypes
spherical_harmonics_cosine_coefficient_block_type: EstimatableParameterTypes
spherical_harmonics_sine_coefficient_block_type: EstimatableParameterTypes