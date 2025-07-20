import numpy
import pybind11_stubgen.typing_ext
from ...dynamics import environment
from ...dynamics.propagation_setup import acceleration
from ...dynamics.propagation_setup import propagator
import typing
__all__ = ['CustomAccelerationPartialSettings', 'EmpiricalAccelerationComponents', 'EmpiricalAccelerationFunctionalShapes', 'EstimatableParameterSettings', 'EstimatableParameterTypes', 'absolute_observation_bias', 'across_track_empirical_acceleration_component', 'along_track_empirical_acceleration_component', 'arc_wise_constant_drag_coefficient_type', 'arc_wise_empirical_acceleration_coefficients_type', 'arc_wise_initial_body_state_type', 'arc_wise_polynomial_clock_corrections', 'arc_wise_polynomial_clock_corrections_type', 'arc_wise_radiation_pressure_coefficient_type', 'arc_wise_time_drift_observation_bias_type', 'arcwise_absolute_observation_bias', 'arcwise_constant_additive_observation_bias_type', 'arcwise_constant_drag_coefficient', 'arcwise_constant_empirical_acceleration_terms', 'arcwise_constant_relative_observation_bias_type', 'arcwise_empirical_accelerations', 'arcwise_radiation_pressure_coefficient', 'arcwise_relative_observation_bias', 'arcwise_time_bias', 'arcwise_time_drift_observation_bias', 'constant_additive_observation_bias_type', 'constant_drag_coefficient', 'constant_drag_coefficient_type', 'constant_empirical', 'constant_empirical_acceleration_terms', 'constant_relative_observation_bias_type', 'constant_rotation_rate', 'constant_rotation_rate_type', 'constant_time_bias', 'constant_time_drift_observation_bias_type', 'core_factor', 'core_factor_type', 'cosine_empirical', 'create_parameter_set', 'custom_analytical_partial', 'custom_numerical_partial', 'custom_parameter', 'desaturation_delta_v_values_type', 'direct_dissipation_tidal_time_lag_type', 'direct_tidal_dissipation_time_lag', 'drag_component_scaling', 'empirical_acceleration_coefficients_type', 'empirical_accelerations', 'equivalence_principle_lpi_violation_parameter_type', 'free_core_nutation_rate', 'free_core_nutation_rate_type', 'full_degree_tidal_love_number_type', 'full_empirical_acceleration_terms', 'global_polynomial_clock_corrections', 'global_polynomial_clock_corrections_type', 'gravitational_parameter', 'gravitational_parameter_type', 'ground_station_position', 'ground_station_position_type', 'initial_body_state_type', 'initial_rotational_body_state_type', 'initial_states', 'inverse_tidal_quality_factor', 'inverse_tidal_quality_factor_type', 'lift_component_scaling', 'mean_moment_of_inertia', 'mean_moment_of_inertia_type', 'mode_coupled_k_love_numbers', 'monomial_full_block_gravity_field_variation_amplitudes', 'monomial_gravity_field_variation_amplitudes', 'order_invariant_k_love_number', 'order_varying_k_love_number', 'periodic_gravity_field_variation_amplitudes', 'periodic_spin_variation_type', 'periodic_spin_variations', 'polar_motion_amplitude_type', 'polar_motion_amplitudes', 'polynomial_gravity_field_variation_amplitudes', 'ppn_parameter_beta', 'ppn_parameter_beta_type', 'ppn_parameter_gamma', 'ppn_parameter_gamma_type', 'quasi_impulsive_shots', 'radial_empirical_acceleration_component', 'radiation_pressure_coefficient', 'radiation_pressure_coefficient_type', 'radiation_pressure_target_direction_scaling', 'radiation_pressure_target_perpendicular_direction_scaling', 'reference_point_position', 'relative_observation_bias', 'rotation_pole_position', 'rotation_pole_position_type', 'scaled_longitude_libration_amplitude', 'side_component_scaling', 'sine_empirical', 'single_degree_variable_tidal_love_number_type', 'spherical_harmonics_c_coefficients', 'spherical_harmonics_c_coefficients_block', 'spherical_harmonics_cosine_coefficient_block_type', 'spherical_harmonics_s_coefficients', 'spherical_harmonics_s_coefficients_block', 'spherical_harmonics_sine_coefficient_block_type', 'time_drift_observation_bias', 'yarkovsky_parameter']

class CustomAccelerationPartialSettings:
    """No documentation found."""

class EmpiricalAccelerationComponents:
    """Enumeration of the available empirical acceleration components that are available to estimate.
                
                These are used in the :func:`~tudatpy.dynamics.parameters_setup.empirical_accelerations` function to specify which components of the empirical acceleration are to be estimated.
                
    
    Members:
    
      radial_empirical_acceleration_component
    
      along_track_empirical_acceleration_component
    
      across_track_empirical_acceleration_component"""
    __members__: typing.ClassVar[dict[str, EmpiricalAccelerationComponents]]
    across_track_empirical_acceleration_component: typing.ClassVar[EmpiricalAccelerationComponents]
    along_track_empirical_acceleration_component: typing.ClassVar[EmpiricalAccelerationComponents]
    radial_empirical_acceleration_component: typing.ClassVar[EmpiricalAccelerationComponents]

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
    """Enumeration of the available empirical acceleration shapes that are available per component
                
                These are used in the :func:`~tudatpy.dynamics.parameters_setup.empirical_accelerations` function to specify the signature of the estimated empirical acceleration component.
                .
    
    Members:
    
      constant_empirical
    
      sine_empirical
    
      cosine_empirical"""
    __members__: typing.ClassVar[dict[str, EmpiricalAccelerationFunctionalShapes]]
    constant_empirical: typing.ClassVar[EmpiricalAccelerationFunctionalShapes]
    cosine_empirical: typing.ClassVar[EmpiricalAccelerationFunctionalShapes]
    sine_empirical: typing.ClassVar[EmpiricalAccelerationFunctionalShapes]

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
    Settings of simple parameters types are managed via this class, more complex parameter types are handled by specialised derivatives of this class.
    Instances of either base or derived class can be created via dedicated functions."""
    custom_partial_settings: list[...]

class EstimatableParameterTypes:
    """Enumeration of model parameters that are available for estimation.
             In order to establish a parameter estimation settings for a parameter of a certain type, use the function dedicated to this parameter type.
             Note that not all of the listed types might be accessible via functions in the python interface yet.
    
    
    
    
    
    
          
    
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
    
      global_polynomial_clock_corrections_type
    
      arc_wise_polynomial_clock_corrections_type
    
      inverse_tidal_quality_factor_type"""
    __members__: typing.ClassVar[dict[str, EstimatableParameterTypes]]
    arc_wise_constant_drag_coefficient_type: typing.ClassVar[EstimatableParameterTypes]
    arc_wise_empirical_acceleration_coefficients_type: typing.ClassVar[EstimatableParameterTypes]
    arc_wise_initial_body_state_type: typing.ClassVar[EstimatableParameterTypes]
    arc_wise_polynomial_clock_corrections_type: typing.ClassVar[EstimatableParameterTypes]
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
    global_polynomial_clock_corrections_type: typing.ClassVar[EstimatableParameterTypes]
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
    """Function for creating parameter settings for an absolute observation bias.
    
    Function for creating parameter settings object for an observation's absolute bias parameter.
    Using the absolute observation bias as estimatable parameter requires:
    
    * The observation model (corresponding to the `link_ends` and `observable_type`) to include an absolute bias (:func:`~tudatpy.estimation.observable_models_setup.biases.absolute_bias`)
    
    
    Parameters
    ----------
    link_ends : Dict[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, Tuple[str, str]
        Set of link ends that define the geometry of the biased observations.
    
    observable_type : ObservableType
        Observable type of the biased observations.
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of the :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.ConstantObservationBiasEstimatableParameterSettings`
        for the specified observation's arc-wise absolute bias."""

def arc_wise_polynomial_clock_corrections(associated_body: str, associated_station: str, correction_powers: list[int], arc_indices: list[int]) -> EstimatableParameterSettings:
    ...

def arcwise_absolute_observation_bias(link_ends: ..., observable_type: ..., arc_start_times: list[float], time_link_end: ...) -> EstimatableParameterSettings:
    """Function for creating parameter settings for arc-wise absolute observation bias.
    
    Function for creating parameter settings object for the arc-wise treatment of an observation's absolute bias parameter.
    Using the arc-wise absolute observation bias as estimatable parameter requires
    
    * The observation model (corresponding to the `link_ends` and `observable_type`) to include an arc-wise absolute bias (:func:`~tudatpy.estimation.observable_models_setup.biases.arcwise_absolute_bias`)
    
    
    Parameters
    ----------
    link_ends : Dict[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, Tuple[str, str]
        Set of link ends that define the geometry of the biased observations.
    
    observable_type : ObservableType
        Observable type of the biased observations.
    arc_start_times : List[ float ]
        List of times at which the arcs over which the bias is to be estimated will start.
    time_link_end : LinkEndType
        The link end type (transmitter, receiver, etc.) at which the arc_start_times is evaluated.
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of the :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.ArcWiseConstantObservationBiasEstimatableParameterSettings`
        for the specified observation's arc-wise absolute bias."""

def arcwise_constant_drag_coefficient(body: str, arc_initial_times: list[float]) -> EstimatableParameterSettings:
    """Function for creating parameter settings for arc-wise constant drag coefficients.
    
    Function for creating parameter settings object for arc-wise constant drag coefficients :math:`C_{D}` 
    (arc-wise version of :func:`~tudatpy.dynamics.parameters_setup.constant_drag_coefficient`).
    Using the arc-wise constant drag coefficient as an estimatable parameter requires:
    
    * A :func:`~tudatpy.dynamics.environment_setup.aerodynamic_coefficients.constant` aerodynamic interface to be defined for the body specified by the ``body`` parameter
    * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.aerodynamic` acceleration
    
    When using this parameter, whenever :math:`C_{D}` is required at a time :math:`t`, the index math:`i` in the ``arc_initial_times`` ordered list is
    found for which :math:`t_{i}\\le t<t_{i+1}` (or, if :math:`t` is larger than the largest value in the list, :math:`i` is set to be last index of the list),
    and the parameter entry representing :math:`C_{D,i}` will be used.
    
    .. note:: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the drag coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose drag acceleration model the estimatable parameter is associated.
    arc_initial_times : List[ float ]
        Ordered list of times at which the arcs over which the drag coefficient is to be estimated will start.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.ArcWiseDragCoefficientEstimatableParameterSettings` class
        for arc-wise treatment of the specified body's constant drag coefficient."""

def arcwise_constant_empirical_acceleration_terms(body: str, centralBody: str, arc_start_times: list[float]) -> EstimatableParameterSettings:
    """Function for creating parameter settings for arc-wise constant empirical acceleration terms.
    
    As :func:`~tudatpy.dynamics.parameters_setup.arcwise_empirical_accelerations`, but only using the constant R, S and W components (no sine or cosine term estimation). This function is added as a function of convenience
    
    Parameters
    ----------
    body : str
        Name of the body, with whose empirical acceleration model the estimatable parameter is associated.
    
    centralBody : str
        Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration
    
    arc_initial_times : List[ float ]
        List of times at which the arcs over which the empirical accelerations are to be estimated will start.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.EmpiricalAccelerationEstimatableParameterSettings` class
        for the specified body's arc-wise constant empirical acceleration terms."""

def arcwise_empirical_accelerations(body: str, centralBody: str, acceleration_components: dict[EmpiricalAccelerationComponents, list[EmpiricalAccelerationFunctionalShapes]], arc_start_times: list[float]) -> EstimatableParameterSettings:
    """Function for creating parameter settings for arc-wise empirical acceleration magnitudes.
    
    Function for creating parameter settings object for arc-wise empirical acceleration magnitudes (arc-wise version of :func:`~tudatpy.dynamics.parameters_setup.empirical_accelerations`).
    Using the empirical acceleration terms as estimatable parameters requires:
    
    * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.empirical` acceleration, which include constant (in RSW frame) terms
    
    When using this parameter, whenever an empirical acceleration is required at a time :math:`t`, the index math:`i` in the ``arc_initial_times`` ordered list is
    found for which :math:`t_{i}\\le t<t_{i+1}` (or, if :math:`t` is larger than the largest value in the list, :math:`i` is set to be last index of the list),
    and the parameter values representing empirical acceleration components in arc :math:`i` will be used.
    
    .. note:: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the radiation pressure coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose empirical acceleration model the estimatable parameter is associated.
    centralBody : str
        Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration
    acceleration_components : Dict[ EmpiricalAccelerationComponents, List[ EmpiricalAccelerationFunctionalShapes] ]
        Dictionary of components of the empirical acceleration which are to be estimated. There are two 'degrees of freedom' in these components: the direction of the acceleration (e.g. R, S or W direction) and the temporal signature (constant, sine of true anomaly or cosine of true anomaly). With this input, any subset may be selected. This parameter is a dictionary, with the key denoting the direction of the acceleration, and the value a list of the temporal signatures to estimate for this empirical acceleration direction.
    arc_initial_times : List[ float ]
        List of times at which the arcs over which the empirical accelerations are to be estimated will start.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.EmpiricalAccelerationEstimatableParameterSettings` class
        for the specified body's arc-wise empirical acceleration terms."""

def arcwise_radiation_pressure_coefficient(body: str, arc_initial_times: list[float]) -> EstimatableParameterSettings:
    """Function for creating parameter settings for arc-wise radiation pressure coefficients.
    
    Function for creating parameter settings object for arc-wise radiation pressure coefficients :math:`C_{r}` (arc-wise version of :func:`~tudatpy.dynamics.parameters_setup.radiation_pressure_coefficient`).
    Using the radiation pressure coefficient as an estimatable parameter requires:
    
    * A :func:`~tudatpy.dynamics.environment_setup.radiation_pressure.cannonball_radiation_target` target model to be defined for the body specified by the ``body`` parameter
    * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.radiation_pressure` acceleration
    
    When using this parameter, whenever :math:`C_{r}` is required at a time :math:`t`, the index math:`i` in the ``arc_initial_times`` ordered list is
    found for which :math:`t_{i}\\le t<t_{i+1}` (or, if :math:`t` is larger than the largest value in the list, :math:`i` is set to be last index of the list),
    and the parameter entry representing :math:`C_{r,i}` will be used.
    
    .. note:: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the radiation pressure coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.
    
    Parameters
    ----------
    body : str
        Name of the body, with whose radiation pressure acceleration model the estimatable parameter is associated.
    
    arc_initial_times : List[ float ]
        List of times at which the arcs over which the radiation pressure coefficient is to be estimated will start.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.ArcWiseRadiationPressureCoefficientEstimatableParameterSettings` class
        for arc-wise treatment of the specified body's radiation pressure coefficient."""

def arcwise_relative_observation_bias(link_ends: ..., observable_type: ..., arc_start_times: list[float], time_link_end: ...) -> EstimatableParameterSettings:
    """Function for creating parameter settings for arc-wise absolute observation bias.
    
    Function for creating parameter settings object for the arc-wise treatment of an observation's relative bias parameter.
    Using the arc-wise relative observation bias as estimatable parameter requires
    
    * The observation model (corresponding to the `link_ends` and `observable_type`) to include an arc-wise relative bias (:func:`~tudatpy.estimation.observable_models_setup.biases.arcwise_relative_bias`)
    
    .. note:: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the bias may, but need not, correspond to the arcs used for a multi-arc propagation.
    
    
    Parameters
    ----------
    link_ends : Dict[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, Tuple[str, str]
        Set of link ends that define the geometry of the biased observations.
    
    observable_type : ObservableType
        Observable type of the biased observations.
    arc_start_times : List[ float ]
        List of times at which the arcs over which the bias is to be estimated will start.
    time_link_end : LinkEndType
        The link end type (transmitter, receiver, etc.) at which the arc_start_times is evaluated.
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of the :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.ArcWiseConstantObservationBiasEstimatableParameterSettings`
        for the specified observation's arc-wise relative bias."""

def arcwise_time_bias(link_ends: dict[..., ...], observable_type: ..., arc_start_times: list[float], reference_link_end: ...) -> EstimatableParameterSettings:
    ...

def arcwise_time_drift_observation_bias(link_ends: dict[..., ...], observable_type: ..., arc_start_times: list[float], ref_epochs: list[float], time_link_end: ...) -> EstimatableParameterSettings:
    ...

def constant_drag_coefficient(body: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for constant drag coefficients.
    
    Function for creating parameter settings object for a constant drag coefficient parameter :math:`C_{D}`.
    Using the constant drag coefficient as an estimatable parameter requires:
    
    * A :func:`~tudatpy.dynamics.environment_setup.aerodynamic_coefficients.constant` aerodynamic interface to be defined for the body specified by the ``body`` parameter
    * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.aerodynamic` acceleration
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose drag acceleration model the estimatable parameter is associated.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's constant drag coefficient."""

def constant_empirical_acceleration_terms(body: str, centralBody: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for constant empirical acceleration terms.
    
    As :func:`~tudatpy.dynamics.parameters_setup.empirical_accelerations`, but only using the constant R, S and W components (no sine or cosine term estimation). This function is added as a function of convenience
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose empirical acceleration model the estimatable parameter is associated.
    
    centralBody : str
        Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.EmpiricalAccelerationEstimatableParameterSettings` class
        for the specified body's empirical acceleration terms."""

def constant_rotation_rate(body: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a body's constant rotation rate.
    
    Function for creating parameter settings object for a body's constant rotation rate parameter.
    Using the constant rotation rate as estimatable parameter requires:
    
    * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.simple` or :func:`~tudatpy.dynamics.environment_setup.rotation_model.simple_from_spice` rotation model specified by the ``body`` parameter
    * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose rotation model the estimatable parameter is associated.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's constant spin rate."""

def constant_time_bias(link_ends: dict[..., ...], observable_type: ..., reference_link_end: ...) -> EstimatableParameterSettings:
    ...

def core_factor(body: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a body's core factor.
    
    Function for creating parameter settings object for a body's core factor.
    Using the core factor as estimatable parameter requires
    
    * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
    * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose rotation model the estimatable parameter is associated.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's core factor."""

def create_parameter_set(parameter_settings: list[EstimatableParameterSettings], bodies: environment.SystemOfBodies, propagator_settings: propagator.PropagatorSettings=None, consider_parameters_names: list[EstimatableParameterSettings]=[]) -> ...:
    """Function for creating a consolidated parameter from the given estimatable parameter settings.
    
    Function for creating a consolidated parameter from the given estimatable parameter settings.
    The function checks for consistency between the parameter settings and the models contained in the simulation setup (given by the `bodies` & and `propagator_settings` parameters).
    
    
    Parameters
    ----------
    parameter_settings : list( :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` )
        List of objects that define the settings for the parameters that are to be created. Each entry in this list is typically created by a call to a function in the :ref:`parameter` module
    
    bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
        Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.
    
    propagator_settings : :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`
        Object containing the consolidated propagation settings of the simulation.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`
        Instance of :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet` class, consolidating all estimatable parameters and simulation models.
    
    Examples
    --------
    .. code-block:: python
    
        # Create bodies
        bodies = ...
        # Define parameters settings
        parameter_settings = ...
        # Create the parameters that will be estimated
        parameters_to_estimate = dynamics.parameters_setup.create_parameter_set(parameter_settings, bodies)
    
    This code snippet closely follows what is done in: `Full Estimation Example <https://github.com/tudat-team/tudatpy-examples/blob/master/estimation/full_estimation_example.ipynb>`_."""

def custom_analytical_partial(analytical_partial_function: typing.Callable[[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]], body_undergoing_acceleration: str, body_exerting_acceleration: str, acceleration_type: acceleration.AvailableAcceleration) -> CustomAccelerationPartialSettings:
    """No documentation found."""

def custom_numerical_partial(parameter_perturbation: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)], body_undergoing_acceleration: str, body_exerting_acceleration: str, acceleration_type: acceleration.AvailableAcceleration, environment_updates: dict[..., list[str]]={}) -> CustomAccelerationPartialSettings:
    """No documentation found."""

def custom_parameter(custom_id: str, parameter_size: int, get_parameter_function: typing.Callable[[], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], set_parameter_function: typing.Callable[[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], None]) -> EstimatableParameterSettings:
    """No documentation found."""

@typing.overload
def direct_tidal_dissipation_time_lag(body: str, deforming_body: str) -> EstimatableParameterSettings:
    """No documentation found."""

@typing.overload
def direct_tidal_dissipation_time_lag(body: str, deforming_body: list[str]) -> EstimatableParameterSettings:
    """No documentation found."""

def drag_component_scaling(body: str) -> EstimatableParameterSettings:
    """No documentation."""

def empirical_accelerations(body: str, centralBody: str, acceleration_components: dict[EmpiricalAccelerationComponents, list[EmpiricalAccelerationFunctionalShapes]]) -> EstimatableParameterSettings:
    """Function for creating parameter settings for empirical acceleration magnitudes.
    
    Function for creating parameter settings object for empirical acceleration magnitudes.
    Using the empirical acceleration terms as estimatable parameters requires:
    
    * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.empirical` acceleration, which include constant (in RSW frame) terms
    
    Any subset of the directions and functional shapes  can be estimated. The values in the parameter
    vector are ordered first by functional shape (constant, sine, cosine) and then by component (radial, normal, cross-track)
    For instance, if all nine coefficients are estimated, they will be ordered as :math:`\\mathbf{a}_{R,\\text{const.}},\\mathbf{a}_{R,\\text{sine}},\\mathbf{a}_{R,\\text{cosine}},\\mathbf{a}_{S,\\text{const.}},\\mathbf{a}_{S,\\text{sine}},\\mathbf{a}_{S,\\text{cosine}},\\mathbf{a}_{W,\\text{const.}},\\mathbf{a}_{W,\\text{sine}},\\mathbf{a}_{W,\\text{cosine}}`
    Any non-estimated components will be left to the values at which they were initialized.
    
    Parameters
    ----------
    body : str
        Name of the body, with whose empirical acceleration model the estimatable parameter is associated.
    
    centralBody : str
        Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration
    
    acceleration_components : dict[ EmpiricalAccelerationComponents, list[ EmpiricalAccelerationFunctionalShapes] ]
        Dictionary of components of the empirical acceleration which are to be estimated. There are two 'degrees of freedom' in these components: the direction of the acceleration (e.g. R, S or W direction) and the temporal signature (constant, sine of true anomaly or cosine of true anomaly). With this input, any subset may be selected. This parameter is a dictionary, with the key denoting the direction of the acceleration, and the value a list of the temporal signatures to estimate for this empirical acceleration direction.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.EmpiricalAccelerationEstimatableParameterSettings` class
        for the specified body's empirical acceleration terms."""

def free_core_nutation_rate(body: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a body's free core nutation rate.
    
    Function for creating parameter settings object for a body's free core nutation rate.
    Using the free core nutation rate as estimatable parameter requires
    
    * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
    * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose rotation model the estimatable parameter is associated.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's free core nutation rate."""

def full_empirical_acceleration_terms(body: str, centralBody: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for empirical acceleration magnitudes for all components.
    
    As :func:`~tudatpy.dynamics.parameters_setup.empirical_accelerations`, but using selecting all nine components. This function is added as a function of convenience
    
    Parameters
    ----------
    body : str
        Name of the body, with whose empirical acceleration model the estimatable parameter is associated.
    
    centralBody : str
        Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.EmpiricalAccelerationEstimatableParameterSettings` class
        for the specified body's empirical acceleration terms."""

def global_polynomial_clock_corrections(associated_body: str, associated_station: str, correction_powers: list[int]) -> EstimatableParameterSettings:
    ...

def gravitational_parameter(body: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a massive body's gravitational parameter.
    
    Function for creating parameter settings object for the gravitational parameter of massive bodies.
    Using the gravitational parameter as estimatable parameter requires:
    
    * The body specified by the ``body`` parameter to be endowed with a gravity field (see :ref:`gravity_field` module for options)
    * Any dynamical or observational model to depend on the gravitational parameter of the body specified by the ``body`` parameter
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose gravitational model the estimatable parameter is associated.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's gravitational parameter."""

def ground_station_position(body: str, ground_station_name: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for ground station position bias.
    
    Function for creating parameter settings object for a ground station's body-fixed Cartesian position.
    Using the ground station position bias as estimatable parameter requires:
    
    * At least one observation model to rely on the specified ground station
    
    
    Parameters
    ----------
    body : str
        Body name identifying the body, with which the ground station is associated.
    ground_station_name : str
        Name which identifies the position-biased ground station.
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified ground station's position bias."""

def initial_states(propagator_settings: propagator.PropagatorSettings, bodies: environment.SystemOfBodies, arc_initial_times: list[float]=[]) -> list[EstimatableParameterSettings]:
    """Function for creating parameter settings for initial state parameters.
    
    Function for creating a parameter settings object for initial state parameters.
    The function uses the propagator settings to determine which type of initial state parameter (single/multi/hybrid-arc; translational/rotational/... dynamics) is to be estimated,
    e.g. if a single-arc translational state propagator is defined, the function will automatically create the parameters for the associated initial state parameter
    
    .. note::
    
       This function return lists of parameter settings objects.
       This means that the return of this function cannot simply be added to the parameter settings objects of single parameters in a list creation statement.
       Instead, list concatenation is recommended. Please see the following example:
    
    .. code-block:: python
    
       # define single estimatable parameters
       single_parameter_1 = ...
       single_parameter_2 = ...
       ...
    
       # bad: list creation statement --> will result in nested list, undesired!
       list_of_all_parameters = [dynamics.parameters_setup.initial_states(...), single_parameter_1, single_parameter_2, ...]
    
       # better: list concatenation --> will result in simple list, desired!
       list_of_all_parameters = dynamics.parameters_setup.initial_states(...) + [single_parameter_1, single_parameter_2, ...]
    
    
    Parameters
    ----------
    propagator_settings : :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`
        Object containing the consolidated propagation settings of the simulation in the context of which the given model parameters are to be estimated.
    
    bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
        Object consolidating all bodies and environment models that constitute the physical environment.
    
    arc_initial_times : List[ float ] = []
        Initial times of arcs, only required if arc-wise propagation settings are passed via the `propagator_settings` object.
    
    Returns
    -------
    List[ :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` ]
        List of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` objects, one per component of each initial state in the simulation."""

@typing.overload
def inverse_tidal_quality_factor(body: str, deforming_body: str) -> EstimatableParameterSettings:
    """No documentation found."""

@typing.overload
def inverse_tidal_quality_factor(body: str, deforming_body: list[str]) -> EstimatableParameterSettings:
    """No documentation found."""

def lift_component_scaling(body: str) -> EstimatableParameterSettings:
    """No documentation."""

def mean_moment_of_inertia(body: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a body's mean moment of inertia.
    
    Function for creating parameter settings object for a body's mean moment of inertia.
    In most cases, the mean moment of inertia will not influence the dynamics/observation directly and sensitivity to this parameter will not be included. The dynamics/observation will be sensitive to this parameter if the rotational dynamics of a relevant body is estimated.
    Using the mean moment of inertia as estimatable parameter requires:
    
    * The estimation of an initial rotational state of the body specified by the ``body`` parameter
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose body model the estimatable parameter is associated.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's mean moment of inertia."""

def mode_coupled_k_love_numbers(deformed_body: str, love_number_indices: dict[tuple[int, int], list[tuple[int, int]]], deforming_bodies: list[str]) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a body's mode-coupled :math:`k_{lm}^{l'm'}` Love numbers.
    
    Function for creating parameter settings for a body's :math:`k_{lm}^{l'm'}` Love numbers (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.mode_coupled_solid_body_tide model` details). The estimation of
    the Love numbers can be limited to a subset of the bodies that raise a (mode-coupled) tide on the body undergoing tidal deformation.
    
    Using the :math:`k_{lm}` Love number as estimatable parameter requires:
    
    * A :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.mode_coupled_solid_body_tide` gravity field variation model in the ``deformed_body``.
    * Any dynamical model to depend on the gravity field of the body specified by the ``deformed_body`` parameter
    
    
    Parameters
    ----------
    deformed_body : str
        Name of the body that is undergoing tidal deformation
    love_number_per_degree : dict[tuple[int, int], list[int,int]]]
        Dictionary of Love number indices for each combination for forcing and response degree and order.
        The first tuple (key) is the forcing degree and order :math:`l,m`, the list of tuples (key) is the list of associated response degrees and orders :math:`l',m'`
        for which the Love numbers are to be estimated (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.mode_coupled_solid_body_tide` for mathematical definition))
    deforming_bodies : list[str]
        List of bodies that raise a tide on ``deformed_body`` for which the Love numbers defined by this setting is to be used. If the list is left empty, all tide-raising bodies will be used. By using this parameter, the values of :math:`k_{lm}` will be indentical for the tides raised by each body in this list once parameter values are reset, even if they were different upon environment initialization
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Object for the specified body's mode-coupled Love numbers :math:`k_{lm}^{l'm'}` for the tides raised by the specified bodies"""

def monomial_full_block_gravity_field_variation_amplitudes(body_name: str, power: int, minimum_degree: int, minimum_order: int, maximum_degree: int, maximum_order: int) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a body's polynomial gravity field amplitudes at a single power.
    
    Identical to :func:`~polynomial`, but for only a single power, and a full block of spherical harminic coefficient degrees :math:`l` and orders :math`m`
    For each degree :math:`l_{\\text{min}}\\le l \\le l_{\\text{max}}`, variations are estimated for all orders :math:`m_{\\text{min}}\\le m \\le \\left( \\text{min}(m_{\\text{max}},l) \\right)`
    
    Parameters
    ----------
    body_name : str
        Name of the body that is undergoing gravity field variation
    power: int
        Power :math:`i` for which to estimate polynomial gravity field variations
    minimum_degree: int
        Minimum degree :math:`l_{\\text{min}}` for which :math:`K_{i,\\bar{C}_{lm}}` and :math:`K_{i,\\bar{S}_{lm}}` are to be estimated (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)
    minimum_order: int
        Minimum order :math:`m_{\\text{min}}` for which :math:`K_{i,\\bar{C}_{lm}}` and :math:`K_{i,\\bar{S}_{lm}}` are to be estimated (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)
    maximum_degree: int
        Maximum degree :math:`l_{\\text{max}}` for which :math:`K_{i,\\bar{C}_{lm}}` and :math:`K_{i,\\bar{S}_{lm}}` are to be estimated (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)
    maximum_order: int
        Maximum degree :math:`m_{\\text{max}}` for which :math:`K_{i,\\bar{C}_{lm}}` and :math:`K_{i,\\bar{S}_{lm}}` are to be estimated (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Object for the specified body's polynomial gravity field variations"""

def monomial_gravity_field_variation_amplitudes(body_name: str, power: int, cosine_indices: list[tuple[int, int]], sine_indices: list[tuple[int, int]]) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a body's polynomial gravity field amplitudes at a single power.
    
    Identical to :func:`~polynomial`, but for only a single power.
    
    Parameters
    ----------
    body_name : str
        Name of the body that is undergoing gravity field variation
    power: int
        Power :math:`i` for which to estimate polynomial gravity field variations
    cosine_indices: list[int,int]
        List of combinations of degrees :math:`l` and orders :math:`m` for which to estimate :math:`K_{i,\\bar{C}_{lm}}` (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)
    sine_indices: list[int,int]
        List of combinations of degrees :math:`l` and orders :math:`m` for which to estimate :math:`K_{i,\\bar{S}_{lm}}` (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Object for the specified body's polynomial gravity field variations"""

def order_invariant_k_love_number(deformed_body: str, degree: int, deforming_bodies: list[str], use_complex_love_number: bool=0) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a body's :math:`k_{l}` Love number.
    
    Function for creating parameter settings for a body's :math:`k_{l}` Love number. When using this function,
    we assume (for the case of degree 2 Love number) that :math:`k_{20}=k_{21}=k_{22}`. The estimation of
    the Love number can be limited to a subset of the bodies that raise a tide on the body undergoing tidal deformation.
    
    Using the :math:`k_{l}` Love number as estimatable parameter requires:
    
    * A :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide` gravity field variation model in the ``deformed_body`` (or one the more complex ones such as :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k`). The parameter settings have to match the specifics of the variation model. For instance, if ``use_complex_love_number`` is set to true, the gravity field variation has to have been created using a complex Love number
    * Any dynamical model to depend on the gravity field of the body specified by the ``deformed_body`` parameter
    
    
    Parameters
    ----------
    deformed_body : str
        Name of the body that is undergoing tidal deformation
    degree : int
        Degree :math:`l` of the Love number :math:`k_{l}` that is to be estimated
    deforming_bodies : list[str]
        List of bodies that raise a tide on ``deformed_body`` for which the single Love number defined by this setting is to be used. If the list is left empty, all tide-raising bodies will be used. By using this parameter, the value of :math:`k_{l}` will be indentical for the tides raised by each body in this list once parameter values are reset, even if they were different upon environment initialization
    use_complex_love_number: bool
        Boolean defining whether the estimated Love number is real or imaginary
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Object for the specified body's Love number :math:`k_{l}` for the tides raised by the specified bodies"""

def order_varying_k_love_number(deformed_body: str, degree: int, orders: list[int], deforming_bodies: list[str], use_complex_love_number: bool) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a body's :math:`k_{lm}` Love numbers.
    
    Function for creating parameter settings for a body's :math:`k_{lm}` Love numbers. When using this function,
    we assume (for the case of degree 2 Love number) that :math:`k_{20}\\neq k_{21}\\neq k_{22}`. The estimation of
    the Love numbers can be limited to a subset of the bodies that raise a tide on the body undergoing tidal deformation.
    
    Using the :math:`k_{lm}` Love number as estimatable parameter requires:
    
    * A :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k` gravity field variation model in the ``deformed_body`` (or :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_complex_k` if ``use_complex_love_number`` is set to true)
    * Any dynamical model to depend on the gravity field of the body specified by the ``deformed_body`` parameter
    
    
    Parameters
    ----------
    deformed_body : str
        Name of the body that is undergoing tidal deformation
    degree : int
        Degree :math:`l` of the Love numbers :math:`k_{lm}` that are to be estimated
    degree : list[int]
        Orders :math:`m` of the Love numbers :math:`k_{lm}` that are to be estimated
    deforming_bodies : list[str]
        List of bodies that raise a tide on ``deformed_body`` for which the Love numbers defined by this setting is to be used. If the list is left empty, all tide-raising bodies will be used. By using this parameter, the values of :math:`k_{lm}` will be indentical for the tides raised by each body in this list once parameter values are reset, even if they were different upon environment initialization
    use_complex_love_number: bool
        Boolean defining whether the estimated Love number is real or imaginary
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Object for the specified body's Love numbers :math:`k_{lm}` for the tides raised by the specified bodies"""

def periodic_gravity_field_variation_amplitudes(body_name: str, cosine_indices_per_period: dict[int, list[tuple[int, int]]], sine_indices_per_period: dict[int, list[tuple[int, int]]]) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a body's polynomial gravity field variation amplitudes
    
    Function for creating parameter settings for a body's polynomial gravity field variation amplitudes
    :math:`A_{i,\\bar{C}_{lm}}`, :math:`B_{i,\\bar{C}_{lm}}`, :math:`A_{i,\\bar{S}_{lm}}` and :math:`B_{i,\\bar{S}_{lm}}`
    as defined in :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.single_period_periodic`.
    
    Using this settings as estimatable parameter requires:
    
    * A :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.periodic` (or :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.single_period_periodic`) gravity field variation model in the ``body_name``.
    * Any dynamical model to depend on the gravity field of the body specified by the ``deformed_body`` parameter
    
    When using this parameter, a subset of all the variation amplitudes defined in the gravity field variation model can be estimated.
    These are defined in the ``cosine_indices_per_period`` and ``sine_indices_per_period`` inputs
    
    Parameters
    ----------
    body_name : str
        Name of the body that is undergoing gravity field variation
    cosine_indices_per_period : dict[int, list[int,int]]
        Dictionary of frequency index :math:`i` (as keys; corresponding to frequency :math:`f_{i}`) with list of combinations of degrees :math:`l` and orders :math:`m` for which to estimate :math:`A_{i,\\bar{C}_{lm}}` and :math:`B_{i,\\bar{C}_{lm}}` as values (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.periodic` for mathematical definition)
    sine_indices_per_period : dict[int, list[int,int]]
        Dictionary of frequency index :math:`i` (as keys; corresponding to frequency :math:`f_{i}`) with list of combinations of degrees :math:`l` and orders :math:`m` for which to estimate :math:`A_{i,\\bar{S}_{lm}}` and :math:`B_{i,\\bar{S}_{lm}}` as values (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.periodic` for mathematical definition)
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Object for the specified body's periodic gravity field variations"""

def periodic_spin_variations(body: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a body's periodic spin variations.
    
    Function for creating parameter settings object for a body's periodic spin variation parameters.
    Using the mean moment of inertia as estimatable parameter requires:
    
    * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
    * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose rotation model the estimatable parameter is associated.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's periodic spin variations."""

def polar_motion_amplitudes(body: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a body's polar motion amplitudes.
    
    Function for creating parameter settings object for a body's polar motion amplitudes.
    Using the polar motion amplitudes as estimatable parameter requires
    
    * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
    * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose rotation model the estimatable parameter is associated.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's polar motion amplitudes."""

def polynomial_gravity_field_variation_amplitudes(body_name: str, cosine_indices_per_power: dict[int, list[tuple[int, int]]], sine_indices_per_power: dict[int, list[tuple[int, int]]]) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a body's polynomial gravity field amplitudes.
    
    Function for creating parameter settings for a body's polynomial gravity field amplitudes :math:`K_{i,\\bar{C}_{lm}}` and :math:`K_{i,\\bar{S}_{lm}}`,
    as defined in :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial`.
    
    Using this settings as estimatable parameter requires:
    
    * A :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` (or :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.single_power_polynomial`) gravity field variation model in the ``body_name``.
    * Any dynamical model to depend on the gravity field of the body specified by the ``deformed_body`` parameter
    
    When using this parameter, a subset of all the variation amplitudes defined in the gravity field variation model can be estimated. These are defined in the ``cosine_indices_per_power`` and ``sine_indices_per_power`` inputs
    
    Parameters
    ----------
    body_name : str
        Name of the body that is undergoing gravity field variation
    cosine_indices_per_power : dict[int, list[int,int]]
        Dictionary of powers :math:`i` (as keys) with list of combinations of degrees :math:`l` and orders :math:`m` for which to estimate :math:`K_{i,\\bar{C}_{lm}}` as values (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)
    sine_indices_per_power : dict[int, list[int,int]]
        Dictionary of powers :math:`i` (as keys) with list of combinations of degrees :math:`l` and orders :math:`m` for which to estimate :math:`K_{i,\\bar{S}_{lm}}` as values (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)
    
    Returns
    -------                                                                       -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Object for the specified body's polynomial gravity field variations"""

def ppn_parameter_beta() -> EstimatableParameterSettings:
    """Function for creating parameter settings for post-newtonian beta parameter.
    
    Function for creating parameter settings object for a global PPN :math:`\\beta` parameter.
    Using the post-newtonian gamma parameter as estimatable parameter requires at least one of the following:
    
    * An acceleration model depending on this parameter, such as :func:`~tudatpy.dynamics.propagation_setup.acceleration.relativistic_correction`
    * An observation model with a light-time correction depending on this parameter (none yet implemented)
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for a global post-newtonian :math:`\\beta` parameter."""

def ppn_parameter_gamma() -> EstimatableParameterSettings:
    """Function for creating parameter settings for post-newtonian gamma parameter.
    
    Function for creating parameter settings object for a global PPN :math:`\\gamma` parameter.
    Using the post-newtonian gamma parameter as estimatable parameter requires at least one of the following:
    
    * An acceleration model depending on this parameter, such as :func:`~tudatpy.dynamics.propagation_setup.acceleration.relativistic_correction`
    * An observation model with a light-time correction depending on this parameter, such as :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.first_order_relativistic_light_time_correction`
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for a global post-newtonian :math:`\\gamma` parameter."""

def quasi_impulsive_shots(body: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for quasi-impulsive shots.
    
    Function for creating parameter settings object for so-called 'quasi-impulsive shots', such as desaturation maneuvers.
    With this parameter, the total :math:`\\Delta \\mathbf{V}` vector of a set of such maneuvers can be estimated (see :func:`~tudatpy.dynamics.propagation_setup.acceleration.quasi_impulsive_shots_acceleration for mathematical details).
    Using the quasi-impulsive shots as an estimatable parameter requires:
    
    * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.quasi_impulsive_shots_acceleration` acceleration
    
    .. note:: this parameter considers *all* shots/maneuvers used in the above acceleration model, and estimates the value of the 'delta_v_values' input of this acceleration.
    
    Parameters
    ----------
    body : str
        Name of the body, with which the quasi-impulsive shot estimatable parameter is associated.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's quasi-impulsive shots"""

def radiation_pressure_coefficient(body: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for radiation pressure coefficients.
    
    Function for creating parameter settings object for a radiation pressure coefficient  :math:`C_{r}`.
    Using the radiation pressure coefficient as an estimatable parameter requires:
    
    * A :func:`~tudatpy.dynamics.environment_setup.radiation_pressure.cannonball` radiation pressure interface to be defined for the body specified by the ``body`` parameter
    * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.cannonball_radiation_pressure` acceleration
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose radiation pressure acceleration model the estimatable parameter is associated.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's radiation pressure coefficient."""

def radiation_pressure_target_direction_scaling(target_body: str, source_body: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a radiation pressure acceleration scaling factor in target direction.
    
    Function for creating parameter settings for scaling the radiation pressure acceleration component in the direction from the body
    undergoing the acceleration to the source model. When using this parameter, the radiation pressure :math:`\\mathbf{a}` is decomposed into
    a component :math:`\\mathbf{a}_{\\parallel}` and :math:`\\mathbf{a}_{\\perp}, such that :math:`\\mathbf{a}=\\mathbf{a}_{\\parallel}+\\mathbf{a}_{\\perp}`,
    where the parallel direction is computed as the component parallel with the vector from the center of mass of the source direction to the center of mass of the target direction.
    The radiation pressure model has parameters :math:`c_{\\parallel}` and :math:`c_{\\perp}` (nominally set to unity) that modify the acceleration as:
              
    .. math::
       \\mathbf{a}=c_{\\parallel}\\mathbf{a}_{\\parallel}+c_{\\perp}\\mathbf{a}_{\\perp}
    
    The present function creates settings for a parameter defining :math:`c_{\\parallel}` 
    
    Using this parameter requires:
    
    * The body specified by the ``target_body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.radiation_pressure` acceleration exerted by ``source_body``
              
    Parameters
    ----------
    target_body : str
        Name of the body on which radiation pressure is exerted
    source_body : str
        Name of the body exerting the radiation pressure
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` class dedining parallel radiation pressure scaling"""

def radiation_pressure_target_perpendicular_direction_scaling(target_body: str, source_body: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a radiation pressure acceleration scaling factor perpendicular to target direction.
    
    Function for creating parameter settings for scaling the radiation pressure acceleration component perpenedicular to the direction from the body
    undergoing the acceleration to the source model. The present function creates settings for a parameter defining :math:`c_{\\perp}`,
    see :func:`~radiation_pressure_target_direction_scaling`
    
    Using this parameter requires:
    
    * The body specified by the ``target_body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.radiation_pressure` acceleration exerted by ``source_body``
    
    Parameters
    ----------
    target_body : str
        Name of the body on which radiation pressure is exerted
    source_body : str
        Name of the body exerting the radiation pressure
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` class dedining parallel radiation pressure scaling"""

def reference_point_position(body: str, reference_point_name: str) -> EstimatableParameterSettings:
    """No documentation found."""

def relative_observation_bias(link_ends: ..., observable_type: ...) -> EstimatableParameterSettings:
    """Function for creating parameter settings for an relative observation bias.
    
    Function for creating parameter settings object for an observation's relative bias parameter.
    Using the relative observation bias as estimatable parameter requires
    
    * The observation model (corresponding to the `link_ends` and `observable_type`) to include a relative bias (:func:`~tudatpy.estimation.observable_models_setup.biases.relative_bias`)
    
    Parameters
    ----------
    link_ends : Dict[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, Tuple[str, str]
        Set of link ends that define the geometry of the biased observations.
    
    observable_type : ObservableType
        Observable type of the biased observations.
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of the :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.ConstantObservationBiasEstimatableParameterSettings`
        for the specified observation's arc-wise relative bias."""

def rotation_pole_position(body: str) -> EstimatableParameterSettings:
    """Function for creating parameter settings for a body's rotation pole position.
    
    Function for creating parameter settings object for a body's rotation pole position, parameterized by the constant pole rotation angles (:math:`\\alpha` and :math:`\\delta`).
    Using the rotation pole position as estimatable parameter requires:
    
    * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.simple` or :func:`~tudatpy.dynamics.environment_setup.rotation_model.simple_from_spice` rotation model specified by the ``body`` parameter
    * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose rotation model the estimatable parameter is associated.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's rotation pole position."""

def scaled_longitude_libration_amplitude(body_name: str) -> EstimatableParameterSettings:
    """No documentation found."""

def side_component_scaling(body: str) -> EstimatableParameterSettings:
    """No documentation."""

def spherical_harmonics_c_coefficients(body: str, minimum_degree: int, minimum_order: int, maximum_degree: int, maximum_order: int) -> EstimatableParameterSettings:
    """Function for creating parameter settings for the cosine coefficients of body's spherical harmonics gravitational model.
    
    Function for creating parameter settings object for the spherical harmonics cosine-coefficients (:math:`\\bar{C}_{lm}`) of a body with a spherical harmonic gravity field. Using this function, a 'full' set of spherical harmonic coefficients between an minimum/maximum degree/order are estimated. For instance, for minimum degree/order of 2/0, and maximum degree/order 4/4, all spherical harmonic cosine coefficients of degrees 2, 3 and 4 are estimated. If the maximum degree/order is set to 4/2, only coefficients with an order of 0, 1 and 2 are included. The entries in the parameter are sorted first by degree, and then by order (both in ascending order)
    Using the spherical harmonics cosine coefficients as estimatable parameter requires:
    
    * A :func:`~tudatpy.dynamics.environment_setup.gravity_field.spherical_harmonic` (or derived) gravity model to be defined for the body specified by the ``body`` parameter
    * Any dynamical or observational model to depend on the estimated cosine coefficients of the body specified by the ``body`` parameter. Typically, this dependency will be a :func:`~tudatpy.dynamics.propagation_setup.acceleration.spherical_harmonic` acceleration
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose gravitational model the estimatable parameters are associated.
    
    minimum_degree : int
        Minimum degree of c-coefficients to be included.
    minimum_order : int
        Minimum order of c-coefficients to be included.
    maximum_degree : int
        Maximum degree of c-coefficients to be included.
    maximum_order : int
        Maximum order of c-coefficients to be included.
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.SphericalHarmonicEstimatableParameterSettings` class
        for the applicable spherical harmonics c-coefficients of the specified body's gravitational model."""

def spherical_harmonics_c_coefficients_block(body: str, block_indices: list[tuple[int, int]]) -> EstimatableParameterSettings:
    """Function for creating parameter settings for the cosine coefficients of body's spherical harmonics gravitational model.
    
    As :class:`~tudatpy.dynamics.parameters_setup.spherical_harmonics_c_coefficients`, but with a manually defined set of coefficients.
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose gravitational model the estimatable parameters are associated.
    
    block_indices : List[ Tuple[int, int] ]
        List of block indices. The length of this list can be arbitrary, as long as the pairs are unique.
        For each pair, the first value is the degree and the second the order of the coefficient to be included.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.SphericalHarmonicEstimatableParameterSettings` class
        for the applicable spherical harmonics c-coefficients of the specified body's gravitational model."""

def spherical_harmonics_s_coefficients(body: str, minimum_degree: int, minimum_order: int, maximum_degree: int, maximum_order: int) -> EstimatableParameterSettings:
    """Function for creating parameter settings for the sine coefficients of body's spherical harmonics gravitational model.
    
    Function for creating parameter settings object for the spherical harmonics sine-coefficients (:math:`\\bar{S}_{lm}`) of a body with a spherical harmonic gravity field. Using this function, a 'full' set of spherical harmonic coefficients between an minimum/maximum degree/order are estimated. For instance, for minimum degree/order of 2/1 (there is no order 0 sine coefficient), and maximum degree/order 4/4, all spherical harmonic sine coefficients of degrees 2, 3 and 4 are estimated. If the maximum degree/order is set to 4/2, only coefficients with an order of 1 and 2 are included. The entries in the parameter are sorted first by degree, and then by order (both in ascending order)
    Using the spherical harmonics cosine coefficients as estimatable parameter requires:
    
    * A :func:`~tudatpy.dynamics.environment_setup.gravity_field.spherical_harmonic` (or derived) gravity model to be defined for the body specified by the ``body`` parameter
    * Any dynamical or observational model to depend on the estimated cosine coefficients of the body specified by the ``body`` parameter. Typically, this dependency will be a :func:`~tudatpy.dynamics.propagation_setup.acceleration.spherical_harmonic` acceleration
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose gravitational model the estimatable parameters are associated.
    
    minimum_degree : int
        Minimum degree of s-coefficients to be included.
    minimum_order : int
        Minimum order of s-coefficients to be included.
    maximum_degree : int
        Maximum degree of s-coefficients to be included.
    maximum_order : int
        Maximum order of s-coefficients to be included.
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.SphericalHarmonicEstimatableParameterSettings` class
        for the applicable spherical harmonics s-coefficients of the specified body's gravitational model."""

def spherical_harmonics_s_coefficients_block(body: str, block_indices: list[tuple[int, int]]) -> EstimatableParameterSettings:
    """Function for creating parameter settings for the sine coefficients of body's spherical harmonics gravitational model.
    
    As :class:`~tudatpy.dynamics.parameters_setup.spherical_harmonics_s_coefficients`, but with a manually defined set of coefficients.
    
    
    Parameters
    ----------
    body : str
        Name of the body, with whose gravitational model the estimatable parameters are associated.
    
    block_indices : List[ Tuple[int, int] ]
        List of block indices. The length of this list can be arbitrary, as long as the pairs are unique.
        For each pair, the first value is the degree and the second the order of the coefficient to be included.
    
    Returns
    -------
    :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
        Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.SphericalHarmonicEstimatableParameterSettings` class
        for the applicable spherical harmonics s-coefficients of the specified body's gravitational model."""

def time_drift_observation_bias(link_ends: dict[..., ...], observable_type: ..., ref_epoch: float, time_link_end: ...) -> EstimatableParameterSettings:
    ...

def yarkovsky_parameter(body_name: str, central_body_name: str='Sun') -> EstimatableParameterSettings:
    """No documentation found."""
across_track_empirical_acceleration_component: EmpiricalAccelerationComponents
along_track_empirical_acceleration_component: EmpiricalAccelerationComponents
arc_wise_constant_drag_coefficient_type: EstimatableParameterTypes
arc_wise_empirical_acceleration_coefficients_type: EstimatableParameterTypes
arc_wise_initial_body_state_type: EstimatableParameterTypes
arc_wise_polynomial_clock_corrections_type: EstimatableParameterTypes
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
global_polynomial_clock_corrections_type: EstimatableParameterTypes
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