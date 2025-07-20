import numpy
import pybind11_stubgen.typing_ext
from ....astro import time_representation
import typing
__all__ = ['AdamsBashforthMoultonSettings', 'AvailableIntegrators', 'BulirschStoerIntegratorSettings', 'CoefficientSets', 'ExtrapolationMethodStepSequences', 'IntegratorSettings', 'IntegratorStepSizeControlSettings', 'IntegratorStepSizeValidationSettings', 'MinimumIntegrationTimeStepHandling', 'OrderToIntegrate', 'RungeKuttaFixedStepSizeSettings', 'RungeKuttaVariableStepSizeBaseSettings', 'RungeKuttaVariableStepSizeSettingsScalarTolerances', 'RungeKuttaVariableStepSizeSettingsVectorTolerances', 'SSPRK3', 'adams_bashforth_moulton', 'adams_bashforth_moulton_fixed_order', 'adams_bashforth_moulton_fixed_step', 'adams_bashforth_moulton_fixed_step_fixed_order', 'adams_bashforth_moulton_type', 'bulirsch_stoer', 'bulirsch_stoer_fixed_step', 'bulirsch_stoer_sequence', 'bulirsch_stoer_type', 'bulirsch_stoer_variable_step', 'deufelhard_sequence', 'euler', 'euler_forward', 'explicit_mid_point', 'explicit_trapezoid_rule', 'heun_euler', 'higher', 'lower', 'print_butcher_tableau', 'ralston', 'ralston_3', 'ralston_4', 'rk_3', 'rk_4', 'rkdp_87', 'rkf_108', 'rkf_12', 'rkf_1210', 'rkf_1412', 'rkf_45', 'rkf_56', 'rkf_78', 'rkf_89', 'rkv_89', 'runge_kutta_4', 'runge_kutta_fixed_step', 'runge_kutta_fixed_step_size', 'runge_kutta_fixed_step_size_type', 'runge_kutta_variable_step', 'runge_kutta_variable_step_size', 'runge_kutta_variable_step_size_type', 'runge_kutta_variable_step_size_vector_tolerances', 'set_to_minimum_step_every_time_warning', 'set_to_minimum_step_silently', 'set_to_minimum_step_single_warning', 'standard_cartesian_state_element_blocks', 'standard_rotational_state_element_blocks', 'step_size_control_blockwise_matrix_tolerance', 'step_size_control_blockwise_scalar_tolerance', 'step_size_control_custom_blockwise_matrix_tolerance', 'step_size_control_custom_blockwise_scalar_tolerance', 'step_size_control_elementwise_matrix_tolerance', 'step_size_control_elementwise_scalar_tolerance', 'step_size_validation', 'three_eight_rule_rk_4', 'throw_exception_below_minimum']

class AdamsBashforthMoultonSettings(IntegratorSettings):
    """`IntegratorSettings`-derived class to define settings for Adams-Bashforth-Moulton integrator settings."""

class AvailableIntegrators:
    """Enumeration of integrators available with tudat.
    
    
    
    
    
          
    
    Members:
    
      runge_kutta_fixed_step_size_type
    
      runge_kutta_variable_step_size_type
    
      bulirsch_stoer_type
    
      adams_bashforth_moulton_type"""
    __members__: typing.ClassVar[dict[str, AvailableIntegrators]]
    adams_bashforth_moulton_type: typing.ClassVar[AvailableIntegrators]
    bulirsch_stoer_type: typing.ClassVar[AvailableIntegrators]
    runge_kutta_fixed_step_size_type: typing.ClassVar[AvailableIntegrators]
    runge_kutta_variable_step_size_type: typing.ClassVar[AvailableIntegrators]

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

class BulirschStoerIntegratorSettings(IntegratorSettings):
    """`IntegratorSettings`-derived class to define settings for Bulirsch-Stoer integrator settings."""

class CoefficientSets:
    """Coefficient sets for Runge-Kutta-type integrators.
    
    Coefficient sets for Runge-Kutta-type integrators. The coefficients are defined
    in a Butcher Tableau, with an coefficient set yielding an x(y) method yielding an integrator
    with global truncation error of :math:`O(\\Delta t^{x})`. Some of these coefficients also contain an embedded integrator of :math:`O(\\Delta t^{y})`
    for step size control.
    
          
    
    Members:
    
      euler_forward : 
    
    Coefficients for the classic forward Euler method
    
    
    
      rk_4 : 
    
    Coefficients for the original Runge-Kutta method of order 4
    
    
    
      explicit_mid_point : 
    
    Coefficients for the explicit midpoint method
    
    
    
      explicit_trapezoid_rule : 
    
    Coefficients for the explicit trapezoid rule, also called Heun's method or improved Euler's method
    
    
    
      ralston : 
    
    Coefficients for Ralston's method
    
    
    
      rk_3 : 
    
    Coefficients for the Runge-Kutta method of order 3
    
    
    
      ralston_3 : 
    
    Coefficients for Ralston's third-order method
    
    
    
      SSPRK3 : 
    
    Coefficients for the Strong Stability Preserving Runge-Kutta third-order method
    
    
    
      ralston_4 : 
    
    Coefficients for Ralston's fourth-order method
    
    
    
      three_eight_rule_rk_4 : 
    
    Coefficients for the classic Runge Kutta 3/8-rule fourth-order method
    
    
    
      heun_euler : 
    
    Coefficients for the Heun's method of order 2 with an embedded Euler method of order 1
    
    
    
      rkf_12 : 
    
    Coefficients for the Runge-Kutta-Fehlberg method of order 2 with an embedded 1st order
    
    
    
      rkf_45 : 
    
    Coefficients for the Runge-Kutta-Fehlberg method of order 5 with an embedded 4th order
    
    
    
      rkf_56 : 
    
    Coefficients for the Runge-Kutta-Fehlberg method of order 6 with an embedded 5th order
    
    
    
      rkf_78 : 
    
    Coefficients for the Runge-Kutta-Fehlberg method of order 8 with an embedded 7th order
    
    
    
      rkdp_87 : 
    
    Coefficients for the Dormand-Prince method of order 7 with an embedded 8th order
    
    
    
      rkf_89 : 
    
    Coefficients for the Runge-Kutta-Fehlberg method of order 9 with an embedded 8th order
    
    
    
      rkv_89 : 
    
    Coefficients for the Runge-Kutta-Verner method of order 9 with an embedded 8th order
    
    
    
      rkf_108 : 
    
    Coefficients for the Runge-Kutta-Feagin method of order 8 with an embedded 10th order
    
    
    
      rkf_1210 : 
    
    Coefficients for the Runge-Kutta-Feagin method of order 10 with an embedded 12ve order
    
    
    
      rkf_1412 : 
    
    Coefficients for the Runge-Kutta-Feagin method of order 12 with an embedded 14th order"""
    SSPRK3: typing.ClassVar[CoefficientSets]
    __members__: typing.ClassVar[dict[str, CoefficientSets]]
    euler_forward: typing.ClassVar[CoefficientSets]
    explicit_mid_point: typing.ClassVar[CoefficientSets]
    explicit_trapezoid_rule: typing.ClassVar[CoefficientSets]
    heun_euler: typing.ClassVar[CoefficientSets]
    ralston: typing.ClassVar[CoefficientSets]
    ralston_3: typing.ClassVar[CoefficientSets]
    ralston_4: typing.ClassVar[CoefficientSets]
    rk_3: typing.ClassVar[CoefficientSets]
    rk_4: typing.ClassVar[CoefficientSets]
    rkdp_87: typing.ClassVar[CoefficientSets]
    rkf_108: typing.ClassVar[CoefficientSets]
    rkf_12: typing.ClassVar[CoefficientSets]
    rkf_1210: typing.ClassVar[CoefficientSets]
    rkf_1412: typing.ClassVar[CoefficientSets]
    rkf_45: typing.ClassVar[CoefficientSets]
    rkf_56: typing.ClassVar[CoefficientSets]
    rkf_78: typing.ClassVar[CoefficientSets]
    rkf_89: typing.ClassVar[CoefficientSets]
    rkv_89: typing.ClassVar[CoefficientSets]
    three_eight_rule_rk_4: typing.ClassVar[CoefficientSets]

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

class ExtrapolationMethodStepSequences:
    """Enumeration of available extrapolation method substep sequences, with :math:`n_{j}` defining the number of substeps in iteration :math:`j`.
    
    
    
    
    
          
    
    Members:
    
      bulirsch_stoer_sequence : 
    
    Sequence for which :math:`n_{j}=2n_{j-2}` (2, 4, 6, 8, 12, 16, 24, ....)
    
    
    
      deufelhard_sequence : 
    
    Sequence for which :math:`n_{j}=2(j+1)` (2, 4, 6, 8, 10, 12, 14, ....)"""
    __members__: typing.ClassVar[dict[str, ExtrapolationMethodStepSequences]]
    bulirsch_stoer_sequence: typing.ClassVar[ExtrapolationMethodStepSequences]
    deufelhard_sequence: typing.ClassVar[ExtrapolationMethodStepSequences]

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

class IntegratorSettings:
    """Functional base class to define settings for integrators.
    
    Class to define settings for numerical integrators, for instance for use in numerical integration of equations of motion/
    variational equations. This class can be used for simple integrators such as fixed step RK and Euler. Integrators that
    require more settings to define have their own derived class."""
    initial_time: time_representation.Time

class IntegratorStepSizeControlSettings:
    """Base class to define settings for step-size control algorithm.
    
    Base class to define settings for step-size control algorithm, typically created by one of the functions provided in this module"""
    maximum_step_decrease: float
    minimum_step_decrease: float
    safety_factor: float

class IntegratorStepSizeValidationSettings:
    """Base class to define settings for step-size validation algorithm.
    
    Base class to define settings for step-size validation algorithm, typically created by one of the functions provided in this module"""
    maximum_step: float
    minimum_step: float
    minimum_step_handling: MinimumIntegrationTimeStepHandling

class MinimumIntegrationTimeStepHandling:
    '''Enumeration defining possible behaviours when :math:`\\Delta t_{rec}<\\Delta t_{\\min}`. in step-size control (e.g. recommended time step is smaller than minimum time step)
    
    
    
    
    
          
    
    Members:
    
      throw_exception_below_minimum : 
    
    The propagation is terminated and a :class:`tudatpy.exceptions.MinimumStepSizeViolatedError` is thrown.
    
    
    
      set_to_minimum_step_silently : 
    
    The final time step is set to :math:`\\Delta t=\\Delta t_{\\min}`, violating requirements of step-size control algorithm, without any message to user"
    
    
    
      set_to_minimum_step_single_warning : 
    
    The final time step is set to :math:`\\Delta t=\\Delta t_{\\min}`, violating requirements of step-size control algorithm, a warning is printed to the terminal the first time this happens during a propagation"
    
    
    
      set_to_minimum_step_every_time_warning : 
    
    The final time step is set to :math:`\\Delta t=\\Delta t_{\\min}`, violating requirements of step-size control algorithm, a warning is printed to the terminal every time this happens during a propagation"'''
    __members__: typing.ClassVar[dict[str, MinimumIntegrationTimeStepHandling]]
    set_to_minimum_step_every_time_warning: typing.ClassVar[MinimumIntegrationTimeStepHandling]
    set_to_minimum_step_silently: typing.ClassVar[MinimumIntegrationTimeStepHandling]
    set_to_minimum_step_single_warning: typing.ClassVar[MinimumIntegrationTimeStepHandling]
    throw_exception_below_minimum: typing.ClassVar[MinimumIntegrationTimeStepHandling]

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

class OrderToIntegrate:
    """Enumeration defining Which integrator order needs to be integrated, only used for coefficient sets with an embedded order.
    
    
    
    
    
          
    
    Members:
    
      lower : 
    
    For a method of order :math:`p`, with embedded method of order :math:`q`, the step is taken using the method with order :math:`\\min(p,q)`
    
    
    
      higher : 
    
    For a method of order :math:`p`, with embedded method of order :math:`q`, the step is taken using the method with order :math:`\\max(p,q)`"""
    __members__: typing.ClassVar[dict[str, OrderToIntegrate]]
    higher: typing.ClassVar[OrderToIntegrate]
    lower: typing.ClassVar[OrderToIntegrate]

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

class RungeKuttaFixedStepSizeSettings(IntegratorSettings):
    """`IntegratorSettings`-derived class to define settings for Runge Kutta integrators with a fixed step size"""

class RungeKuttaVariableStepSizeBaseSettings(IntegratorSettings):
    """No documentation found."""

class RungeKuttaVariableStepSizeSettingsScalarTolerances(RungeKuttaVariableStepSizeBaseSettings):
    """No documentation found."""

class RungeKuttaVariableStepSizeSettingsVectorTolerances(RungeKuttaVariableStepSizeBaseSettings):
    """No documentation found."""

@typing.overload
def adams_bashforth_moulton(initial_time_step: time_representation.Time, minimum_step_size: time_representation.Time, maximum_step_size: time_representation.Time, relative_error_tolerance: float=1e-12, absolute_error_tolerance: float=1e-12, minimum_order: int=6, maximum_order: int=11, assess_termination_on_minor_steps: bool=False, bandwidth: time_representation.Time=200.0) -> IntegratorSettings:
    """Creates the settings for the Adams-Bashforth-Moulton integrator.
    
    Function to create settings for the Adams-Bashforth-Moulton multistep integrator.
    For this integrator, the step size and order are both according to a control algorithm
    similar to :func:`~tudatpy.dynamics.propagation_setup.integrator.step_size_control_elementwise_scalar_tolerance`.
    The integrator is initialized using an RKF7(8) integrator.
    
    NOTE: this integrator's step-size and order control algorithm work in a method that is overly simplistic,
    when increasing/decreasing the order, existing function evaluations are re-used, without any recomputations.
    Similarly, when halving or doubling the time-step, the existing interpolating polynomial is evaluated at the relevant points.
    This can lead to unwanted behaviour, where the time-step reduces to unrealistically low values. It is strongly
    recommended that a reasonable minimum step is provided to this function, to partially mitigate this behaviour.
    
    Parameters
    ----------
    initial_time_step : float
        Initial time step to be used.
    minimum_step_size : float
        Minimum time step to be used during the integration.
    maximum_step_size : float
        Maximum time step to be used during the integration.
    relative_error_tolerance : float, default=1.0E-12
        Relative tolerance to adjust the time step.
    absolute_error_tolerance : float, default=1.0E-12
        Relative tolerance to adjust the time step.
    minimum_order
        Minimum order of the integrator.
    maximum_order
        Maximum order of the integrator.
    assess_termination_on_minor_steps : bool, default=false
        Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
    bandwidth : float, default=200.0
        Maximum error factor for doubling the step size.
    Returns
    -------
    IntegratorSettings
        Object containing settings for the integrator."""

@typing.overload
def adams_bashforth_moulton(initial_time: time_representation.Time, initial_time_step: time_representation.Time, minimum_step_size: time_representation.Time, maximum_step_size: time_representation.Time, relative_error_tolerance: float=1e-12, absolute_error_tolerance: float=1e-12, minimum_order: int=6, maximum_order: int=11, assess_termination_on_minor_steps: bool=False, bandwidth: time_representation.Time=200.0) -> IntegratorSettings:
    ...

def adams_bashforth_moulton_fixed_order(initial_time_step: time_representation.Time, minimum_step_size: time_representation.Time, maximum_step_size: time_representation.Time, relative_error_tolerance: float=1e-12, absolute_error_tolerance: float=1e-12, order: int=6, assess_termination_on_minor_steps: bool=False, bandwidth: time_representation.Time=200.0) -> IntegratorSettings:
    """Creates the settings for the Adams-Bashforth-Moulton integrator of fixed order.
    
    Same as :func:`~tudatpy.dynamics.propagation_setup.integrator.adams_bashforth_moulton`, but
    with fixed order and variable step
    
    Parameters
    ----------
    initial_time_step : float
        Initial time step to be used.
    minimum_step_size : float
        Minimum time step to be used during the integration.
    maximum_step_size : float
        Maximum time step to be used during the integration.
    relative_error_tolerance : float, default=1.0E-12
        Relative tolerance to adjust the time step.
    absolute_error_tolerance : float, default=1.0E-12
        Relative tolerance to adjust the time step.
    order
        Order of the integrator.
    assess_termination_on_minor_steps : bool, default=false
        Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
    bandwidth : float, default=200.0
        Maximum error factor for doubling the step size.
    Returns
    -------
    IntegratorSettings
        Object containing settings for the integrator."""

def adams_bashforth_moulton_fixed_step(time_step: time_representation.Time, relative_error_tolerance: float=1e-12, absolute_error_tolerance: float=1e-12, minimum_order: int=6, maximum_order: int=11, assess_termination_on_minor_steps: bool=False, bandwidth: time_representation.Time=200.0) -> IntegratorSettings:
    """Creates the settings for the Adams-Bashforth-Moulton fixed-step integrator.
    
    Same as :func:`~tudatpy.dynamics.propagation_setup.integrator.adams_bashforth_moulton`, but
    with fixed step and variable order
    
    
    Parameters
    ----------
    time_step : float
        Initial time step to be used.
    relative_error_tolerance : float, default=1.0E-12
        Relative tolerance to adjust the time step.
    absolute_error_tolerance : float, default=1.0E-12
        Relative tolerance to adjust the time step.
    minimum_order
        Minimum order of the integrator.
    maximum_order
        Maximum order of the integrator.
    assess_termination_on_minor_steps : bool, default=false
        Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
    bandwidth : float, default=200.0
        Maximum error factor for doubling the step size.
    Returns
    -------
    IntegratorSettings
        Object containing settings for the integrator."""

def adams_bashforth_moulton_fixed_step_fixed_order(time_step: time_representation.Time, order: int=6, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    """Creates the settings for the Adams-Bashforth-Moulton fixed-step, fixed-order integrator.
    
    Same as :func:`~tudatpy.dynamics.propagation_setup.integrator.adams_bashforth_moulton`, but
    with fixed step and fixed order
    
    
    Parameters
    ----------
    time_step : float
        Initial time step to be used.
    order
        Order of the integrator.
    assess_termination_on_minor_steps : bool, default=false
        Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
    Returns
    -------
    IntegratorSettings
        Object containing settings for the integrator."""

@typing.overload
def bulirsch_stoer(initial_time: time_representation.Time, initial_time_step: time_representation.Time, extrapolation_sequence: ExtrapolationMethodStepSequences, maximum_number_of_steps: int, minimum_step_size: time_representation.Time, maximum_step_size: time_representation.Time, relative_error_tolerance: float=1e-12, absolute_error_tolerance: float=1e-12, assess_termination_on_minor_steps: bool=False, safety_factor: time_representation.Time=0.7, maximum_factor_increase: time_representation.Time=10.0, minimum_factor_increase: time_representation.Time=0.1) -> IntegratorSettings:
    ...

@typing.overload
def bulirsch_stoer(initial_time_step: time_representation.Time, extrapolation_sequence: ExtrapolationMethodStepSequences, maximum_number_of_steps: int, minimum_step_size: time_representation.Time, maximum_step_size: time_representation.Time, relative_error_tolerance: float=1e-12, absolute_error_tolerance: float=1e-12, assess_termination_on_minor_steps: bool=False, safety_factor: time_representation.Time=0.7, maximum_factor_increase: time_representation.Time=10.0, minimum_factor_increase: time_representation.Time=0.1) -> IntegratorSettings:
    """Creates the settings for the Bulirsch-Stoer integrator.
    
    
    NOTE: THIS FUNCTION IS DEPRECATED, IT IS RECOMMENDED TO USE THE NEW :func:`~tudatpy.dynamics.propagation_setup.integrator.bulirsch_stoer_variable_step` INTERFACE INSTEAD
    
    Function to create settings for the Bulirsch-Stoer integrator.
    For this integrator, the step size is varied based on the tolerances and safety factor provided.
    The tolerance is composed of an absolute and a relative part.
    Different extrapolation sequences can be used (see the `ExtrapolationMethodStepSequences` enum).
    
    
    Parameters
    ----------
    initial_time_step : float
        Initial time step to be used.
    extrapolation_sequence : ExtrapolationMethodStepSequences
        Extrapolation sequence to be used in the integration.
    maximum_number_of_steps : int
        Number of entries in the sequence (e.g., number of integrations used for a single extrapolation).
    minimum_step_size : float
        Minimum time step to be used during the integration.
    maximum_step_size : float
        Maximum time step to be used during the integration.
    relative_error_tolerance : float, default=1.0E-12
        Relative tolerance to adjust the time step.
    absolute_error_tolerance : float, default=1.0E-12
        Relative tolerance to adjust the time step.
    assess_termination_on_minor_steps : bool, default=false
        Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
    safety_factor : float, default=0.7
        Safety factor used in the step size control.
    maximum_factor_increase : float, default=10.0
        Maximum increase between consecutive time steps, expressed as the factor between new and old step size.
    minimum_factor_increase : float, default=0.1
        Minimum increase between consecutive time steps, expressed as the factor between new and old step size.
    Returns
    -------
    BulirschStoerIntegratorSettings
        BulirschStoerIntegratorSettings object."""

def bulirsch_stoer_fixed_step(time_step: time_representation.Time, extrapolation_sequence: ExtrapolationMethodStepSequences, maximum_number_of_steps: int, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    """Creates the settings for the fixed time-step Bulirsch-Stoer integrator.
    
    Function to create settings for the fixed time-step Bulirsch-Stoer integrator. The
    underlying method is the same as :func:`~tudatpy.dynamics.propagation_setup.integrator.bulirsch_stoer_variable_step`,
    but using a fixed, user-defined, time step.
    
    
    Parameters
    ----------
    time_step : float
        Time step to be used.
    extrapolation_sequence : ExtrapolationMethodStepSequences
        Extrapolation sequence to be used for the integration (defining the number of substeps in iteration :math:`i`).
    maximum_number_of_steps : int
        Number of entries from the sequence to be used (e.g., total number of iterations used for a single extrapolation and time step).
    assess_termination_on_minor_steps : bool, default=false
        Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
    Returns
    -------
    IntegratorSettings
        Object containing settings for the integrator.
    
    
    
    
    
    Examples
    --------
    In this example, settings for the Bulirsch-Stoer integrator with 600 second time step are created, using the typical
    sequence, using 6 iterations of the same step. By using the Bulirsch-Stoer sequence, this means that the same step is done
    using 2, 4, 6, 8, 12 and 16 substeps
    
    .. code-block:: python
    
      # Create BS settings
      integrator_settings = integrator.bulirsch_stoer_fixed_step(
          time_step = 300.0,
          extrapolation_sequence = integrator.bulirsch_stoer_sequence,
          maximum_number_of_steps = 6 )"""

def bulirsch_stoer_variable_step(initial_time_step: time_representation.Time, extrapolation_sequence: ExtrapolationMethodStepSequences, maximum_number_of_steps: int, step_size_control_settings: IntegratorStepSizeControlSettings, step_size_validation_settings: IntegratorStepSizeValidationSettings, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    """Creates the settings for the variable time-step Bulirsch-Stoer integrator.
    
    Function to create settings for the variable time-step Bulirsch-Stoer integrator. This integrator
    works by performing the same (typically very large) step multiple times, using an ever increasing number of substeps.
    Each substep is performed using the modified midpoint method. The successive integrations from :math:`t_{i}` to :math:`t_{i}+\\Delta t`
    are (in principle) done using ever increasing accuracy, as the size of the substep decreases. This integrator works
    by extrapolating the behaviour to a substep length of 0 (e.g. an infinite number of substeps), at which the solution should be perfect.
    The number of substeps on the :math:`i^{t}` iteration are done using the number of substeps defined by  entry :math:`i` of the
    ``extrapolation_sequence`` input. The number of iterations for a single step is defined by the ``maximum_number_of_steps`` entry.
    For instance, using the ``bulirsch_stoer_sequence`` sequence, and 5 iterations, the same step is done using 2, 4, 6, 8 and 12 substeps,
    and the results are then extrapolated to an infinite number of steps. Different extrapolation sequences can be used (see the `ExtrapolationMethodStepSequences` enum).
    
    The step-size control algorithm is defined by a :class:`~tudatpy.dynamics.propagation_setup.integrator.IntegratorStepSizeControlSettings` and
    :class:`~tudatpy.dynamics.propagation_setup.integrator.IntegratorStepSizeValidationSettings` object, created using one of the functions
    listed above. The time step control uses the result from the final, and second to final iteration to generate an error estimate of the current step.
    
    
    Parameters
    ----------
    time_step : float
        Initial time step to be used.
    extrapolation_sequence : ExtrapolationMethodStepSequences
        Extrapolation sequence to be used for the integration (defining the number of substeps in iteration :math:`i`).
    maximum_number_of_steps : int
        Number of entries from the sequence to be used (e.g., total number of iterations used for a single extrapolation and time step).
    step_size_control_settings : IntegratorStepSizeControlSettings
        Object used to control the step size, by computing a new step size :math:`\\Delta t_{rec.}`, from the embedded Runge-Kutta integrator pair,
        and recommending whether the steps is to be accepted, or recomputed with a different time step.
    
    step_size_validation_settings : IntegratorStepSizeValidationSettings
        Object used to validate whether the :math:`\\Delta t_{rec.}` provided by model defined by the ``step_size_control_settings`` meets with user-defined
        criteria (minimum, maximum values, etc.)
    
    assess_termination_on_minor_steps : bool, default=false
        Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
    Returns
    -------
    IntegratorSettings
        Object containing settings for the integrator.
    
    
    
    
    
    Examples
    --------
    In this example, settings for the Bulirsch-Stoer integrator with 600 second initial time step are created, using the typical
    sequence, using 6 iterations of the same step. By using the Bulirsch-Stoer sequence, this means that the same step is done
    using 2, 4, 6, 8, 12 and 16 substeps. The same tolerances (:math:`10^{-10}`)
    applied element-wise on the propagated state. The minimum and maximum time steps are set to 0.1 and 10000 seconds,
    and the initial step is set to 600 seconds. All other inputs are left on their defaults
    
    .. code-block:: python
    
      # Create BS settings
      control_settings = integrator.step_size_control_elementwise_scalar_tolerance( 1.0E-10, 1.0E-10 )
      validation_settings = integrator.step_size_validation( 0.1, 10000.0 )
      integrator_settings = integrator.bulirsch_stoer_variable_step(
          initial_time_step = 600.0,
          extrapolation_sequence = integrator.bulirsch_stoer_sequence,
          maximum_number_of_steps = 6,
          step_size_control_settings = control_settings,
          step_size_validation_settings = validation_settings )"""

@typing.overload
def euler(initial_time: time_representation.Time, initial_time_step: time_representation.Time, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

@typing.overload
def euler(initial_time_step: time_representation.Time, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

def print_butcher_tableau(coefficient_set: CoefficientSets) -> None:
    """Print the Butcher tableau of a given coefficient set.
    
    
    Parameters
    ----------
    coefficient_set : CoefficientSets
        Coefficient set of which the Butcher tableau will be printed."""

@typing.overload
def runge_kutta_4(initial_time: time_representation.Time, initial_time_step: time_representation.Time, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

@typing.overload
def runge_kutta_4(initial_time_step: time_representation.Time, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

def runge_kutta_fixed_step(time_step: time_representation.Time, coefficient_set: CoefficientSets, order_to_use: OrderToIntegrate=..., assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    """Creates the settings for the Runge-Kutta fixed step size integrator.
    
    Function to create settings for the Runge-Kutta integrator with a constant step size.
    Different coefficient sets (Butcher's tableau) can be used (see the `CoefficientSets` enum).
    
    
    Parameters
    ----------
    time_step : float
        Initial time step to be used.
    coefficient_set : CoefficientSets
        Coefficient set (Butcher's tableau) to be used in the integration.
    order_to_use : OrderToIntegrate, default=OrderToIntegrate.lower
        If the coefficient set is supposed to be for variable step sizes (with an embedded method of a different order),
        this parameter can be used to set the order that will be used.
    
    assess_termination_on_minor_steps : bool, default=false
        Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
        integrator (true) or only at the end of each integration step (false).
    
    Returns
    -------
    IntegratorSettings
        Object containing settings for the integrator.
    
    
    
    
    
    Examples
    --------
    In this example, settings for the classical RK4 integrator with 30 second time step are created
    
    .. code-block:: python
    
      # Create RK4 settings
      integrator_settings = integrator.runge_kutta_fixed_step(
          time_step = 30.0,
          coefficient_set = integrator.CoefficientSets.rk_4 )
    
    In this example, settings for fixed-step integration using the higher-order (8th-order) of the two
    embedded propagators of the RKF7(8) method are created, with a time-step of 120 seconds.
    
    .. code-block:: python
    
      # Create 8th order RKF settings
      integrator_settings = integrator.runge_kutta_fixed_step(
          time_step = 120.0,
          coefficient_set = integrator.CoefficientSets.rkf_78,
          order_to_use = integrator.OrderToIntegrate.higher )"""

@typing.overload
def runge_kutta_fixed_step_size(initial_time: time_representation.Time, initial_time_step: time_representation.Time, coefficient_set: CoefficientSets, order_to_use: OrderToIntegrate=..., assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

@typing.overload
def runge_kutta_fixed_step_size(initial_time_step: time_representation.Time, coefficient_set: CoefficientSets, order_to_use: OrderToIntegrate=..., assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

def runge_kutta_variable_step(initial_time_step: time_representation.Time, coefficient_set: CoefficientSets, step_size_control_settings: IntegratorStepSizeControlSettings, step_size_validation_settings: IntegratorStepSizeValidationSettings, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    """Creates the settings for the Runge-Kutta variable step size integrator.
    
    Function to create settings for the Runge-Kutta variable step size integrator.
    Different coefficient sets (Butcher's tableau) can be used (see the :class:`CoefficientSets` enum).
    The step-size control algorithm is defined by a :class:`~tudatpy.dynamics.propagation_setup.integrator.IntegratorStepSizeControlSettings` and
    :class:`~tudatpy.dynamics.propagation_setup.integrator.IntegratorStepSizeValidationSettings` object, created using one of the functions
    listed above.
    
    
    Parameters
    ----------
    initial_time_step : float
        Initial time step to be used.
    coefficient_set : CoefficientSets
        Coefficient set (Butcher's tableau) to be used in the integration.
    step_size_control_settings : IntegratorStepSizeControlSettings
        Object used to control the step size, by computing a new step size :math:`\\Delta t_{rec.}`, from the embedded Runge-Kutta integrator pair,
        and recommending whether the steps is to be accepted, or recomputed with a different time step.
    
    step_size_validation_settings : IntegratorStepSizeValidationSettings
        Object used to validate whether the :math:`\\Delta t_{rec.}` provided by model defined by the ``step_size_control_settings`` meets with user-defined
        criteria (minimum, maximum values, etc.)
    
    assess_termination_on_minor_steps : bool, default=false
        Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
        integrator (true) or only at the end of each integration step (false).
    
    Returns
    -------
    IntegratorSettings
        Object containing settings for the integrator.
    
    
    
    
    
    Examples
    --------
    In this example, settings for the variable step RK4(5) integrator are created, with the same tolerances (:math:`10^{-10}`)
    applied element-wise on the propagated state. The minimum and maximum time steps are set to 0.001 and 1000 seconds,
    and the initial step is set to 30 seconds. All other inputs are left on their defaults
    
    .. code-block:: python
    
      # Create RK4(5) settings
      control_settings = integrator.step_size_control_elementwise_scalar_tolerance( 1.0E-10, 1.0E-10 )
      validation_settings = integrator.step_size_validation( 0.001, 1000.0 )
      integrator_settings = integrator.runge_kutta_variable_step(
          initial_time_step = 30.0,
          coefficient_set = integrator.CoefficientSets.rkf_45,
          step_size_control_settings = control_settings,
          step_size_validation_settings = validation_settings )
    
    In this example, the above is modified such that step-size control is applied on position and velocity
    element blocks.
    
    .. code-block:: python
    
      # Create RK4(5) settings
      control_settings = integrator.step_size_control_custom_blockwise_scalar_tolerance(
          integrator.standard_cartesian_state_element_blocks
          1.0E-10, 1.0E-10 )
      validation_settings = integrator.step_size_validation( 0.001, 1000.0 )
      integrator_settings = integrator.runge_kutta_variable_step(
          initial_time_step = 30.0,
          coefficient_set = integrator.CoefficientSets.rkf_45,
          step_size_control_settings = control_settings,
          step_size_validation_settings = validation_settings )"""

@typing.overload
def runge_kutta_variable_step_size(initial_time_step: time_representation.Time, coefficient_set: CoefficientSets, minimum_step_size: time_representation.Time, maximum_step_size: time_representation.Time, relative_error_tolerance: float, absolute_error_tolerance: float, assess_termination_on_minor_steps: bool=False, safety_factor: time_representation.Time=0.8, maximum_factor_increase: time_representation.Time=4.0, minimum_factor_increase: time_representation.Time=0.1, throw_exception_if_minimum_step_exceeded: bool=True) -> IntegratorSettings:
    """Creates the settings for the Runge-Kutta variable step size integrator with scalar tolerances.
    
    NOTE: THIS FUNCTION IS DEPRECATED, IT IS RECOMMENDED TO USE THE NEW :func:`~tudatpy.dynamics.propagation_setup.integrator.runge_kutta_variable_step` INTERFACE INSTEAD
    
    Function to create settings for the Runge-Kutta variable step size integrator with scalar tolerances.
    For this integrator, the step size is varied based on the tolerances and safety factor provided.
    The tolerance is composed of an absolute and a relative part.
    Different coefficient sets (Butcher's tableau) can be used (see the `CoefficientSets` enum).
    
    
    Parameters
    ----------
    initial_time_step : float
        Initial time step to be used.
    coefficient_set : CoefficientSets
        Coefficient set (Butcher's tableau) to be used in the integration.
    minimum_step_size : float
        Minimum time step to be used during the integration.
    maximum_step_size : float
        Maximum time step to be used during the integration.
    relative_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
        Relative vector tolerance to adjust the time step.
    absolute_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
        Absolute vector tolerance to adjust the time step.
    assess_termination_on_minor_steps : bool, default=false
        Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
        integrator (true) or only at the end of each integration step (false).
    
    safety_factor : float, default=0.8
        Safety factor used in the step size control.
    maximum_factor_increase : float, default=4.0
        Maximum increase between consecutive time steps, expressed as the factor between new and old step size.
    
    minimum_factor_increase : float, default=0.1
        Minimum increase between consecutive time steps, expressed as the factor between new and old step size.
    
    throw_exception_if_minimum_step_exceeded : bool, default=true
        If set to false, the variable step integrator will use the minimum step size specified when the algorithm
        computes the optimum one to be lower, instead of throwing an exception.
    
    Returns
    -------
    RungeKuttaVariableStepSettingsScalarTolerances
        RungeKuttaVariableStepSettingsScalarTolerances object."""

@typing.overload
def runge_kutta_variable_step_size(initial_time: time_representation.Time, initial_time_step: time_representation.Time, coefficient_set: CoefficientSets, minimum_step_size: time_representation.Time, maximum_step_size: time_representation.Time, relative_error_tolerance: float, absolute_error_tolerance: float, assess_termination_on_minor_steps: bool=False, safety_factor: time_representation.Time=0.8, maximum_factor_increase: time_representation.Time=4.0, minimum_factor_increase: time_representation.Time=0.1, throw_exception_if_minimum_step_exceeded: bool=True) -> IntegratorSettings:
    ...

@typing.overload
def runge_kutta_variable_step_size_vector_tolerances(initial_time_step: time_representation.Time, coefficient_set: CoefficientSets, minimum_step_size: time_representation.Time, maximum_step_size: time_representation.Time, relative_error_tolerance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], absolute_error_tolerance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], assess_termination_on_minor_steps: bool=False, safety_factor: time_representation.Time=0.8, maximum_factor_increase: time_representation.Time=4.0, minimum_factor_increase: time_representation.Time=0.1, throw_exception_if_minimum_step_exceeded: bool=True) -> IntegratorSettings:
    """Creates the settings for the Runge-Kutta variable step size integrator with vector tolerances.
    
    NOTE: THIS FUNCTION IS DEPRECATED, IT IS RECOMMENDED TO USE THE NEW :func:`~tudatpy.dynamics.propagation_setup.integrator.runge_kutta_variable_step` INTERFACE INSTEAD
    
    Function to create settings for the Runge-Kutta variable step size integrator with vector tolerances.
    For this integrator, the step size is varied based on the tolerances and safety factor provided.
    The tolerance is composed of an absolute and a relative part.
    Different coefficient sets (Butcher's tableau) can be used (see the `CoefficientSets` enum).
    
    
    Parameters
    ----------
    initial_time_step : float
        Initial time step to be used.
    coefficient_set : CoefficientSets
        Coefficient set (Butcher's tableau) to be used in the integration.
    minimum_step_size : float
        Minimum time step to be used during the integration.
    maximum_step_size : float
        Maximum time step to be used during the integration.
    relative_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
        Relative vector tolerance to adjust the time step.
    absolute_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
        Absolute vector tolerance to adjust the time step.
    assess_termination_on_minor_steps : bool, default=false
        Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
        integrator (true) or only at the end of each integration step (false).
    
    safety_factor : float, default=0.8
        Safety factor used in the step size control.
    maximum_factor_increase : float, default=4.0
        Maximum increase between consecutive time steps, expressed as the factor between new and old step size.
    
    minimum_factor_increase : float, default=0.1
        Minimum increase between consecutive time steps, expressed as the factor between new and old step size.
    
    throw_exception_if_minimum_step_exceeded : bool, default=true
        If set to false, the variable step integrator will use the minimum step size specified when the algorithm
        computes the optimum one to be lower, instead of throwing an exception.
    
    Returns
    -------
    RungeKuttaVariableStepSizeSettingsVectorTolerances
        RungeKuttaVariableStepSizeSettingsVectorTolerances object."""

@typing.overload
def runge_kutta_variable_step_size_vector_tolerances(initial_time: time_representation.Time, initial_time_step: time_representation.Time, coefficient_set: CoefficientSets, minimum_step_size: time_representation.Time, maximum_step_size: time_representation.Time, relative_error_tolerance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], absolute_error_tolerance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], assess_termination_on_minor_steps: bool=False, safety_factor: time_representation.Time=0.8, maximum_factor_increase: time_representation.Time=4.0, minimum_factor_increase: time_representation.Time=0.1, throw_exception_if_minimum_step_exceeded: bool=True) -> IntegratorSettings:
    ...

def standard_cartesian_state_element_blocks(number_of_rows: int, number_of_columns: int) -> list[tuple[int, int, int, int]]:
    """Function to generate step size control blocks on position and velocity elements for numerical integration
    
    Function to generate step size control blocks on position and velocity elements for numerical integration, typically provided
    to the :func:`~tudatpy.dynamics.propagation_setup.integrator.step_size_control_custom_blockwise_scalar_tolerance` or
    :func:`~tudatpy.dynamics.propagation_setup.integrator.step_size_control_custom_blockwise_matrix_tolerance` function.
    By providing this function to one of these step-size control functions, the final column of the state vector is taken (such that
    it works  both for state-only, and variational equations and state propagation) and combined into :math:`N` blocks of size 3.
    The step-size control is then done on each of these blocks, which will represent the position and velocity blocks.
    
    
    Parameters
    ----------
    number_of_rows : int
        Number of rows in state vector
    number_of_columns : int
        Number of columns in state vector
    Returns
    -------
    list[tuple[int,int,int,int]]
        List of matrix blocks over which the step size control is to be done (see ``block_indices_function`` input to :func:`~tudatpy.dynamics.propagation_setup.integrator.step_size_control_custom_blockwise_scalar_tolerance`)"""

def standard_rotational_state_element_blocks(number_of_rows: int, number_of_columns: int) -> list[tuple[int, int, int, int]]:
    ...

def step_size_control_blockwise_matrix_tolerance(block_indices: list[tuple[int, int, int, int]], relative_error_tolerance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], absolute_error_tolerance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], safety_factor: float=0.8, minimum_factor_increase: float=0.1, maximum_factor_increase: float=4.0) -> IntegratorStepSizeControlSettings:
    """Creates settings for integrator step-size control, using block-wise analysis for the propagated states.
    
    Function to create settings for integrator step-size control, using block-wise analysis for the propagated states. This function
    is similar to :func:`~tudatpy.dynamics.propagation_setup.integrator.step_size_control_blockwise_scalar_tolerance`,
    with the differences that the tolerances are provided as a list (which must be of equal size as the number of state blocks used), such that
    different tolerances can be provided for each state block.
    
    If the size of the tolerances used as input differ from one another, or differ from the number of blocks, an exception is thrown
    
    
    .. note::
    
       If you would like to create block indices that group the position and velocity elements, take a look at the :func:`~tudatpy.dynamics.propagation_setup.integrator.standard_cartesian_state_element_blocks` function.
       The function will create a list of blocks that can be used as input to the `block_indices` argument of this function.
    
    Parameters
    ----------
    block_indices : list[tuple[int,int,int,int]]
        List of matrix blocks over which the norms are to be taken (with entries of the tuple denoting :math:`i,j,k,l`, in order.
    relative_error_tolerance : numpy.ndarray[numpy.float64[m, 1]]
        Values of relative error tolerance :math:`\\boldsymbol{\\epsilon}_{r}`.
    absolute_error_tolerance : numpy.ndarray[numpy.float64[m, 1]]
        Values of absolute error tolerance :math:`\\boldsymbol{\\epsilon}_{a}`.
    safety_factor : float, default = 0.8
        Safety factor :math:`K` for step size control
    minimum_factor_increase : float, default = 0.1
        Minimum permissible value for :math:`\\Delta t_{rec.}/\\Delta t`
    maximum_factor_increase : float, default = 4.0
        Maximum permissible value for :math:`\\Delta t_{rec.}/\\Delta t`
    Returns
    -------
    IntegratorStepSizeControlSettings
        Object containing settings for per-element step-size control.
    
    Examples
    --------
    In this example, step size control settings are created for a Cartesian state vector, which group the position and velocity elements for the step size validation. Note, these block indices can also be conveniently created using the :func:`~tudatpy.dynamics.propagation_setup.integrator.standard_cartesian_state_element_blocks` function.
    Here we will create them manually for demonstration purposes.
    We would like to create integrator settings with the following settings:
    
    - Relative error tolerance of 1e-12 for position and velocity blocks
    - Absolute error tolerance of 1e-9 for position and 1e-12 for velocity blocks
    
    .. code-block:: python
    
        import numpy as np
        from tudatpy.dynamics import propagation_setup
        ...
        # Define integrator step settings
        initial_time_step = 10.0
        minimum_step_size = 1.0e-12
        maximum_step_size = 60.0
        relative_tolerance_pos = 1.0e-12
        relative_tolerance_vel = 1.0e-12
        absolute_tolerance_pos = 1.0e-9
        absolute_tolerance_vel = 1.0e-12
        \"\"\"
        Our Cartesian state y is a 2D vector, with dimensions [6, 1].
        # column-index
            0
        y = [[ x ], # row-index 0
            [ y ], # row-index 1
            [ z ], # row-index 2
            [ vx ], # row-index 3
            [ vy ], # row-index 4
            [ vz ]] # row-index 5
        The block indices are denoted as i, j, k, l, where:
        - i: start row index of the block
        - j: start column index of the block
        - k: number of rows in the block
        - l: number of columns in the block
        The corresponding block indices for the position block are therefore
        - i = 0 (start the block at row-index 0)
        - j = 0 (start the block at column-index 0)
        - k = 3 (the block has 3 rows: x, y, z)
        - l = 1 (the block has 1 column)
        which gives us the block indices (i=0, j=0, k=3, l=1).
        For the velocity block, the indices are:
        - i = 3 (start the block at row-index 3)
        - j = 0 (start the block at column-index 0)
        - k = 3 (the block has 3 rows: vx, vy, vz)
        - l = 1 (the block has 1 column)
        which gives us the block indices (i=3, j=0, k=3, l=1).
        \"\"\"
        # Manually define block indices for position and velocity,
        # which is equivalent to the standard block indices using:
        # block_indices = propagation_setup.integrator.standard_cartesian_state_element_blocks(6, 1)
        block_indices = [(0, 0, 3, 1), (3, 0, 3, 1)]
        # Different from the scalar tolerance, the matrix tolerance is defined as
        # the relative and absolute tolerances for each block.
        relative_tolerances = np.array([[relative_tolerance_pos], [relative_tolerance_vel]])
        absolute_tolerances = np.array([[absolute_tolerance_pos], [absolute_tolerance_vel]])
        step_size_control_settings = (
            propagation_setup.integrator.step_size_control_blockwise_matrix_tolerance(
                block_indices, relative_tolerances, absolute_tolerances
            )
        )
        step_size_validation_settings = propagation_setup.integrator.step_size_validation(
            minimum_step=minimum_step_size, maximum_step=maximum_step_size
        )
        # Retrieve coefficient set
        coefficient_set = propagation_setup.integrator.rkf_78
        variable_step_integrator_settings = (
            propagation_setup.integrator.runge_kutta_variable_step(
                initial_time_step,
                coefficient_set,
                step_size_control_settings,
                step_size_validation_settings,
            )
        )
        ..."""

def step_size_control_blockwise_scalar_tolerance(block_indices: list[tuple[int, int, int, int]], relative_error_tolerance: float, absolute_error_tolerance: float, safety_factor: float=0.8, minimum_factor_increase: float=0.1, maximum_factor_increase: float=4.0) -> IntegratorStepSizeControlSettings:
    """Creates settings for integrator step-size control, using block-wise analysis for the propagated states.
    
    Function to create settings for integrator step-size control, using block-wise analysis for the propagated states. This function
    is similar to :func:`~tudatpy.dynamics.propagation_setup.integrator.step_size_control_elementwise_scalar_tolerance`,
    with the difference that the error estimation :math:`\\boldsymbol{\\epsilon}` is not used on an element-by-element basis, but using the norms
    of user defined matrix blocks. This is for instance very useful when propagating Cartesian states, where the tolerances are then
    typically applied twice: once to the norm of the position error, and once to the norm of the velocity error.
    
    The algorithm is then run, using the modification
    that :math:`\\epsilon_{i}\\rightarrow||\\boldsymbol{\\epsilon_{[i,k],[j,l]}}||`. Where the indices on the right-hand side denote start row :math:`i`,
    start column :math:`j`, number of rows :math:`k` and number of columns :math:`l`. over
    which the state error norm is to be taken. For a single Cartesian state vector, the norm is taken on blocks :math:`[0,3],[0,1]` and :math:`[3,3],[0,1]`
    
    .. note::
    
        If you would like to create block indices that group the position and velocity elements, take a look at the :func:`~tudatpy.dynamics.propagation_setup.integrator.standard_cartesian_state_element_blocks` function.
        The function will create a list of blocks that can be used as input to the `block_indices` argument of this function.
    
    Parameters
    ----------
    block_indices : list[tuple[int,int,int,int]]
        List of matrix blocks over which the norms are to be taken (with entries of the tuple denoting :math:`i,j,k,l`, in order.
    relative_error_tolerance : float
        Value of relative error tolerance :math:`\\epsilon_{r}`.
    absolute_error_tolerance : float
        Value of absolute error tolerance :math:`\\epsilon_{a}`.
    safety_factor : float, default = 0.8
        Safety factor :math:`K` for step size control
    minimum_factor_increase : float, default = 0.1
        Minimum permissible value for :math:`\\Delta t_{rec.}/\\Delta t`
    maximum_factor_increase : float, default = 4.0
        Maximum permissible value for :math:`\\Delta t_{rec.}/\\Delta t`
    Returns
    -------
    IntegratorStepSizeControlSettings
        Object containing settings for per-element step-size control.
    
    
    Examples
    --------
    In this example, step size control settings are created for a Cartesian state vector, which group the position and velocity elements for the step size validation. Note, these block indices can also be conveniently created using the :func:`~tudatpy.dynamics.propagation_setup.integrator.standard_cartesian_state_element_blocks` function.
    Here we will create them manually for demonstration purposes.
    We would like to create integrator settings with the following settings:
    
    - Relative error tolerance of 1e-12
    - Absolute error tolerance of 1e-12
    - Relative and absolute error tolerances are applied to the position and velocity blocks of a single Cartesian state vector
    
    .. code-block:: python
    
        from tudatpy.dynamics import propagation_setup
        ...
        # Define integrator step settings
        initial_time_step = 10.0
        minimum_step_size = 1.0e-12
        maximum_step_size = 60.0
        relative_tolerance = 1.0e-12
        absolute_tolerance = 1.0e-12
        \"\"\"
        Our Cartesian state y is a 2D vector, with dimensions [6, 1].
        # column-index
               0
        y = [[ x ], # row-index 0
             [ y ], # row-index 1
             [ z ], # row-index 2
             [ vx ], # row-index 3
             [ vy ], # row-index 4
             [ vz ]] # row-index 5
        The block indices are denoted as i, j, k, l, where:
        - i: start row index of the block
        - j: start column index of the block
        - k: number of rows in the block
        - l: number of columns in the block
        The corresponding block indices for the position block are therefore
        - i = 0 (start the block at row-index 0)
        - j = 0 (start the block at column-index 0)
        - k = 3 (the block has 3 rows: x, y, z)
        - l = 1 (the block has 1 column)
        which gives us the block indices (i=0, j=0, k=3, l=1).
    
        For the velocity block, the indices are:
        - i = 3 (start the block at row-index 3)
        - j = 0 (start the block at column-index 0)
        - k = 3 (the block has 3 rows: vx, vy, vz)
        - l = 1 (the block has 1 column)
        which gives us the block indices (i=3, j=0, k=3, l=1).
        \"\"\"
        # Manually define block indices for position and velocity,
        # which is equivalent to the standard block indices using:
        # block_indices = propagation_setup.integrator.standard_cartesian_state_element_blocks(6, 1)
        block_indices = [(0, 0, 3, 1), (3, 0, 3, 1)]
        step_size_control_settings = (
            propagation_setup.integrator.step_size_control_blockwise_scalar_tolerance(
                block_indices, relative_tolerance, absolute_tolerance
            )
        )
        step_size_validation_settings = propagation_setup.integrator.step_size_validation(
            minimum_step=minimum_step_size, maximum_step=maximum_step_size
        )
        # Retrieve coefficient set
        coefficient_set = propagation_setup.integrator.rkf_78
        variable_step_integrator_settings = (
            propagation_setup.integrator.runge_kutta_variable_step(
                initial_time_step,
                coefficient_set,
                step_size_control_settings,
                step_size_validation_settings,
            )
        )
        ..."""

def step_size_control_custom_blockwise_matrix_tolerance(block_indices_function: typing.Callable[[int, int], list[tuple[int, int, int, int]]], relative_error_tolerance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], absolute_error_tolerance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], safety_factor: float=0.8, minimum_factor_increase: float=0.1, maximum_factor_increase: float=4.0) -> IntegratorStepSizeControlSettings:
    """Creates settings for integrator step-size control, using block-wise analysis for the propagated states.
    
    Function to create settings for integrator step-size control, using block-wise analysis for the propagated states. This function
    is similar to :func:`~tudatpy.dynamics.propagation_setup.integrator.step_size_control_custom_blockwise_scalar_tolerance`,
    but uses blockwise tolerances (as in :func:`~tudatpy.dynamics.propagation_setup.integrator.step_size_control_blockwise_matrix_tolerance`)
    
    
    
    Parameters
    ----------
    block_indices_function : Callable[[int,int],list[tuple[int,int,int,int]]]
        Function returning list of matrix blocks over which the norms are to be taken (with entries of the tuple denoting :math:`i,j,k,l`, in order, with number of rows and columns of propagated state as input.
    relative_error_tolerance : numpy.ndarray[numpy.float64[m, 1]]
        Values of relative error tolerance :math:`\\boldsymbol{\\epsilon}_{r}`.
    absolute_error_tolerance : numpy.ndarray[numpy.float64[m, 1]]
        Values of absolute error tolerance :math:`\\boldsymbol{\\epsilon}_{a}`.
    safety_factor : float, default = 0.8
        Safety factor :math:`K` for step size control
    minimum_factor_increase : float, default = 0.1
        Minimum permissible value for :math:`\\Delta t_{rec.}/\\Delta t`
    maximum_factor_increase : float, default = 4.0
        Maximum permissible value for :math:`\\Delta t_{rec.}/\\Delta t`
    Returns
    -------
    IntegratorStepSizeControlSettings
        Object containing settings for per-element step-size control."""

def step_size_control_custom_blockwise_scalar_tolerance(block_indices_function: typing.Callable[[int, int], list[tuple[int, int, int, int]]], relative_error_tolerance: float, absolute_error_tolerance: float, safety_factor: float=0.8, minimum_factor_increase: float=0.1, maximum_factor_increase: float=4.0) -> IntegratorStepSizeControlSettings:
    """Creates settings for integrator step-size control, using block-wise analysis for the propagated states.
    
    Function to create settings for integrator step-size control, using block-wise analysis for the propagated states. This function
    is similar to :func:`~tudatpy.dynamics.propagation_setup.integrator.step_size_control_blockwise_scalar_tolerance`,
    but rather than providing the ``block_indices`` directly, a function to determine the block indices, based on the size of the
    propagated state, is provided. For instance, the :func:`~tudatpy.dynamics.propagation_setup.integrator.standard_cartesian_state_element_blocks`
    can be provided to this function (as ``block_indices_function``), which will adapt the block indices depending on the size of the propagated state
    (e.g. regardless of how many bodies are propagated, step size control will always be done on position and velocity element blocks)
    
    
    
    
    Parameters
    ----------
    block_indices_function : Callable[[int,int],list[tuple[int,int,int,int]]]
        Function returning list of matrix blocks over which the norms are to be taken (with entries of the tuple denoting :math:`i,j,k,l`, in order, with number of rows and columns of propagated state as input.
    relative_error_tolerance : float
        Value of relative error tolerance :math:`\\epsilon_{r}`.
    absolute_error_tolerance : float
        Value of absolute error tolerance :math:`\\epsilon_{a}`.
    safety_factor : float, default = 0.8
        Safety factor :math:`K` for step size control
    minimum_factor_increase : float, default = 0.1
        Minimum permissible value for :math:`\\Delta t_{rec.}/\\Delta t`
    maximum_factor_increase : float, default = 4.0
        Maximum permissible value for :math:`\\Delta t_{rec.}/\\Delta t`
    Returns
    -------
    IntegratorStepSizeControlSettings
        Object containing settings for per-element step-size control."""

def step_size_control_elementwise_matrix_tolerance(relative_error_tolerance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], absolute_error_tolerance: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], safety_factor: float=0.8, minimum_factor_increase: float=0.1, maximum_factor_increase: float=4.0) -> IntegratorStepSizeControlSettings:
    """Creates settings for integrator step-size control, using element-wise analysis for the propagated states.
    
    Function to create settings for integrator step-size control, using element-wise analysis for the propagated states. This function
    is similar to :func:`~tudatpy.dynamics.propagation_setup.integrator.step_size_control_elementwise_scalar_tolerance`,
    with the differences that the tolerances are provided as a vector/matrix (which must be of equal size as the propagates state), such that
    different tolerances can be provided for each state element. The behaviour of the algorithm is then such that
    :math:`\\epsilon_{r}\\rightarrow\\epsilon_{r,i}` and :math:`\\epsilon_{a}\\rightarrow\\epsilon_{a,i}`.
    
    If the size of the tolerances used as input differ from one another, or differ from the size of the state vector, an exception is thrown
    
    
    Parameters
    ----------
    relative_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
        Values of relative error tolerance :math:`\\boldsymbol{\\epsilon}_{r}`.
    absolute_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
        Values of absolute error tolerance :math:`\\boldsymbol{\\epsilon}_{a}`.
    safety_factor : float, default = 0.8
        Safety factor :math:`K` for step size control
    minimum_factor_increase : float, default = 0.1
        Minimum permissible value for :math:`\\Delta t_{rec.}/\\Delta t`
    maximum_factor_increase : float, default = 4.0
        Maximum permissible value for :math:`\\Delta t_{rec.}/\\Delta t`
    Returns
    -------
    IntegratorStepSizeControlSettings
        Object containing settings for per-element step-size control."""

def step_size_control_elementwise_scalar_tolerance(relative_error_tolerance: float, absolute_error_tolerance: float, safety_factor: float=0.8, minimum_factor_increase: float=0.1, maximum_factor_increase: float=4.0) -> IntegratorStepSizeControlSettings:
    """Creates settings for integrator step-size control, using element-wise analysis for the propagated states.
    
    Function to create settings for integrator step-size control, using element-wise analysis for the propagated states. For a propagated
    state :math:`\\mathbf{x}` with entries :math:`x_{i}`, and an estimate :math:`\\boldsymbol{\\epsilon}` for the current local error, the following
    algorithm is performed per element :math:`i` to calculate the required error :math:`\\epsilon_{i,req}` on this element:
    
    .. math::
    
       \\epsilon_{i,req}=\\epsilon_{r}x_{i}+\\epsilon_{a}
    
    A proposed modification to the step size is then computed, using the most constraining of all state elements
    
    .. math::
    
       \\bar{\\Delta t_{rec.}}&=\\Delta t\\left(\\min_{i}\\left(\\frac{\\epsilon_{i,req}}{\\epsilon_{i}}\\right)\\right)^{p}\\\\
       \\Delta t_{rec.}&=K\\bar{\\Delta t_{rec.}}
    
    with :math:`p` the order of the local truncation error of the method for which step-size control is being applied,
    :math:`\\Delta t_{rec.}` the new, recommended step size, and :math:`\\Delta t` the current step size. The factor :math:`K` is a safety factor
    used make the time step slightly smaller than strictly required.
    
    A minimum and maximum change in time step may be provided by the user, such that if :math:`\\Delta t_{rec.}/\\Delta t` is too large or too small,
    the proposed increase/decrease to the step size is constrained to this limit value. That is, if :math:`\\Delta t_{rec.}/\\Delta t` proposed by the algorithm is
    1000, and the ``maximum_factor_increase`` input is equal to 20, the algorithm will use :math:`\\Delta t_{rec.}/\\Delta t=20` in what follows.
    
    For cases where :math:`\\bar{\\Delta t_{rec.}}/\\Delta t < 1`, the step is recommended to be recomputed with the new proposed step size (e.g. the current step
    is not accepted, and will be re-attempted with a smaller step size). For cases where :math:`\\bar{\\Delta t_{rec.}}/\\Delta t > 1`, the step is accepted, and
    the next step will be performed with the new, higher, step size.
    
    
    Parameters
    ----------
    relative_error_tolerance : float
        Value of relative error tolerance :math:`\\epsilon_{r}`.
    absolute_error_tolerance : float
        Value of absolute error tolerance :math:`\\epsilon_{a}`.
    safety_factor : float, default = 0.8
        Safety factor :math:`K` for step size control
    minimum_factor_increase : float, default = 0.1
        Minimum permissible value for :math:`\\Delta t_{rec.}/\\Delta t`
    maximum_factor_increase : float, default = 4.0
        Maximum permissible value for :math:`\\Delta t_{rec.}/\\Delta t`
    Returns
    -------
    IntegratorStepSizeControlSettings
        Object containing settings for per-element step-size control."""

def step_size_validation(minimum_step: float, maximum_step: float, minimum_step_size_handling: MinimumIntegrationTimeStepHandling=..., accept_infinity_step: bool=False, accept_nan_step: bool=False) -> IntegratorStepSizeValidationSettings:
    """Creates settings step size validation in a variable step-size integrator.
    
    Function to create settings step size validation in a variable step-size integrator. The validation
    model takes the proposed new step size  :math:`\\Delta t_{rec}` as input, and checks if it meets predefined conditions, specifically
    whether the proposed time step falls in a given predefined range :math:`[\\Delta t_{\\min}, \\Delta t_{\\max}]`.
    This function also provides the option of handling recommended step sizes below :math:`\\Delta t_{\\min}` in various ways,
    and control on how to deal with recommend Inf/NaN step sizes.
    
    
    Parameters
    ----------
    minimum_step : float
        Value of minimum permitted time step :math:`\\Delta t_{\\min}`.
    maximum_step : float
        Value of maximum permitted time step :math:`\\Delta t_{\\max}`.
    minimum_step_size_handling : MinimumIntegrationTimeStepHandling, default = throw_exception_below_minimum
        Entry defining the behaviour when :math:`\\Delta t_{rec}<\\Delta t_{\\min}`.
    accept_infinity_step : bool, default = False
        Entry defining whether to accept a step size of infinity (if False, exception is throw in such cases)
    accept_nan_step : bool, default = False
        Entry defining whether to accept a step size of NaN (if False, exception is throw in such cases)
    Returns
    -------
    IntegratorStepSizeValidationSettings
        Object containing settings for step-size validation."""
SSPRK3: CoefficientSets
adams_bashforth_moulton_type: AvailableIntegrators
bulirsch_stoer_sequence: ExtrapolationMethodStepSequences
bulirsch_stoer_type: AvailableIntegrators
deufelhard_sequence: ExtrapolationMethodStepSequences
euler_forward: CoefficientSets
explicit_mid_point: CoefficientSets
explicit_trapezoid_rule: CoefficientSets
heun_euler: CoefficientSets
higher: OrderToIntegrate
lower: OrderToIntegrate
ralston: CoefficientSets
ralston_3: CoefficientSets
ralston_4: CoefficientSets
rk_3: CoefficientSets
rk_4: CoefficientSets
rkdp_87: CoefficientSets
rkf_108: CoefficientSets
rkf_12: CoefficientSets
rkf_1210: CoefficientSets
rkf_1412: CoefficientSets
rkf_45: CoefficientSets
rkf_56: CoefficientSets
rkf_78: CoefficientSets
rkf_89: CoefficientSets
rkv_89: CoefficientSets
runge_kutta_fixed_step_size_type: AvailableIntegrators
runge_kutta_variable_step_size_type: AvailableIntegrators
set_to_minimum_step_every_time_warning: MinimumIntegrationTimeStepHandling
set_to_minimum_step_silently: MinimumIntegrationTimeStepHandling
set_to_minimum_step_single_warning: MinimumIntegrationTimeStepHandling
three_eight_rule_rk_4: CoefficientSets
throw_exception_below_minimum: MinimumIntegrationTimeStepHandling