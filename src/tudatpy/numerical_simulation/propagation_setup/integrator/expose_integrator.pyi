import numpy
import typing
__all__ = ['AdamsBashforthMoultonSettings', 'AvailableIntegrators', 'BulirschStoerIntegratorSettings', 'CoefficientSets', 'ExtrapolationMethodStepSequences', 'IntegratorSettings', 'IntegratorStepSizeControlSettings', 'IntegratorStepSizeValidationSettings', 'MinimumIntegrationTimeStepHandling', 'OrderToIntegrate', 'RungeKuttaFixedStepSizeSettings', 'RungeKuttaVariableStepSizeBaseSettings', 'RungeKuttaVariableStepSizeSettingsScalarTolerances', 'RungeKuttaVariableStepSizeSettingsVectorTolerances', 'SSPRK3', 'adams_bashforth_moulton', 'adams_bashforth_moulton_fixed_order', 'adams_bashforth_moulton_fixed_step', 'adams_bashforth_moulton_fixed_step_fixed_order', 'adams_bashforth_moulton_type', 'bulirsch_stoer', 'bulirsch_stoer_fixed_step', 'bulirsch_stoer_sequence', 'bulirsch_stoer_type', 'bulirsch_stoer_variable_step', 'deufelhard_sequence', 'euler', 'euler_forward', 'explicit_mid_point', 'explicit_trapezoid_rule', 'heun_euler', 'higher', 'lower', 'print_butcher_tableau', 'ralston', 'ralston_3', 'ralston_4', 'rk_3', 'rk_4', 'rkdp_87', 'rkf_108', 'rkf_12', 'rkf_1210', 'rkf_1412', 'rkf_45', 'rkf_56', 'rkf_78', 'rkf_89', 'rkv_89', 'runge_kutta_4', 'runge_kutta_fixed_step', 'runge_kutta_fixed_step_size', 'runge_kutta_fixed_step_size_type', 'runge_kutta_variable_step', 'runge_kutta_variable_step_size', 'runge_kutta_variable_step_size_type', 'runge_kutta_variable_step_size_vector_tolerances', 'set_to_minimum_step_every_time_warning', 'set_to_minimum_step_silently', 'set_to_minimum_step_single_warning', 'standard_cartesian_state_element_blocks', 'standard_rotational_state_element_blocks', 'step_size_control_blockwise_matrix_tolerance', 'step_size_control_blockwise_scalar_tolerance', 'step_size_control_custom_blockwise_matrix_tolerance', 'step_size_control_custom_blockwise_scalar_tolerance', 'step_size_control_elementwise_matrix_tolerance', 'step_size_control_elementwise_scalar_tolerance', 'step_size_validation', 'three_eight_rule_rk_4', 'throw_exception_below_minimum']

class AdamsBashforthMoultonSettings(IntegratorSettings):
    """`IntegratorSettings`-derived class to define settings for Adams-Bashforth-Moulton integrator settings.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class AvailableIntegrators:
    """Enumeration of integrators available with tudat.
	
	
	:member runge_kutta_fixed_step_size_type:
	:member runge_kutta_variable_step_size_type:
	:member bulirsch_stoer_type:
	:member adams_bashforth_moulton_type:
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, AvailableIntegrators]]
    adams_bashforth_moulton_type: typing.ClassVar[AvailableIntegrators]
    bulirsch_stoer_type: typing.ClassVar[AvailableIntegrators]
    runge_kutta_fixed_step_size_type: typing.ClassVar[AvailableIntegrators]
    runge_kutta_variable_step_size_type: typing.ClassVar[AvailableIntegrators]

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

class BulirschStoerIntegratorSettings(IntegratorSettings):
    """`IntegratorSettings`-derived class to define settings for Bulirsch-Stoer integrator settings.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class CoefficientSets:
    """Coefficient sets for Runge-Kutta-type integrators.
	
	Coefficient sets for Runge-Kutta-type integrators. The coefficients are defined
	in a Butcher Tableau, with an coefficient set yielding an x(y) method yielding an integrator
	with global truncation error of :math:`O(\\Delta t^{x})`. Some of these coefficients also contain an embedded integrator of :math:`O(\\Delta t^{y})`
	for step size control.
	
	
	:member euler_forward:
	:member rk_4:
	:member explicit_mid_point:
	:member explicit_trapezoid_rule:
	:member ralston:
	:member rk_3:
	:member ralston_3:
	:member SSPRK3:
	:member ralston_4:
	:member three_eight_rule_rk_4:
	:member rkf_12:
	:member heun_euler:
	:member rkf_45:
	:member rkf_56:
	:member rkf_78:
	:member rkdp_87:
	:member rkf_89:
	:member rkv_89:
	:member rkf_108:
	:member rkf_1210:
	:member rkf_1412:
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
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

class ExtrapolationMethodStepSequences:
    """Enumeration of available extrapolation method step sequences.
	
	
	:member bulirsch_stoer_sequence:
	:member deufelhard_sequence:
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, ExtrapolationMethodStepSequences]]
    bulirsch_stoer_sequence: typing.ClassVar[ExtrapolationMethodStepSequences]
    deufelhard_sequence: typing.ClassVar[ExtrapolationMethodStepSequences]

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

class IntegratorSettings:
    """Functional base class to define settings for integrators.
	
	Class to define settings for numerical integrators, for instance for use in numerical integration of equations of motion/
	variational equations. This class can be used for simple integrators such as fixed step RK and Euler. Integrators that
	require more settings to define have their own derived class.
	"""
    initial_time: float

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class IntegratorStepSizeControlSettings:
    """
		"""
    maximum_step_decrease: float
    minimum_step_decrease: float
    safety_factor: float

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class IntegratorStepSizeValidationSettings:
    """
		"""
    maximum_step: float
    minimum_step: float
    minimum_step_handling: MinimumIntegrationTimeStepHandling

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class MinimumIntegrationTimeStepHandling:
    """Members:
	
	throw_exception_below_minimum :
	
	set_to_minimum_step_silently :
	
	set_to_minimum_step_single_warning :
	
	set_to_minimum_step_every_time_warning :
	"""
    __members__: typing.ClassVar[dict[str, MinimumIntegrationTimeStepHandling]]
    set_to_minimum_step_every_time_warning: typing.ClassVar[MinimumIntegrationTimeStepHandling]
    set_to_minimum_step_silently: typing.ClassVar[MinimumIntegrationTimeStepHandling]
    set_to_minimum_step_single_warning: typing.ClassVar[MinimumIntegrationTimeStepHandling]
    throw_exception_below_minimum: typing.ClassVar[MinimumIntegrationTimeStepHandling]

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

class OrderToIntegrate:
    """Which integrator order needs to be integrated, only for coefficient sets with an embedded order.
	
	
	:member lower:
	:member higher:
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, OrderToIntegrate]]
    higher: typing.ClassVar[OrderToIntegrate]
    lower: typing.ClassVar[OrderToIntegrate]

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

class RungeKuttaFixedStepSizeSettings(IntegratorSettings):
    """`IntegratorSettings`-derived class to define settings for Runge Kutta integrators with a fixed step size
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class RungeKuttaVariableStepSizeBaseSettings(IntegratorSettings):
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class RungeKuttaVariableStepSizeSettingsScalarTolerances(RungeKuttaVariableStepSizeBaseSettings):
    """`IntegratorSettings`-derived class to define settings for Runge Kutta integrators with scalar tolerances.
	
	This
	class is actually derived from an intermediate class (`RungeKuttaVariableStepSizeBaseSettings`, not documented),
	which is derived directly from `IntegratorSettings`.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class RungeKuttaVariableStepSizeSettingsVectorTolerances(RungeKuttaVariableStepSizeBaseSettings):
    """`IntegratorSettings`-derived class to define settings for Runge Kutta integrators with vector tolerances.
	
	This class is actually derived from an intermediate class (`RungeKuttaVariableStepSizeBaseSettings`, not documented),
	which is derived directly from `IntegratorSettings`.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

@typing.overload
def adams_bashforth_moulton(initial_time_step: float, minimum_step_size: float, maximum_step_size: float, relative_error_tolerance: float=1e-12, absolute_error_tolerance: float=1e-12, minimum_order: int=6, maximum_order: int=11, assess_termination_on_minor_steps: bool=False, bandwidth: float=200.0) -> IntegratorSettings:
    """Creates the settings for the Adams-Bashforth-Moulton integrator.
	
	Factory function to create settings for the Adams-Bashforth-Moulton integrator.
	For this integrator, the step size is varied based on the tolerances and safety factor provided.
	The tolerance is composed of an absolute and a relative part.
	
	
	:param initial_time_step:
			Initial time step to be used.
	:param minimum_step_size:
			Minimum time step to be used during the integration.
	:param maximum_step_size:
			Maximum time step to be used during the integration.
	:param relative_error_tolerance:
			Relative tolerance to adjust the time step.
	:param absolute_error_tolerance:
			Relative tolerance to adjust the time step.
	:param minimum_order:
			Minimum order of the integrator.
	:param maximum_order:
			Maximum order of the integrator.
	:param save_frequency:
			Frequency at which to save the numerical integrated states (expressed per unit integration time step, with n = saveFrequency, so n = 2 means that the state is saved every two integration steps).
	:param assess_termination_on_minor_steps:
			Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
	:param bandwidth:
			Maximum error factor for doubling the step size.
	:return:
			AdamsBashforthMoultonSettings object.
	"""

@typing.overload
def adams_bashforth_moulton(initial_time: float, initial_time_step: float, minimum_step_size: float, maximum_step_size: float, relative_error_tolerance: float=1e-12, absolute_error_tolerance: float=1e-12, minimum_order: int=6, maximum_order: int=11, assess_termination_on_minor_steps: bool=False, bandwidth: float=200.0) -> IntegratorSettings:
    ...

def adams_bashforth_moulton_fixed_order(initial_time_step: float, minimum_step_size: float, maximum_step_size: float, relative_error_tolerance: float=1e-12, absolute_error_tolerance: float=1e-12, order: int=6, assess_termination_on_minor_steps: bool=False, bandwidth: float=200.0) -> IntegratorSettings:
    ...

def adams_bashforth_moulton_fixed_step(time_step: float, relative_error_tolerance: float=1e-12, absolute_error_tolerance: float=1e-12, minimum_order: int=6, maximum_order: int=11, assess_termination_on_minor_steps: bool=False, bandwidth: float=200.0) -> IntegratorSettings:
    ...

def adams_bashforth_moulton_fixed_step_fixed_order(time_step: float, order: int=6, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

@typing.overload
def bulirsch_stoer(initial_time: float, initial_time_step: float, extrapolation_sequence: ExtrapolationMethodStepSequences, maximum_number_of_steps: int, minimum_step_size: float, maximum_step_size: float, relative_error_tolerance: float=1e-12, absolute_error_tolerance: float=1e-12, assess_termination_on_minor_steps: bool=False, safety_factor: float=0.7, maximum_factor_increase: float=10.0, minimum_factor_increase: float=0.1) -> IntegratorSettings:
    ...

@typing.overload
def bulirsch_stoer(initial_time_step: float, extrapolation_sequence: ExtrapolationMethodStepSequences, maximum_number_of_steps: int, minimum_step_size: float, maximum_step_size: float, relative_error_tolerance: float=1e-12, absolute_error_tolerance: float=1e-12, assess_termination_on_minor_steps: bool=False, safety_factor: float=0.7, maximum_factor_increase: float=10.0, minimum_factor_increase: float=0.1) -> IntegratorSettings:
    """Creates the settings for the Bulirsch-Stoer integrator.
	
	Factory function to create settings for the Bulirsch-Stoer integrator.
	For this integrator, the step size is varied based on the tolerances and safety factor provided.
	The tolerance is composed of an absolute and a relative part.
	Different extrapolation sequences can be used (see the `ExtrapolationMethodStepSequences` enum).
	
	
	:param initial_time_step:
			Initial time step to be used.
	:param extrapolation_sequence:
			Extrapolation sequence to be used in the integration.
	:param maximum_number_of_steps:
			Number of entries in the sequence (e.g., number of integrations used for a single extrapolation).
	:param minimum_step_size:
			Minimum time step to be used during the integration.
	:param maximum_step_size:
			Maximum time step to be used during the integration.
	:param relative_error_tolerance:
			Relative tolerance to adjust the time step.
	:param absolute_error_tolerance:
			Relative tolerance to adjust the time step.
	:param save_frequency:
			Frequency at which to save the numerical integrated states (expressed per unit integration time step, with n = saveFrequency, so n = 2 means that the state is saved every two integration steps).
	:param assess_termination_on_minor_steps:
			Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
	:param safety_factor:
			Safety factor used in the step size control.
	:param maximum_factor_increase:
			Maximum increase between consecutive time steps, expressed as the factor between new and old step size.
	:param minimum_factor_increase:
			Minimum increase between consecutive time steps, expressed as the factor between new and old step size.
	:return:
			BulirschStoerIntegratorSettings object.
	"""

def bulirsch_stoer_fixed_step(time_step: float, extrapolation_sequence: ExtrapolationMethodStepSequences, maximum_number_of_steps: int, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

def bulirsch_stoer_variable_step(initial_time_step: float, extrapolation_sequence: ExtrapolationMethodStepSequences, maximum_number_of_steps: int, step_size_control_settings: IntegratorStepSizeControlSettings, step_size_validation_settings: IntegratorStepSizeValidationSettings, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

@typing.overload
def euler(initial_time: float, initial_time_step: float, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

@typing.overload
def euler(initial_time_step: float, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

def print_butcher_tableau(arg0: CoefficientSets) -> None:
    """Print the Butcher tableau of a given coefficient set.
	
	:param coefficient_set:
			Coefficient set of which the Butcher tableau will be printed.
	"""

@typing.overload
def runge_kutta_4(initial_time: float, initial_time_step: float, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

@typing.overload
def runge_kutta_4(initial_time_step: float, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

def runge_kutta_fixed_step(time_step: float, coefficient_set: CoefficientSets, order_to_use: OrderToIntegrate=..., assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

@typing.overload
def runge_kutta_fixed_step_size(initial_time: float, initial_time_step: float, coefficient_set: CoefficientSets, order_to_use: OrderToIntegrate=..., assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

@typing.overload
def runge_kutta_fixed_step_size(initial_time_step: float, coefficient_set: CoefficientSets, order_to_use: OrderToIntegrate=..., assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

def runge_kutta_variable_step(initial_time_step: float, coefficient_set: CoefficientSets, step_size_control_settings: IntegratorStepSizeControlSettings, step_size_validation_settings: IntegratorStepSizeValidationSettings, assess_termination_on_minor_steps: bool=False) -> IntegratorSettings:
    ...

@typing.overload
def runge_kutta_variable_step_size(initial_time_step: float, coefficient_set: CoefficientSets, minimum_step_size: float, maximum_step_size: float, relative_error_tolerance: float, absolute_error_tolerance: float, assess_termination_on_minor_steps: bool=False, safety_factor: float=0.8, maximum_factor_increase: float=4.0, minimum_factor_increase: float=0.1, throw_exception_if_minimum_step_exceeded: bool=True) -> IntegratorSettings:
    """Creates the settings for the Runge-Kutta variable step size integrator with scalar tolerances.
	
	Factory function to create settings for the Runge-Kutta variable step size integrator with scalar tolerances.
	For this integrator, the step size is varied based on the tolerances and safety factor provided.
	The tolerance is composed of an absolute and a relative part.
	Different coefficient sets (Butcher's tableau) can be used (see the `CoefficientSets` enum).
	
	
	:param initial_time_step:
			Initial time step to be used.
	:param coefficient_set:
			Coefficient set (Butcher's tableau) to be used in the integration.
	:param minimum_step_size:
			Minimum time step to be used during the integration.
	:param maximum_step_size:
			Maximum time step to be used during the integration.
	:param relative_error_tolerance:
			Relative vector tolerance to adjust the time step.
	:param absolute_error_tolerance:
			Absolute vector tolerance to adjust the time step.
	:param save_frequency:
			Frequency at which to save the numerical integrated states (expressed per unit integration time step,
			with n = saveFrequency, so n = 2 means that the state is saved every two integration steps).
	
	:param assess_termination_on_minor_steps:
			Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
			integrator (true) or only at the end of each integration step (false).
	
	:param safety_factor:
			Safety factor used in the step size control.
	:param maximum_factor_increase:
			Maximum increase between consecutive time steps, expressed as the factor between new and old step size.
	
	:param minimum_factor_increase:
			Minimum increase between consecutive time steps, expressed as the factor between new and old step size.
	
	:param throw_exception_if_minimum_step_exceeded:
			If set to false, the variable step integrator will use the minimum step size specified when the algorithm
			computes the optimum one to be lower, instead of throwing an exception.
	
	:return:
			RungeKuttaVariableStepSettingsScalarTolerances object.
	"""

@typing.overload
def runge_kutta_variable_step_size(initial_time: float, initial_time_step: float, coefficient_set: CoefficientSets, minimum_step_size: float, maximum_step_size: float, relative_error_tolerance: float, absolute_error_tolerance: float, assess_termination_on_minor_steps: bool=False, safety_factor: float=0.8, maximum_factor_increase: float=4.0, minimum_factor_increase: float=0.1, throw_exception_if_minimum_step_exceeded: bool=True) -> IntegratorSettings:
    ...

@typing.overload
def runge_kutta_variable_step_size_vector_tolerances(initial_time_step: float, coefficient_set: CoefficientSets, minimum_step_size: float, maximum_step_size: float, relative_error_tolerance: numpy.ndarray, absolute_error_tolerance: numpy.ndarray, assess_termination_on_minor_steps: bool=False, safety_factor: float=0.8, maximum_factor_increase: float=4.0, minimum_factor_increase: float=0.1, throw_exception_if_minimum_step_exceeded: bool=True) -> IntegratorSettings:
    """Creates the settings for the Runge-Kutta variable step size integrator with vector tolerances.
	
	Factory function to create settings for the Runge-Kutta variable step size integrator with vector tolerances.
	For this integrator, the step size is varied based on the tolerances and safety factor provided.
	The tolerance is composed of an absolute and a relative part.
	Different coefficient sets (Butcher's tableau) can be used (see the `CoefficientSets` enum).
	
	
	:param initial_time_step:
			Initial time step to be used.
	:param coefficient_set:
			Coefficient set (Butcher's tableau) to be used in the integration.
	:param minimum_step_size:
			Minimum time step to be used during the integration.
	:param maximum_step_size:
			Maximum time step to be used during the integration.
	:param relative_error_tolerance:
			Relative vector tolerance to adjust the time step.
	:param absolute_error_tolerance:
			Absolute vector tolerance to adjust the time step.
	:param save_frequency:
			Frequency at which to save the numerical integrated states (expressed per unit integration time step,
			with n = saveFrequency, so n = 2 means that the state is saved every two integration steps).
	
	:param assess_termination_on_minor_steps:
			Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
			integrator (true) or only at the end of each integration step (false).
	
	:param safety_factor:
			Safety factor used in the step size control.
	:param maximum_factor_increase:
			Maximum increase between consecutive time steps, expressed as the factor between new and old step size.
	
	:param minimum_factor_increase:
			Minimum increase between consecutive time steps, expressed as the factor between new and old step size.
	
	:param throw_exception_if_minimum_step_exceeded:
			If set to false, the variable step integrator will use the minimum step size specified when the algorithm
			computes the optimum one to be lower, instead of throwing an exception.
	
	:return:
			RungeKuttaVariableStepSizeSettingsVectorTolerances object.
	"""

@typing.overload
def runge_kutta_variable_step_size_vector_tolerances(initial_time: float, initial_time_step: float, coefficient_set: CoefficientSets, minimum_step_size: float, maximum_step_size: float, relative_error_tolerance: numpy.ndarray, absolute_error_tolerance: numpy.ndarray, assess_termination_on_minor_steps: bool=False, safety_factor: float=0.8, maximum_factor_increase: float=4.0, minimum_factor_increase: float=0.1, throw_exception_if_minimum_step_exceeded: bool=True) -> IntegratorSettings:
    ...

def standard_cartesian_state_element_blocks(number_of_rows: int, number_of_columns: int) -> list[tuple[int, int, int, int]]:
    ...

def standard_rotational_state_element_blocks(number_of_rows: int, number_of_columns: int) -> list[tuple[int, int, int, int]]:
    ...

def step_size_control_blockwise_matrix_tolerance(block_indices: list[tuple[int, int, int, int]], relative_error_tolerance: numpy.ndarray, absolute_error_tolerance: numpy.ndarray, safety_factor: float=0.8, minimum_factor_increase: float=0.1, maximum_factor_increase: float=4.0) -> IntegratorStepSizeControlSettings:
    ...

def step_size_control_blockwise_scalar_tolerance(block_indices: list[tuple[int, int, int, int]], relative_error_tolerance: float, absolute_error_tolerance: float, safety_factor: float=0.8, minimum_factor_increase: float=0.1, maximum_factor_increase: float=4.0) -> IntegratorStepSizeControlSettings:
    ...

def step_size_control_custom_blockwise_matrix_tolerance(block_indices_function: typing.Callable[[int, int], list[tuple[int, int, int, int]]], relative_error_tolerance: numpy.ndarray, absolute_error_tolerance: numpy.ndarray, safety_factor: float=0.8, minimum_factor_increase: float=0.1, maximum_factor_increase: float=4.0) -> IntegratorStepSizeControlSettings:
    ...

def step_size_control_custom_blockwise_scalar_tolerance(block_indices_function: typing.Callable[[int, int], list[tuple[int, int, int, int]]], relative_error_tolerance: float, absolute_error_tolerance: float, safety_factor: float=0.8, minimum_factor_increase: float=0.1, maximum_factor_increase: float=4.0) -> IntegratorStepSizeControlSettings:
    ...

def step_size_control_elementwise_matrix_tolerance(relative_error_tolerance: numpy.ndarray, absolute_error_tolerance: numpy.ndarray, safety_factor: float=0.8, minimum_factor_increase: float=0.1, maximum_factor_increase: float=4.0) -> IntegratorStepSizeControlSettings:
    ...

def step_size_control_elementwise_scalar_tolerance(relative_error_tolerance: float, absolute_error_tolerance: float, safety_factor: float=0.8, minimum_factor_increase: float=0.1, maximum_factor_increase: float=4.0) -> IntegratorStepSizeControlSettings:
    ...

def step_size_validation(minimum_step: float, maximum_step: float, minimum_step_size_handling: MinimumIntegrationTimeStepHandling=..., accept_infinity_step: bool=False, accept_nan_step: bool=False) -> IntegratorStepSizeValidationSettings:
    ...
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