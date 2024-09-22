import numpy
import tudatpy.math.interpolators.expose_interpolators
import tudatpy.numerical_simulation.environment.expose_environment
import typing
__all__ = ['AerodynamicCoefficientSettings', 'ConstantAerodynamicCoefficientSettings', 'ControlSurfaceIncrementAerodynamicCoefficientSettings', 'constant', 'custom_aerodynamic_force_and_moment_coefficients', 'custom_aerodynamic_force_coefficients', 'custom_control_surface', 'scaled_by_constant', 'scaled_by_vector', 'scaled_by_vector_function', 'tabulated', 'tabulated_force_only', 'tabulated_force_only_from_files', 'tabulated_from_files', 'tabulated_from_files_control_surface']

class AerodynamicCoefficientSettings:
    """
		"""
    add_force_contribution_to_moments: bool
    moment_reference_point: numpy.ndarray

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def add_single_control_surface(self, control_surface_settings: ..., control_surface_name: str) -> None:
        ...

class ConstantAerodynamicCoefficientSettings(AerodynamicCoefficientSettings):
    """Class for defining model settings from constant aerodynamic coefficients.
	
	`AerodynamicCoefficientSettings` derived class for aerodynamic interface model settings using only constant aerodynamic coefficients.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class ControlSurfaceIncrementAerodynamicCoefficientSettings:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

def constant(reference_area: float, constant_force_coefficient: numpy.ndarray, force_coefficients_frame: tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientFrames=...) -> AerodynamicCoefficientSettings:
    """Factory function for creating aerodynamic interface model settings entirely from constant coefficients.
	
	Factory function for settings object, defining aerodynamic interface model entirely from constant aerodynamic coefficients,
	i.e. coefficients are not a function of any independent variables.
	
	
	:param reference_area:
			Reference area with which aerodynamic forces and moments are non-dimensionalized.
	:param constant_force_coefficient:
			Constant force coefficients.
	:param are_coefficients_in_aerodynamic_frame:
			Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
			(drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).
	
	:param are_coefficients_in_negative_axis_direction:
			Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or
			aerodynamic frame (see arg are_coefficients_in_aerodynamic_frame).
			Note that for drag, side and lift force, the coefficients are typically defined in negative direction.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ConstantAerodynamicCoefficientSettings` class
	"""

def custom_aerodynamic_force_and_moment_coefficients(force_coefficient_function: typing.Callable[[list[float]], numpy.ndarray], moment_coefficient_function: typing.Callable[[list[float]], numpy.ndarray], reference_length: float, reference_area: float, independent_variable_names: list[tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientsIndependentVariables], force_coefficients_frame: tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientFrames=..., moment_coefficients_frame: tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientFrames=..., moment_reference_point: numpy.ndarray=...) -> AerodynamicCoefficientSettings:
    """Factory function for creating aerodynamic interface model settings from custom coefficients.
	
	Factory function for settings object, defining aerodynamic interface model via a custom force coefficient function
	(function of independent variable).
	
	
	:param force_coefficient_function:
			Function that is defining the aerodynamic force coefficients as function of an independent variable (see arg independent_variable_names).
	:param moment_coefficient_function:
			Function that is defining the aerodynamic moment coefficients as function of an independent variable (see arg independent_variable_names).
	:param reference_area:
			Reference area with which aerodynamic forces and moments are non-dimensionalized.
	:param reference_length:
			Reference length with which aerodynamic moments are non-dimensionalized.
	:param moment_reference_point:
			Reference point in the body-fixed frame w.r.t. which the moments are provided.
	:param independent_variable_name:
			Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
	:param are_coefficients_in_aerodynamic_frame:
			Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
			(drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).
	
	:param are_coefficients_in_negative_axis_direction:
			Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or
			aerodynamic frame (see arg are_coefficients_in_aerodynamic_frame).
			Note that for drag, side and lift force, the coefficients are typically defined in negative direction.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.CustomAerodynamicCoefficientSettings` class
	"""

def custom_aerodynamic_force_coefficients(force_coefficient_function: typing.Callable[[list[float]], numpy.ndarray], reference_area: float, independent_variable_names: list[tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientsIndependentVariables], force_coefficients_frame: tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientFrames=...) -> AerodynamicCoefficientSettings:
    """Factory function for creating aerodynamic interface model settings from custom coefficients.
	
	Factory function for settings object, defining aerodynamic interface model via a custom force coefficient function
	(function of independent variable).
	
	
	:param force_coefficient_function:
			Function that is defining the aerodynamic coefficients as function of an independent variable (see arg independent_variable_names).
	:param reference_area:
			Reference area with which aerodynamic forces and moments are non-dimensionalized.
	:param independent_variable_name:
			Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
	:param are_coefficients_in_aerodynamic_frame:
			Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
			(drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).
	
	:param are_coefficients_in_negative_axis_direction:
			Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or
			aerodynamic frame (see arg are_coefficients_in_aerodynamic_frame).
			Note that for drag, side and lift force, the coefficients are typically defined in negative direction.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.CustomAerodynamicCoefficientSettings` class
	"""

def custom_control_surface(force_and_moment_coefficient_function: typing.Callable[[list[float]], numpy.ndarray], independent_variable_names: list[tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientsIndependentVariables]) -> ...:
    ...

def scaled_by_constant(unscaled_coefficient_settings: AerodynamicCoefficientSettings, force_scaling_constant: float, moment_scaling_constant: float, is_scaling_absolute: bool=False) -> AerodynamicCoefficientSettings:
    """Factory function for creating aerodynamic interface model settings by applying one constant scaling factor/value to all coefficients of an existing model settings object.
	
	Factory function for settings object, defining aerodynamic interface based on scaling the coefficients of an existing model settings object by one constant factor or value.
	Via the ``is_scaling_absolute``
	boolean, the user can apply a constant scaling factor or an absolute value to the resulting force and moment coefficients (for instance for an uncertainty analysis).
	
	
	:param unscaled_coefficient_settings:
			Existing aerodynamic interface model settings object that is used as the base for the scaled settings object.
	:param force_scaling_constant:
			Constant scaling factor to be applied to all aerodynamic force coefficients.
	:param moment_scaling_constant:
			Constant scaling factor to be applied to all aerodynamic moment coefficients.
	:param is_scaling_absolute, default = False:
			Boolean indicating whether aerodynamic coefficient scaling is absolute.
			Setting this boolean to true will add the scaling value to the base value,
			instead of the default behaviour of multiplying the base value by the scaling factor.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ScaledAerodynamicCoefficientInterfaceSettings` class
	"""

def scaled_by_vector(unscaled_coefficient_settings: AerodynamicCoefficientSettings, force_scaling_vector: numpy.ndarray, moment_scaling_vector: numpy.ndarray, is_scaling_absolute: bool=False) -> AerodynamicCoefficientSettings:
    """Factory function for creating aerodynamic interface model settings by applying constant scaling factors/values to the coefficients of an existing model settings object.
	
	Factory function for settings object, defining aerodynamic interface based on scaling the coefficients of an existing model settings object by constant factors or values.
	Via the ``is_scaling_absolute`` boolean, the user can apply one constant scaling factor or an absolute value to each resulting force and moment coefficient (for instance for an uncertainty analysis).
	
	
	:param unscaled_coefficient_settings:
			Existing aerodynamic interface model settings object that is used as the base for the scaled settings object.
	:param force_scaling_vector:
			Constant scaling factors to be applied to each aerodynamic force coefficient.
	:param moment_scaling_vector:
			Constant scaling factors to be applied to each aerodynamic moment coefficient.
	:param is_scaling_absolute, default = False:
			Boolean indicating whether aerodynamic coefficient scaling is absolute.
			Setting this boolean to true will add the scaling value to the base value,
			instead of the default behaviour of multiplying the base value by the scaling factor.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ScaledAerodynamicCoefficientInterfaceSettings` class
	"""

def scaled_by_vector_function(unscaled_coefficient_settings: AerodynamicCoefficientSettings, force_scaling_vector_function: typing.Callable[[float], numpy.ndarray], moment_scaling_vector_function: typing.Callable[[float], numpy.ndarray], is_scaling_absolute: bool=False) -> AerodynamicCoefficientSettings:
    """Factory function for creating aerodynamic interface model settings by applying custom scaling factors/values to the coefficients of an existing model settings object.
	
	Factory function for settings object, defining aerodynamic interface based on scaling the coefficients of an existing model settings object by custom factors or values.
	Via the ``is_scaling_absolute`` boolean, the user can apply the scaling factors or absolute values to each resulting force and moment coefficient (for instance for an uncertainty analysis).
	
	
	:param unscaled_coefficient_settings:
			Existing aerodynamic interface model settings object that is used as the base for the scaled settings object.
	:param force_scaling_vector_function:
			Custom scaling factors to be applied to each aerodynamic force coefficient.
	:param moment_scaling_vector_function:
			Custom scaling factors to be applied to each aerodynamic moment coefficient.
	:param is_scaling_absolute, default = False:
			Boolean indicating whether aerodynamic coefficient scaling is absolute.
			Setting this boolean to true will add the scaling value to the base value,
			instead of the default behaviour of multiplying the base value by the scaling factor.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ScaledAerodynamicCoefficientInterfaceSettings` class
	"""

def tabulated(independent_variables: list[float], force_coefficients: list[numpy.ndarray], moment_coefficients: list[numpy.ndarray], reference_length: float, reference_area: float, independent_variable_name: tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientsIndependentVariables, force_coefficients_frame: tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientFrames=..., moment_coefficients_frame: tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientFrames=..., moment_reference_point: numpy.ndarray=..., interpolator_settings: tudatpy.math.interpolators.expose_interpolators.InterpolatorSettings=None) -> AerodynamicCoefficientSettings:
    """Factory function for creating aerodynamic interface model settings from user-defined, 1-d tabulated coefficients.
	
	Factory function for settings object, defining aerodynamic interface model via user-defined, 1-dimensional, tabulated aerodynamic force and moment coefficients
	(tabulated w.r.t. independent variable).
	
	
	:param independent_variables:
			Values of independent variables at which the coefficients in the input multi vector are defined (size 1).
	:param force_coefficients:
			Values of force coefficients at independent variables defined by independent_variables.
	:param moment_coefficients:
			Values of moment coefficients at independent variables defined by independent_variables.
	:param reference_length:
			Reference length with which aerodynamic moments about x- and z- axes are non-dimensionalized.
	:param reference_area:
			Reference area with which aerodynamic forces and moments are non-dimensionalized.
	:param lateral_reference_length:
			Reference length with which aerodynamic moment about y-axis is non-dimensionalized.
	:param moment_reference_point:
			Point w.r.t. aerodynamic moment is calculated.
	:param independent_variable_name:
			Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
	:param are_coefficients_in_aerodynamic_frame:
			Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
			(drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).
	
	:param are_coefficients_in_negative_axis_direction:
			Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or
			aerodynamic frame (see arg areCoefficientsInAerodynamicFrame).
			Note that for drag, side and lift force, the coefficients are typically defined in negative direction.
	
	:param interpolator_settings:
			Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class (via :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettingsBase` class)
	"""

def tabulated_force_only(independent_variables: list[float], force_coefficients: list[numpy.ndarray], reference_area: float, independent_variable_name: tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientsIndependentVariables, force_coefficients_frame: tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientFrames=..., interpolator_settings: tudatpy.math.interpolators.expose_interpolators.InterpolatorSettings=None) -> AerodynamicCoefficientSettings:
    """Factory function for creating aerodynamic interface model settings from user-defined, 1-d tabulated force coefficients.
	
	Factory function for settings object, defining aerodynamic interface model via user-defined, 1-dimensional, tabulated aerodynamic force coefficients
	(tabulated w.r.t. independent variable).
	
	
	:param independent_variables:
			Values of independent variables at which the coefficients in the input multi vector are defined (size 1)
	:param force_coefficients:
			Values of force coefficients at independent variables defined by independent_variables.
	:param reference_area:
			Reference area with which aerodynamic forces and moments are non-dimensionalized.
	:param independent_variable_name:
			Identifier of the independent variable w.r.t. which the aerodynamic coefficients are defined.
	:param are_coefficients_in_aerodynamic_frame:
			Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
			(drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).
	
	:param are_coefficients_in_negative_axis_direction:
			Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or
			aerodynamic frame (see arg areCoefficientsInAerodynamicFrame).
			Note that for drag, side and lift force, the coefficients are typically defined in negative direction.
	
	:param interpolator_settings:
			Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.
			Pointer to an interpolator settings object where the conditions for interpolation of tabulated inputs are saved.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class
	"""

def tabulated_force_only_from_files(force_coefficient_files: dict[int, str], reference_area: float, independent_variable_names: list[tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientsIndependentVariables], force_coefficients_frame: tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientFrames=..., interpolator_settings: tudatpy.math.interpolators.expose_interpolators.InterpolatorSettings=None) -> AerodynamicCoefficientSettings:
    """Factory function for creating aerodynamic interface model settings from tabulated force coefficients from files.
	
	Factory function for settings object, defining aerodynamic interface model via user-defined, tabulated aerodynamic force coefficients
	(tabulated w.r.t. independent variable), obtained from data files.
	
	
	:param force_coefficient_files:
			Path of the aerodynamic coefficient files corresponding to the force coefficient of the given dict key.
	:param reference_area:
			Reference area with which aerodynamic forces and moments are non-dimensionalized.
	:param independent_variable_names:
			Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
	:param are_coefficients_in_aerodynamic_frame:
			Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
			(drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).
	
	:param are_coefficients_in_negative_axis_direction:
			Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or
			aerodynamic frame (see arg areCoefficientsInAerodynamicFrame).
			Note that for drag, side and lift force, the coefficients are typically defined in negative direction.
	
	:param interpolator_settings:
			Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class
	"""

def tabulated_from_files(force_coefficient_files: dict[int, str], moment_coefficient_files: dict[int, str], reference_length: float, reference_area: float, independent_variable_names: list[tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientsIndependentVariables], force_coefficients_frame: tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientFrames=..., moment_coefficients_frame: tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientFrames=..., moment_reference_point: numpy.ndarray=..., interpolator_settings: tudatpy.math.interpolators.expose_interpolators.InterpolatorSettings=None) -> AerodynamicCoefficientSettings:
    """Factory function for creating aerodynamic interface model settings from tabulated coefficients from files.
	
	Factory function for settings object, defining aerodynamic interface model via user-defined, tabulated aerodynamic force and moment coefficients
	(tabulated w.r.t. independent variable), obtained from data files.
	
	
	:param force_coefficient_files:
			Path of the aerodynamic coefficient files corresponding to the force coefficient of the given dict key.
	:param moment_coefficient_files:
			Path of the aerodynamic coefficient files corresponding to the moment coefficient of the given dict key.
	:param reference_length:
			Reference length with which aerodynamic moments about x- and z- axes are non-dimensionalized.
	:param reference_area:
			Reference area with which aerodynamic forces and moments are non-dimensionalized.
	:param lateral_reference_length:
			Reference length with which aerodynamic moment about y-axis is non-dimensionalized.
	:param moment_reference_point:
			Point w.r.t. aerodynamic moment is calculated.
	:param independent_variable_names:
			Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
	:param are_coefficients_in_aerodynamic_frame:
			Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
			(drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).
	
	:param are_coefficients_in_negative_axis_direction:
			Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or
			aerodynamic frame (see arg areCoefficientsInAerodynamicFrame).
			Note that for drag, side and lift force, the coefficients are typically defined in negative direction.
	
	:param interpolator_settings:
			Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.
	
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class
	"""

def tabulated_from_files_control_surface(force_coefficient_files: dict[int, str], moment_coefficient_files: dict[int, str], independent_variable_names: list[tudatpy.numerical_simulation.environment.expose_environment.AerodynamicCoefficientsIndependentVariables]) -> ControlSurfaceIncrementAerodynamicCoefficientSettings:
    ...