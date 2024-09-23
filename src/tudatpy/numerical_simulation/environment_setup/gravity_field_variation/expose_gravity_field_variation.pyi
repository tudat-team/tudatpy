import numpy
import tudatpy.math.interpolators.expose_interpolators
import typing
__all__ = ['BasicSolidBodyGravityFieldVariationSettings', 'BodyDeformationTypes', 'GravityFieldVariationSettings', 'basic_solid_body', 'periodic', 'polynomial', 'single_period_periodic', 'single_power_polynomial', 'solid_body_tide', 'solid_body_tide_complex_k', 'solid_body_tide_degree_order_variable_complex_k', 'solid_body_tide_degree_order_variable_k', 'solid_body_tide_degree_variable_complex_k', 'solid_body_tide_degree_variable_k', 'solid_multi_body_tide_degree_order_variable_k', 'tabulated', 'tabulated_deformation']

class BasicSolidBodyGravityFieldVariationSettings(GravityFieldVariationSettings):
    """Class for providing settings for solid body tidal gravity field variations, derived from GravityFieldVariationSettings.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class BodyDeformationTypes:
    """<no_doc>
	
	Members:
	
	basic_solid_body
	
	tabulated_deformation
	"""
    __members__: typing.ClassVar[dict[str, BodyDeformationTypes]]
    basic_solid_body: typing.ClassVar[BodyDeformationTypes]
    tabulated_deformation: typing.ClassVar[BodyDeformationTypes]

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

class GravityFieldVariationSettings:
    """Base class for providing settings for gravity field variations.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

def periodic(cosine_coefficient_amplitudes_cosine_time: list[numpy.ndarray], cosine_coefficient_amplitudes_sine_time: list[numpy.ndarray], sine_coefficient_amplitudes_cosine_time: list[numpy.ndarray], sine_coefficient_amplitudes_sine_time: list[numpy.ndarray], frequencies: list[float], reference_epoch: float, minimum_degree: int=2, minimum_order: int=0) -> GravityFieldVariationSettings:
    ...

def polynomial(cosine_amplitudes_per_power: dict[int, numpy.ndarray], sine_amplitudes_per_power: dict[int, numpy.ndarray], reference_epoch: float, minimum_degree: int=2, minimum_order: int=0) -> GravityFieldVariationSettings:
    ...

def single_period_periodic(cosine_coefficient_amplitude_cosine_time: numpy.ndarray, cosine_coefficient_amplitude_sine_time: numpy.ndarray, sine_coefficient_amplitude_cosine_time: numpy.ndarray, sine_coefficient_amplitude_sine_time: numpy.ndarray, frequency: float, reference_epoch: float, minimum_degree: int=2, minimum_order: int=0) -> GravityFieldVariationSettings:
    ...

def single_power_polynomial(cosine_amplitudes: numpy.ndarray, sine_amplitudes: numpy.ndarray, polynomial_power: int, reference_epoch: float, minimum_degree: int=2, minimum_order: int=0) -> GravityFieldVariationSettings:
    ...

def solid_body_tide(tide_raising_body: str, love_number: float, degree: int) -> GravityFieldVariationSettings:
    """Factory function for creating solid body tides.
	
	Factory function for creating solid body tides, using a single real Love number at a single degree (e.g. :math:`k_{2}`, :math:`k_{3}`, etc.). This function evaluates Eq. (6.6) from the IERS Conventions 2010, with real :math:`k_{l}=k_{lm}`, a single value of :math:`l` and a single tide-raising body :math:`j`.
	
	
	:param tide_raising_body:
			Name of body raising the tide.
	:param love_number:
			Constant real Love number to use for body undergoing deformation, at the spherical harmonic degree defined by 'degree' input.
	:param degree:
			Degree of the spherical harmonic gravity field, and associated Love number, that is to be considered.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
	"""

def solid_body_tide_complex_k(tide_raising_body: str, love_number: complex, degree: int) -> GravityFieldVariationSettings:
    """Factory function for creating solid body tides.
	
	As :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide`, but with complex value for the Love number.
	
	
	:param tide_raising_body:
			Name of body raising the tide.
	:param love_number:
			Constant real Love number to use for body undergoing deformation, at the spherical harmonic degree defined by 'degree' input.
	:param degree:
			Degree of the spherical harmonic gravity field, and associated Love number, that is to be considered.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
	"""

@typing.overload
def solid_body_tide_degree_order_variable_complex_k(tide_raising_body: str, love_number_per_degree_and_order: dict[int, list[complex]]) -> GravityFieldVariationSettings:
    """Factory function for creating solid body tides.
	
	As :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k`, but with complex values for the Love number.
	
	
	:param tide_raising_body:
			Name of body raising the tide.
	:param love_number_per_degree:
			Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the list of Love numbers :math:`k_{lm}` at this degree. Note that, for Love numbers at degree :math:`l`, the associated list should contain :math:`l+1` entries, representing the Love numbers (in order) :math:`k_{l0}`, :math:`k_{l1}`...:math:`k_{ll}`.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
	"""

@typing.overload
def solid_body_tide_degree_order_variable_complex_k(tide_raising_body: str, love_number_per_degree_and_order: dict[int, list[complex]]) -> GravityFieldVariationSettings:
    """Factory function for creating solid body tides.
	
	As :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k`, but with complex values for the Love number.
	
	
	:param tide_raising_body:
			Name of body raising the tide.
	:param love_number_per_degree:
			Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the list of Love numbers :math:`k_{lm}` at this degree. Note that, for Love numbers at degree :math:`l`, the associated list should contain :math:`l+1` entries, representing the Love numbers (in order) :math:`k_{l0}`, :math:`k_{l1}`...:math:`k_{ll}`.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
	"""

def solid_body_tide_degree_order_variable_k(tide_raising_body: str, love_number_per_degree_and_order: dict[int, list[float]]) -> GravityFieldVariationSettings:
    """Factory function for creating solid body tides.
	
	Factory function for creating solid body tides, using a set of real, separate, Love numbers at any number of degrees and orders (e.g. :math:`k_{20}`, :math:`k_{21}`, :math:`k_{22}`, :math:`k_{30}`, etc.).  This function evaluates Eq. (6.6) from the IERS Conventions 2010, with a set of real values :math:`k_{lm}`, at a set of values of :math:`l` and a single tide-raising body :math:`j`.
	
	
	:param tide_raising_body:
			Name of body raising the tide.
	:param love_number_per_degree_and_order:
			Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the list of Love numbers :math:`k_{lm}` at this degree. Note that, for Love numbers at degree :math:`l`, the associated list should contain :math:`l+1` entries, representing the Love numbers (in order) :math:`k_{l0}`, :math:`k_{l1}`...:math:`k_{ll}`.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
	"""

def solid_body_tide_degree_variable_complex_k(tide_raising_body: str, love_number_per_degree: dict[int, complex]) -> GravityFieldVariationSettings:
    """Factory function for creating solid body tides.
	
	As :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide_degree_variable_k`, but with complex values for the Love numbers.
	
	
	:param tide_raising_body:
			Name of body raising the tide.
	:param love_number_per_degree:
			Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the Love number :math:`k_{l}` itself.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
	"""

def solid_body_tide_degree_variable_k(tide_raising_body: str, love_number_per_degree: dict[int, float]) -> GravityFieldVariationSettings:
    """Factory function for creating solid body tides.
	
	Factory function for creating solid body tides, using a set of real, separate, Love numbers at any number of degrees (e.g. :math:`k_{2}`, :math:`k_{3}`, etc.). This output of this function is effectively identical to a list of outputs to :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide`, with differing degrees and associated Love numbers.  This function evaluates Eq. (6.6) from the IERS Conventions 2010, with a set of real values :math:`k_{l}=k_{lm}`, at a set of values of :math:`l` and a single tide-raising body :math:`j`.
	
	
	:param tide_raising_body:
			Name of body raising the tide.
	:param love_number_per_degree:
			Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the Love number :math:`k_{l}` itself
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
	"""

def solid_multi_body_tide_degree_order_variable_k(tide_raising_bodies: list[str], love_number_per_degree_and_order: dict[int, list[float]]) -> GravityFieldVariationSettings:
    ...

def tabulated(cosine_variations_table: dict[float, numpy.ndarray], sine_variations_table: dict[float, numpy.ndarray], minimum_degree: int, minimum_order: int, interpolation_settings: tudatpy.math.interpolators.expose_interpolators.InterpolatorSettings) -> GravityFieldVariationSettings:
    ...
basic_solid_body: BodyDeformationTypes
tabulated_deformation: BodyDeformationTypes