import numpy
import typing
__all__ = ['AvailableLookupScheme', 'BoundaryInterpolationType', 'InterpolatorGenerationSettings', 'InterpolatorSettings', 'LagrangeInterpolatorBoundaryHandling', 'LagrangeInterpolatorSettings', 'OneDimensionalInterpolatorMatrix', 'OneDimensionalInterpolatorScalar', 'OneDimensionalInterpolatorVector', 'binary_search', 'create_one_dimensional_matrix_interpolator', 'create_one_dimensional_scalar_interpolator', 'create_one_dimensional_vector_interpolator', 'cubic_spline_interpolation', 'extrapolate_at_boundary', 'extrapolate_at_boundary_with_warning', 'hermite_interpolation', 'hermite_spline_interpolation', 'hunting_algorithm', 'lagrange_cubic_spline_boundary_interpolation', 'lagrange_interpolation', 'lagrange_no_boundary_interpolation', 'linear_interpolation', 'piecewise_constant_interpolation', 'throw_exception_at_boundary', 'use_boundary_value', 'use_boundary_value_with_warning']

class AvailableLookupScheme:
    """Enumeration of types of behaviour to be used beyond the edges of the interpolation domain.
	
	When the interpolation is performed, the interpolator scheme will typically start by finding the nearest neighbor of
	the requested value of the independent variable :math:`t` in the data set :math:`[t_{0}..t_{N}]`.
	The choice of lookup scheme can have a significant influence on computational efficiency for large data sets and/or simple
	interpolation algorithms
	
	
	:member hunting_algorithm:
	:member binary_search:
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, AvailableLookupScheme]]
    binary_search: typing.ClassVar[AvailableLookupScheme]
    hunting_algorithm: typing.ClassVar[AvailableLookupScheme]

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

class BoundaryInterpolationType:
    """Enumeration of types of behaviour to be used beyond the edges of the interpolation domain.
	
	Enumeration of types of behaviour to be used beyond the edges of the interpolation domain. For independent variable
	data in the range :math:`[t_{0}..t_{N}]`, this enum is used to define the behaviour of the interpolator at
	:math:`t<t_{0}` and :math:`t>t_{N}
	
	
	:member throw_exception_at_boundary:
	:member use_boundary_value:
	:member use_boundary_value_with_warning:
	:member extrapolate_at_boundary:
	:member extrapolate_at_boundary_with_warning:
	
	
	
	
	
	
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, BoundaryInterpolationType]]
    extrapolate_at_boundary: typing.ClassVar[BoundaryInterpolationType]
    extrapolate_at_boundary_with_warning: typing.ClassVar[BoundaryInterpolationType]
    throw_exception_at_boundary: typing.ClassVar[BoundaryInterpolationType]
    use_boundary_value: typing.ClassVar[BoundaryInterpolationType]
    use_boundary_value_with_warning: typing.ClassVar[BoundaryInterpolationType]

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

class InterpolatorGenerationSettings:
    """
		"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class InterpolatorSettings:
    """Base class to define settings for an interpolator.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

class LagrangeInterpolatorBoundaryHandling:
    """Enumeration of types of behaviour to be used close to the edges of the interpolation domain, for the Lagrange interpolator.
	
	Enumeration of types of behaviour to be used close to the edges of the interpolation domain, for the Lagrange interpolator.
	As explained for :func:`lagrange_interpolation`, the algorithm for the Lagrange interpolation breaks down at the edges of
	the interpolation domain. This enum provides the available options a user has to deal with this.
	
	
	:member lagrange_cubic_spline_boundary_interpolation:
	:member lagrange_no_boundary_interpolation:
	
	
	
	
	
	
	
	"""
    __members__: typing.ClassVar[dict[str, LagrangeInterpolatorBoundaryHandling]]
    lagrange_cubic_spline_boundary_interpolation: typing.ClassVar[LagrangeInterpolatorBoundaryHandling]
    lagrange_no_boundary_interpolation: typing.ClassVar[LagrangeInterpolatorBoundaryHandling]

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

class LagrangeInterpolatorSettings(InterpolatorSettings):
    """:class:`InterpolatorSettings`-derived class to define settings for a Lagrange interpolator.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, interpolate_order: int, use_long_double_time_step: bool=0, selected_lookup_scheme: AvailableLookupScheme=..., lagrange_boundary_handling: LagrangeInterpolatorBoundaryHandling=..., boundary_handling: BoundaryInterpolationType=...) -> None:
        ...

class OneDimensionalInterpolatorMatrix:
    """Object that performs interpolation for matrix independent, and matrix dependent variables.
	
	Object that performs interpolation for matrix independent, and matrix dependent variables. This object is
	not created manually, but is set up using the :func:`create_one_dimensional_matrix_interpolator` function.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def interpolate(self, independent_variable_value: float) -> numpy.ndarray:
        """
        This function performs the interpolation at the requested independent variable value.
        
        	:param independent_variable_value:
        		Value of independent variable at which the interpolation is to be performed.
        
        	:return:
        		Interpolated dependent variable value, using implemented algorithm at requested independent variable value
        """

class OneDimensionalInterpolatorScalar:
    """Object that performs interpolation for scalar independent, and scalar dependent variables.
	
	Object that performs interpolation for scalar independent, and scalar dependent variables. This object is
	not created manually, but is set up using the :func:`create_one_dimensional_scalar_interpolator` function.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def interpolate(self, independent_variable_value: float) -> float:
        """
        This function performs the interpolation at the requested independent variable value.
        
        	:param independent_variable_value:
        		Value of independent variable at which the interpolation is to bse performed.
        
        	:return:
        		Interpolated dependent variable value, using implemented algorithm at requested independent variable value
        """

class OneDimensionalInterpolatorVector:
    """Object that performs interpolation for vector independent, and vector dependent variables.
	
	Object that performs interpolation for vector independent, and vector dependent variables. This object is
	not created manually, but is set up using the :func:`create_one_dimensional_vector_interpolator` function.
	"""

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def interpolate(self, independent_variable_value: float) -> numpy.ndarray:
        """
        This function performs the interpolation at the requested independent variable value.
        
        	:param independent_variable_value:
        		Value of independent variable at which the interpolation is to be performed.
        
        	:return:
        		Interpolated dependent variable value, using implemented algorithm at requested independent variable value
        """

def create_one_dimensional_matrix_interpolator(data_to_interpolate: dict[float, numpy.ndarray], interpolator_settings: InterpolatorSettings, data_first_derivatives: list[numpy.ndarray]=[]) -> ...:
    """Function to create an interpolator for matrix dependent variables.
	
	As :func:`create_one_dimensional_scalar_interpolator`, but with matrices (2-dimensional arrays) as dependent variables
	
	
	:param data_to_interpolate:
			Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
	:param interpolator_settings:
			Settings that define the type of interpolator that is to be used
	:param data_first_derivatives:
			Key-value container with pairs of independent variables (key) and first derivative dependent variables w.r.t. independent variable (value) from which the interpolation is to be performed. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator
	:return:
			Interpolator object
	"""

def create_one_dimensional_scalar_interpolator(data_to_interpolate: dict[float, float], interpolator_settings: InterpolatorSettings, data_first_derivatives: list[float]=[]) -> ...:
    """Function to create an interpolator for scalar dependent variables.
	
	Function to create an interpolator for scalar dependent variables, with a single independent
	variable. This function takes the interpolator settings, and the data that is to be interpolated,
	as input to create the object that can perform the actual interpolation
	
	
	:param data_to_interpolate:
			Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
	:param interpolator_settings:
			Settings that define the type of interpolator that is to be used
	:param data_first_derivatives:
			Key-value container with pairs of independent variables (key) and first derivative dependent variables w.r.t. independent variable (value) from which the interpolation is to be performed. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator
	:return:
			Interpolator object
	"""

def create_one_dimensional_vector_interpolator(data_to_interpolate: dict[float, numpy.ndarray], interpolator_settings: InterpolatorSettings, data_first_derivatives: list[numpy.ndarray]=[]) -> ...:
    """Function to create an interpolator for vector dependent variables.
	
	As :func:`create_one_dimensional_scalar_interpolator`, but with vectors as dependent variables
	
	
	:param data_to_interpolate:
			Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
	:param interpolator_settings:
			Settings that define the type of interpolator that is to be used
	:param data_first_derivatives:
			Key-value container with pairs of independent variables (key) and first derivative dependent variables w.r.t. independent variable (value) from which the interpolation is to be performed. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator).
	:return:
			Interpolator object
	"""

def cubic_spline_interpolation(lookup_scheme: AvailableLookupScheme=..., boundary_interpolation: BoundaryInterpolationType=...) -> InterpolatorSettings:
    """Function to create settings for cubic spline interpolation.
	
	Function to create settings for cubic spline interpolation, where the interpolator
	defines a cubic curve polynomial curve between each two subsequent intervals of the
	independent variable input data. The curve has continuous value, first derivative and
	second derivate between subsequent intervals. As boundary condition, the spline has
	a zero second derivative imposed at the upper and lower boundaries of the interpolation
	domain.
	
	
	:param lookup_scheme:
			Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
	:param boundary_interpolation:
			Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
	:return:
			Cubic spline settings object
	"""

def hermite_interpolation(lookup_scheme: AvailableLookupScheme=..., boundary_interpolation: BoundaryInterpolationType=...) -> InterpolatorSettings:
    ...

def hermite_spline_interpolation(lookup_scheme: AvailableLookupScheme=..., boundary_interpolation: BoundaryInterpolationType=...) -> InterpolatorSettings:
    """Function to create settings for cubic Hermite spline interpolation.
	
	Function to create settings for piecewise cubic Hermite spline interpolation. To use this
	interpolator, a key-value container of values, and a key-value container of first derivatives,
	must be provided to the function creating an interpolator (:func:`create_one_dimensional_scalar_interpolator`,
	:func:`create_one_dimensional_vector_interpolator`,  :func:`create_one_dimensional_matrix_interpolator`). The resulting
	spline uses the value and first derivatives (four piece of information for each interval) at two subsequent nodes to construct
	a cubic polynomial between each two subsequent nodes. The resulting spline has constant values and first derivatives
	
	
	:param lookup_scheme:
			Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
	:param boundary_interpolation:
			Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
	:return:
			Hermite spline interpolation settings object
	"""

def lagrange_interpolation(number_of_points: int, lookup_scheme: AvailableLookupScheme=..., boundary_interpolation: BoundaryInterpolationType=..., lagrange_boundary_handling: LagrangeInterpolatorBoundaryHandling=...) -> InterpolatorSettings:
    """Function to create settings for cubic Lagrange interpolation.
	
	Function to create settings for piecewise cubic Lagrange interpolation.  This is typically the interpolator of highest accuracy that is available.
	The Lagrange interpolator uses :math:`m` consecutive points (input to this function) from the input independent variables :math:`[t_{0}...t_{N}]`
	to create the polynomial of order :math:`m-1` that interpolates these points. From here on, we assume :math:`m` is even.
	The algorithm that is used (see `here <https://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html>`_ for mathematical details
	on interpolating Lagrange polynomials) works as follows:
	
	* The nearest lower neighbor of the data point :math:`t` at which the state :math:`\\mathbf{x}` is to be interpolated is determined, and denoted :math:`t_{i}`.
	* An interpolating Lagrange polynomial is constructed from the consecutive data points :math:`[t_{i-(m/2-1)}...t_{i+m}]`
	* This resulting interpolating polynomial is *only* used in the interval :math:`[t_{i}..t_{i+1}]`, to prevent `Runge's phenomenon <https://en.wikipedia.org/wiki/Runge%27s_phenomenon>`_.
	
	For instance, if :math:`m=8` we use a :math:`7^{th}` order polynomial that interpolates a contiguous set of
	8 data points out of the full data set. Normally, the interpolating polynomial is only used between the
	:math:`4^{th}` and :math:`5^{th}` data point, where it will typically be of good accuracy. Consequently,
	a separate interpolating polynomial (using data over a span of :math:`m` consecutive points) is used for
	each single interval :math:`[t_{i}..t_{i+1}]` (with the exception of the boundaries, see below).
	
	.. warning:: Issues can occur if the data point :math:`t` at which the interpolation is to be
				 performed is close to :math:`t_{0}` or :math:`t_{N}`. In those case, there is not sufficient
				 data to construct the interpolating polynomial *and* to only use this interpolating polynomial
				 between the middle two data points that were used to it. In these cases, the user has a number of
				 options (all defined by an entry of the :class:`LagrangeInterpolatorBoundaryHandling` variable,
				 used as input to this function). In short, interpolation between the first and last :math:`m/2`
				 data points will lead to degraded results, warnings, or termination.
	
	
	:param number_of_points:
			Number of consecutive data points that are used to construct a single interpolating polynomial.
	:param lookup_scheme:
			Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
	:param boundary_interpolation:
			Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
	:param lagrange_boundary_handling:
			Interpolator behaviour that is to be used at the boundaries of the domain, where the regular algorithm cannot be executed.
	:return:
			Lagrange interpolation settings object
	"""

def linear_interpolation(lookup_scheme: AvailableLookupScheme=..., boundary_interpolation: BoundaryInterpolationType=...) -> InterpolatorSettings:
    """Function to create settings for linear interpolation.
	
	Function to create settings for linear interpolation, where the interpolator
	defines a linear curve between each two subsequent intervals of the
	independent variable input data.
	
	
	:param lookup_scheme:
			Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
	:param boundary_interpolation:
			Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
	:return:
			Linear interpolation settings object
	"""

def piecewise_constant_interpolation(lookup_scheme: AvailableLookupScheme=..., boundary_interpolation: BoundaryInterpolationType=...) -> InterpolatorSettings:
    """Function to create settings for piecewise constant interpolation.
	
	Function to create settings for piecewise constant interpolation. If interpolator
	is to return the value at :math:`t`, and :math:`t_{i}\\le t \\< t_{i+1}`, the interpolator
	returns :math:`\\mathbf{x}_{i}`
	
	
	:param lookup_scheme:
			Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
	:param boundary_interpolation:
			Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
	:return:
			Piecewise constant interpolation settings object
	"""
binary_search: AvailableLookupScheme
extrapolate_at_boundary: BoundaryInterpolationType
extrapolate_at_boundary_with_warning: BoundaryInterpolationType
hunting_algorithm: AvailableLookupScheme
lagrange_cubic_spline_boundary_interpolation: LagrangeInterpolatorBoundaryHandling
lagrange_no_boundary_interpolation: LagrangeInterpolatorBoundaryHandling
throw_exception_at_boundary: BoundaryInterpolationType
use_boundary_value: BoundaryInterpolationType
use_boundary_value_with_warning: BoundaryInterpolationType