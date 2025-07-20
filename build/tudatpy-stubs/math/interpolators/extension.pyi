import numpy
import pybind11_stubgen.typing_ext
import typing
__all__ = ['AvailableLookupScheme', 'BoundaryInterpolationType', 'InterpolatorGenerationSettings', 'InterpolatorGenerationSettingsTimeObject', 'InterpolatorSettings', 'LagrangeInterpolatorBoundaryHandling', 'LagrangeInterpolatorSettings', 'OneDimensionalInterpolatorMatrix', 'OneDimensionalInterpolatorMatrixTimeObject', 'OneDimensionalInterpolatorScalar', 'OneDimensionalInterpolatorScalarTimeObject', 'OneDimensionalInterpolatorVector', 'OneDimensionalInterpolatorVectorTimeObject', 'binary_search', 'create_one_dimensional_matrix_interpolator', 'create_one_dimensional_matrix_interpolator_time_object', 'create_one_dimensional_scalar_interpolator', 'create_one_dimensional_scalar_interpolator_time_object', 'create_one_dimensional_vector_interpolator', 'create_one_dimensional_vector_interpolator_time_object', 'cubic_spline_interpolation', 'extrapolate_at_boundary', 'extrapolate_at_boundary_with_warning', 'hermite_interpolation', 'hermite_spline_interpolation', 'hunting_algorithm', 'interpolator_generation_settings', 'interpolator_generation_settings_time_object', 'lagrange_cubic_spline_boundary_interpolation', 'lagrange_interpolation', 'lagrange_no_boundary_interpolation', 'linear_interpolation', 'piecewise_constant_interpolation', 'throw_exception_at_boundary', 'use_boundary_value', 'use_boundary_value_with_warning']

class AvailableLookupScheme:
    """Enumeration of types of behaviour to be used beyond the edges of the interpolation domain.
    
    When the interpolation is performed, the interpolator scheme will typically start by finding the nearest neighbor of
    the requested value of the independent variable :math:`t` in the data set :math:`[t_{0}...t_{N}]`.
    The choice of lookup scheme can have a significant influence on computational efficiency for large data sets and/or simple
    interpolation algorithms
    
    
    
    
    
          
    
    Members:
    
      hunting_algorithm : 
    
    With this option, the interpolator 'remembers' which value of :math:`t_{i}` was the nearest neighbor during the previous call to the interpolate function, and starts looking at/near this entry of the data set :math:`[t_{i}]` to find the nearest neighbor.
    
          
    
      binary_search : 
    
    With this option, the algorithm uses a binary search algorithm to find the nearest neighbor, initially starting with the full data range :math:`[t_{0}...t_{N}]`.
    
          """
    __members__: typing.ClassVar[dict[str, AvailableLookupScheme]]
    binary_search: typing.ClassVar[AvailableLookupScheme]
    hunting_algorithm: typing.ClassVar[AvailableLookupScheme]

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
    data in the range :math:`[t_{0}...t_{N}]`, this enum is used to define the behaviour of the interpolator at
    :math:`t<t_{0}` and :math:`t>t_{N}`
    
    
    
    
    
          
    
    Members:
    
      throw_exception_at_boundary : 
    
    The program will terminate and throw a :class:`~tudatpy.exceptions.InterpolationOutOfBoundsError` error when the interpolator is interrogated beyond the range :math:`[t_{0}...t_{N}]`.
    
          
    
      use_boundary_value : 
    
    The value :math:`\\mathbf{x}_{0}` is returned for :math:`t<t_{0}` (and :math:`\\mathbf{x}_{N}` if :math:`t>t_{N}`).
    
          
    
      use_boundary_value_with_warning : 
    
    Same as ``use_boundary_value``, but a warning is printed to the terminal.
    
          
    
      extrapolate_at_boundary : 
    
    The interpolation scheme is extended beyond the range :math:`t_{0}...t_{N}` without any warning. That is, the mathematical equation used to compute the value of :math:`x` in the range :math:`[t_{0}...t_{1}]` is used without any checks for :math:`t<t_{0}`  (and equivalently for :math:`t>t_{N}`). Warning, using this setting can result in divergent/unrealistic behaviour.
    
          
    
      extrapolate_at_boundary_with_warning : 
    
    Same as ``extrapolate_at_boundary``, but a warning is printed to the terminal.
    
          """
    __members__: typing.ClassVar[dict[str, BoundaryInterpolationType]]
    extrapolate_at_boundary: typing.ClassVar[BoundaryInterpolationType]
    extrapolate_at_boundary_with_warning: typing.ClassVar[BoundaryInterpolationType]
    throw_exception_at_boundary: typing.ClassVar[BoundaryInterpolationType]
    use_boundary_value: typing.ClassVar[BoundaryInterpolationType]
    use_boundary_value_with_warning: typing.ClassVar[BoundaryInterpolationType]

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
    """No documentation found."""

class InterpolatorGenerationSettingsTimeObject:
    """No documentation found."""

class InterpolatorSettings:
    """Base class to define settings for an interpolator."""

class LagrangeInterpolatorBoundaryHandling:
    """Enumeration of types of behaviour to be used close to the edges of the interpolation domain, for the Lagrange interpolator.
    
    Enumeration of types of behaviour to be used close to the edges of the interpolation domain, for the Lagrange interpolator.
    As explained for :func:`lagrange_interpolation`, the algorithm for the Lagrange interpolation breaks down at the edges of
    the interpolation domain. This enum provides the available options a user has to deal with this.
    
    
    
    
    
          
    
    Members:
    
      lagrange_cubic_spline_boundary_interpolation : 
    
    A cubic-spline interpolator is created from the first and last :math:`\\max(m/2-1,4)` data points of the full data set, and these cubic spline interpolators are used when an interpolation at :math:`t<t_{(m/2-1)}` or :math:`t<t_{N-(m/2)}` is called.
    
          
    
      lagrange_no_boundary_interpolation : 
    
    The program will terminate and throw a :class:`~tudatpy.exceptions.LagrangeInterpolationOutOfBoundsError` when the Lagrange interpolator is interrogated beyond its valid range.
    
          """
    __members__: typing.ClassVar[dict[str, LagrangeInterpolatorBoundaryHandling]]
    lagrange_cubic_spline_boundary_interpolation: typing.ClassVar[LagrangeInterpolatorBoundaryHandling]
    lagrange_no_boundary_interpolation: typing.ClassVar[LagrangeInterpolatorBoundaryHandling]

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
    """:class:`InterpolatorSettings`-derived class to define settings for a Lagrange interpolator."""

    def __init__(self, interpolate_order: int, use_long_double_time_step: bool=0, selected_lookup_scheme: AvailableLookupScheme=..., lagrange_boundary_handling: LagrangeInterpolatorBoundaryHandling=..., boundary_handling: BoundaryInterpolationType=...) -> None:
        ...

class OneDimensionalInterpolatorMatrix:
    """Object that performs interpolation for matrix dependent variables.
    
    Object that performs interpolation for matrix dependent variables and float independent variables. This object is
    not created manually, but is set up using the :func:`create_one_dimensional_vector_interpolator` function."""

    def interpolate(self, independent_variable_value: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 This function performs the interpolation at the requested independent variable value.
        
        
                 Parameters
                 ----------
                 independent_variable_value : float
                     Value of independent variable at which the interpolation is to be performed.
        
                 Returns
                 -------
                 np.array
                     Interpolated dependent variable value, using implemented algorithm at requested independent variable value
        """

    @property
    def dependent_values(self) -> list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]:
        """
                 Returns the dependent variable values used by the interpolator.
        
                 Returns the dependent variable values used by the interpolator. This is a read-only property.
        
                 Returns
                 -------
                 list[np.ndarray]
                     Dependent variable values used by the interpolator
        """

    @property
    def independent_values(self) -> list[float]:
        """
                 Returns the independent variable values used by the interpolator.
        
                 Returns the independent variable values used by the interpolator. This is a read-only property.
        
                 Returns
                 -------
                 list[float]
                     Independent variable values used by the interpolator
        """

class OneDimensionalInterpolatorMatrixTimeObject:
    """Same as :func:`~OneDimensionalInterpolatorMatrix`, but using the high-resolution :func:`~Time` type used as independent variable for interpolation; created using :func:`~create_one_dimensional_matrix_interpolator_time_object`"""

    def interpolate(self, independent_variable_value: ...) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                This function performs the interpolation at the requested independent variable value.
        
                Parameters
                ----------
                independent_variable_value : float
                    Value of independent variable at which the interpolation is to be performed.
                Returns
                -------
                np.array
                    Interpolated dependent variable value, using implemented algorithm at requested independent variable value
        """

class OneDimensionalInterpolatorScalar:
    """Object that performs interpolation for scalar dependent variables .
    
    Object that performs interpolation for scalar dependent variables and float independent variables. This object is
    not created manually, but is set up using the :func:`create_one_dimensional_scalar_interpolator` function."""

    def interpolate(self, independent_variable_value: float) -> float:
        """
                 This function performs the interpolation at the requested independent variable value.
        
        
                 Parameters
                 ----------
                 independent_variable_value : float
                     Value of independent variable at which the interpolation is to be performed.
        
                 Returns
                 -------
                 float
                     Interpolated dependent variable value, using implemented algorithm at requested independent variable value
        """

    @property
    def dependent_values(self) -> list[float]:
        """
                 Returns the dependent variable values used by the interpolator.
        
                 Returns the dependent variable values used by the interpolator. This is a read-only property.
        
                 Returns
                 -------
        
                 list[float]
                     Dependent variable values used by the interpolator
        """

    @property
    def independent_values(self) -> list[float]:
        """
                 Returns the independent variable values used by the interpolator.
        
                 Returns the independent variable values used by the interpolator. This is a read-only property.
        
                 Returns
                 -------
                 list[float]
                     Independent variable values used by the interpolator
        """

class OneDimensionalInterpolatorScalarTimeObject:
    """Same as :func:`~OneDimensionalInterpolatorScalar`, but using the high-resolution :func:`~Time` type used as independent variable for interpolation; created using :func:`~create_one_dimensional_scalar_interpolator_time_object`"""

    def interpolate(self, independent_variable_value: ...) -> float:
        """
                This function performs the interpolation at the requested independent variable value.
        
                Parameters
                ----------
                independent_variable_value : Time
                    Value of independent variable at which the interpolation is to be performed.
                Returns
                -------
                float
                    Interpolated dependent variable value, using implemented algorithm at requested independent variable value
        """

class OneDimensionalInterpolatorVector:
    """Object that performs interpolation for vector dependent variables.
    
    Object that performs interpolation for vector dependent variables and float independent variable. This object is
    not created manually, but is set up using the :func:`create_one_dimensional_vector_interpolator` function."""

    def interpolate(self, independent_variable_value: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
                 This function performs the interpolation at the requested independent variable value.
        
        
                 Parameters
                 ----------
                 independent_variable_value : float
                     Value of independent variable at which the interpolation is to be performed.
        
                 Returns
                 -------
                 np.array
                     Interpolated dependent variable value, using implemented algorithm at requested independent variable value
        """

    @property
    def dependent_values(self) -> list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
        """
                 Returns the dependent variable values used by the interpolator.
        
                 Returns the dependent variable values used by the interpolator. This is a read-only property.
        
                 Returns
                 -------
                 list[np.ndarray]
                     Dependent variable values used by the interpolator
        """

    @property
    def independent_values(self) -> list[float]:
        """
                 Returns the independent variable values used by the interpolator.
        
                 Returns the independent variable values used by the interpolator. This is a read-only property.
        
                 Returns
                 -------
                 list[float]
                     Independent variable values used by the interpolator
        """

class OneDimensionalInterpolatorVectorTimeObject:
    """Same as :func:`~OneDimensionalInterpolatorVector`, but using the high-resolution :func:`~Time` type used as independent variable for interpolation; created using :func:`~create_one_dimensional_vector_interpolator_time_object`
    
     """

    def interpolate(self, independent_variable_value: ...) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]:
        """
                This function performs the interpolation at the requested independent variable value.
        
                Parameters
                ----------
                independent_variable_value : float
                    Value of independent variable at which the interpolation is to be performed.
                Returns
                -------
                np.array
                    Interpolated dependent variable value, using implemented algorithm at requested independent variable value
        """

def create_one_dimensional_matrix_interpolator(data_to_interpolate: dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]], interpolator_settings: InterpolatorSettings, data_first_derivatives: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]=[]) -> OneDimensionalInterpolatorMatrix:
    """Function to create an interpolator for matrix dependent variables.
    
    As :func:`create_one_dimensional_scalar_interpolator`, but with matrices (2-dimensional arrays) as dependent variables
    
    
    Parameters
    ----------
    data_to_interpolate : dict[float, np.array]
        Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
    interpolator_settings : InterpolatorSettings
        Settings that define the type of interpolator that is to be used
    data_first_derivatives : list[np.ndarray] = []
        List of first derivative dependent variables w.r.t. independent variable from which the interpolation is to be performed. Must be of the same size as the number of data points in ``data_to_interpolate``. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator).
    Returns
    -------
    OneDimensionalInterpolatorMatrix
        Interpolator object"""

def create_one_dimensional_matrix_interpolator_time_object(data_to_interpolate: dict[..., typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]], interpolator_settings: InterpolatorSettings, data_first_derivatives: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]=[]) -> OneDimensionalInterpolatorMatrixTimeObject:
    """Same as :func:`~create_one_dimensional_matrix_interpolator`, but using the high-resolution :func:`~Time` type used as independent variable for interpolation
    
    Parameters
    ----------
    data_to_interpolate : dict[Time, np.array]
        Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
    interpolator_settings : InterpolatorSettings
        Settings that define the type of interpolator that is to be used
    data_first_derivatives : list[np.ndarray] = []
        List of first derivative dependent variables w.r.t. independent variable from which the interpolation is to be performed. Must be of the same size as the number of data points in ``data_to_interpolate``. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator).
    Returns
    -------
    OneDimensionalInterpolatorMatrixTimeObject
        Interpolator object"""

def create_one_dimensional_scalar_interpolator(data_to_interpolate: dict[float, float], interpolator_settings: InterpolatorSettings, data_first_derivatives: list[float]=[]) -> OneDimensionalInterpolatorScalar:
    """Function to create an interpolator for scalar dependent variables.
    
    Function to create an interpolator for scalar dependent variables, with a single float independent
    variable. This function takes the interpolator settings, and the data that is to be interpolated,
    as input to create the object that can perform the actual interpolation
    
    
    Parameters
    ----------
    data_to_interpolate : dict[float, float]
        Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
    interpolator_settings : InterpolatorSettings
        Settings that define the type of interpolator that is to be used
    data_first_derivatives : list[float] = []
        List of first derivative dependent variables w.r.t. independent variable from which the interpolation is to be performed. Must be of the same size as the number of data points in ``data_to_interpolate``. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator).
    Returns
    -------
    OneDimensionalInterpolatorScalar
        Interpolator object"""

def create_one_dimensional_scalar_interpolator_time_object(data_to_interpolate: dict[float, float], interpolator_settings: InterpolatorSettings, data_first_derivatives: list[float]=[]) -> OneDimensionalInterpolatorScalar:
    """Same as :func:`~create_one_dimensional_scalar_interpolator`, but using the high-resolution :func:`~Time` type used as independent variable for interpolation
    
    Parameters
    ----------
    data_to_interpolate : dict[Time, float]
        Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
    interpolator_settings : InterpolatorSettings
        Settings that define the type of interpolator that is to be used
    data_first_derivatives : list[float] = []
        List of first derivative dependent variables w.r.t. independent variable from which the interpolation is to be performed. Must be of the same size as the number of data points in ``data_to_interpolate``. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator).
    Returns
    -------
    OneDimensionalInterpolatorScalarTimeObject
        Interpolator object"""

def create_one_dimensional_vector_interpolator(data_to_interpolate: dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], interpolator_settings: InterpolatorSettings, data_first_derivatives: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]=[]) -> OneDimensionalInterpolatorVector:
    """Function to create an interpolator for vector dependent variables.
    
    As :func:`create_one_dimensional_scalar_interpolator`, but with vectors as dependent variables
    
    
    Parameters
    ----------
    data_to_interpolate : dict[float, np.array]
        Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
    interpolator_settings : InterpolatorSettings
        Settings that define the type of interpolator that is to be used
    data_first_derivatives : list[np.ndarray] = []
        List of first derivative dependent variables w.r.t. independent variable from which the interpolation is to be performed. Must be of the same size as the number of data points in ``data_to_interpolate``. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator).
    Returns
    -------
    OneDimensionalInterpolatorVector
        Interpolator object"""

def create_one_dimensional_vector_interpolator_time_object(data_to_interpolate: dict[..., typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]], interpolator_settings: InterpolatorSettings, data_first_derivatives: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]=[]) -> OneDimensionalInterpolatorVectorTimeObject:
    """Same as :func:`~create_one_dimensional_vector_interpolator`, but using the high-resolution :func:`~Time` type used as independent variable for interpolation
    
    Parameters
    ----------
    data_to_interpolate : dict[Time, np.array]
        Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
    interpolator_settings : InterpolatorSettings
        Settings that define the type of interpolator that is to be used
    data_first_derivatives : list[np.ndarray] = []
        List of first derivative dependent variables w.r.t. independent variable from which the interpolation is to be performed. Must be of the same size as the number of data points in ``data_to_interpolate``. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator).
    Returns
    -------
    OneDimensionalInterpolatorVectorTimeObject
        Interpolator object"""

def cubic_spline_interpolation(lookup_scheme: AvailableLookupScheme=..., boundary_interpolation: BoundaryInterpolationType=...) -> InterpolatorSettings:
    """Function to create settings for cubic spline interpolation.
    
    Function to create settings for cubic spline interpolation, where the interpolator
    defines a cubic curve polynomial curve between each two subsequent intervals of the
    independent variable input data. The curve has continuous value, first derivative and
    second derivate between subsequent intervals. As boundary condition, the spline has
    a zero second derivative imposed at the upper and lower boundaries of the interpolation
    domain.
    
    
    Parameters
    ----------
    lookup_scheme : AvailableLookupScheme, default = hunting_algorithm
        Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
    boundary_interpolation : BoundaryInterpolationType, default = extrapolate_at_boundary_with_warning
        Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
    Returns
    -------
    InterpolatorSettings
        Cubic spline settings object"""

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
    
    
    Parameters
    ----------
    lookup_scheme : AvailableLookupScheme, default = hunting_algorithm
        Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
    boundary_interpolation : BoundaryInterpolationType, default = extrapolate_at_boundary_with_warning
        Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
    Returns
    -------
    InterpolatorSettings
        Hermite spline interpolation settings object"""

def interpolator_generation_settings(interpolator_settings: InterpolatorSettings, initial_time: float, final_time: float, time_step: float) -> InterpolatorGenerationSettings:
    """No documentation found."""

def interpolator_generation_settings_time_object(interpolator_settings: InterpolatorSettings, initial_time: ..., final_time: ..., time_step: ...) -> InterpolatorGenerationSettingsTimeObject:
    """No documentation found."""

def lagrange_interpolation(number_of_points: int, lookup_scheme: AvailableLookupScheme=..., boundary_interpolation: BoundaryInterpolationType=..., lagrange_boundary_handling: LagrangeInterpolatorBoundaryHandling=...) -> InterpolatorSettings:
    """Function to create settings for cubic Lagrange interpolation.
    
    Function to create settings for piecewise cubic Lagrange interpolation.  This is typically the interpolator of highest accuracy that is available.
    The Lagrange interpolator uses :math:`m` consecutive points (input to this function) from the input independent variables :math:`[t_{0}...t_{N}]`
    to create the polynomial of order :math:`m-1` that interpolates these points. From here on, we assume :math:`m` is even.
    The algorithm that is used (see `here <https://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html>`_ for mathematical details
    on interpolating Lagrange polynomials) works as follows:
    
    * The nearest lower neighbor of the data point :math:`t` at which the state :math:`\\mathbf{x}` is to be interpolated is determined, and denoted :math:`t_{i}`.
    * An interpolating Lagrange polynomial is constructed from the consecutive data points :math:`[t_{i-(m/2-1)}...t_{i+m}]`
    * This resulting interpolating polynomial is *only* used in the interval :math:`[t_{i}...t_{i+1}]`, to prevent `Runge's phenomenon <https://en.wikipedia.org/wiki/Runge%27s_phenomenon>`_.
    
    For instance, if :math:`m=8` we use a :math:`7^{th}` order polynomial that interpolates a contiguous set of
    8 data points out of the full data set. Normally, the interpolating polynomial is only used between the
    :math:`4^{th}` and :math:`5^{th}` data point, where it will typically be of good accuracy. Consequently,
    a separate interpolating polynomial (using data over a span of :math:`m` consecutive points) is used for
    each single interval :math:`[t_{i}...t_{i+1}]` (with the exception of the boundaries, see below).
    
    .. warning:: Issues can occur if the data point :math:`t` at which the interpolation is to be
                 performed is close to :math:`t_{0}` or :math:`t_{N}`. In those case, there is not sufficient
                 data to construct the interpolating polynomial *and* to only use this interpolating polynomial
                 between the middle two data points that were used to it. In these cases, the user has a number of
                 options (all defined by an entry of the :class:`LagrangeInterpolatorBoundaryHandling` variable,
                 used as input to this function). In short, interpolation between the first and last :math:`m/2`
                 data points will lead to degraded results, warnings, or termination.
    
    
    Parameters
    ----------
    number_of_points : int
        Number of consecutive data points that are used to construct a single interpolating polynomial.
    lookup_scheme : AvailableLookupScheme, default = hunting_algorithm
        Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
    boundary_interpolation : BoundaryInterpolationType, default = extrapolate_at_boundary_with_warning
        Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
    lagrange_boundary_handling : LagrangeInterpolatorBoundaryHandling, default = lagrange_cubic_spline_boundary_interpolation
        Interpolator behaviour that is to be used at the boundaries of the domain, where the regular algorithm cannot be executed.
    Returns
    -------
    LagrangeInterpolatorSettings
        Lagrange interpolation settings object"""

def linear_interpolation(lookup_scheme: AvailableLookupScheme=..., boundary_interpolation: BoundaryInterpolationType=...) -> InterpolatorSettings:
    """Function to create settings for linear interpolation.
    
    Function to create settings for linear interpolation, where the interpolator
    defines a linear curve between each two subsequent intervals of the
    independent variable input data.
    
    
    Parameters
    ----------
    lookup_scheme : AvailableLookupScheme, default = hunting_algorithm
        Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
    boundary_interpolation : BoundaryInterpolationType, default=extrapolate_at_boundary_with_warning
        Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
    Returns
    -------
    InterpolatorSettings
        Linear interpolation settings object"""

def piecewise_constant_interpolation(lookup_scheme: AvailableLookupScheme=..., boundary_interpolation: BoundaryInterpolationType=...) -> InterpolatorSettings:
    """Function to create settings for piecewise constant interpolation.
    
    Function to create settings for piecewise constant interpolation. If interpolator
    is to return the value at :math:`t`, and :math:`t_{i}\\le t < t_{i+1}`, the interpolator
    returns :math:`\\mathbf{x}_{i}`
    
    
    Parameters
    ----------
    lookup_scheme : AvailableLookupScheme, default = hunting_algorithm
        Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
    boundary_interpolation : BoundaryInterpolationType, default = extrapolate_at_boundary_with_warning
        Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
    Returns
    -------
    InterpolatorSettings
        Piecewise constant interpolation settings object"""
binary_search: AvailableLookupScheme
extrapolate_at_boundary: BoundaryInterpolationType
extrapolate_at_boundary_with_warning: BoundaryInterpolationType
hunting_algorithm: AvailableLookupScheme
lagrange_cubic_spline_boundary_interpolation: LagrangeInterpolatorBoundaryHandling
lagrange_no_boundary_interpolation: LagrangeInterpolatorBoundaryHandling
throw_exception_at_boundary: BoundaryInterpolationType
use_boundary_value: BoundaryInterpolationType
use_boundary_value_with_warning: BoundaryInterpolationType