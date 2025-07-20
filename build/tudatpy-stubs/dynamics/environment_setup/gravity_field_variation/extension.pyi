import numpy
import pybind11_stubgen.typing_ext
from ....math import interpolators
import typing
__all__ = ['BasicSolidBodyGravityFieldVariationSettings', 'BodyDeformationTypes', 'GravityFieldVariationSettings', 'basic_solid_body', 'iers_2010_tidal', 'mode_coupled_solid_body_tide', 'ocean_tide', 'periodic', 'periodic_variation', 'pole_tide', 'polynomial', 'polynomial_variation', 'single_period_periodic', 'single_power_polynomial', 'solid_body_tide', 'solid_body_tide_complex_k', 'solid_body_tide_degree_order_variable_complex_k', 'solid_body_tide_degree_order_variable_k', 'solid_body_tide_degree_variable_complex_k', 'solid_body_tide_degree_variable_k', 'solid_multi_body_tide_degree_order_variable_k', 'tabulated', 'tabulated_deformation']

class BasicSolidBodyGravityFieldVariationSettings(GravityFieldVariationSettings):
    """Class for providing settings for solid body tidal gravity field variations, derived from GravityFieldVariationSettings."""

class BodyDeformationTypes:
    """Enumeration listing the different types of gravity and/or shape field variation models available in tudat.
    Note that for some types, only one of the two types of variations is available
    
          
    
    Members:
    
      basic_solid_body : 
    
    Basic tidal variation model, assuming a single constant Love number for the variation
    
    
    
      iers_2010_tidal : 
    
    High-fidelity Earth tidal variation model based on IERS 2010 conventions
    
    
    
      tabulated_deformation : 
    
    Variation model using interpolated tabular data for the variation model
    
    
    
      periodic_variation : 
    
    Variation model using purely sinusoidal variations (as a function of time) for gravity field coefficients
    
    
    
      polynomial_variation : 
    
    Variation model using polynomial functions of time for gravity field coefficients variations
    
    
    
      ocean_tide : 
    
    Variation model due to ocean tides
    
    
    
      pole_tide : 
    
    Variation model due to pole tides
    
          """
    __members__: typing.ClassVar[dict[str, BodyDeformationTypes]]
    basic_solid_body: typing.ClassVar[BodyDeformationTypes]
    iers_2010_tidal: typing.ClassVar[BodyDeformationTypes]
    ocean_tide: typing.ClassVar[BodyDeformationTypes]
    periodic_variation: typing.ClassVar[BodyDeformationTypes]
    pole_tide: typing.ClassVar[BodyDeformationTypes]
    polynomial_variation: typing.ClassVar[BodyDeformationTypes]
    tabulated_deformation: typing.ClassVar[BodyDeformationTypes]

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
    """Base class for providing settings for gravity field variations."""

def mode_coupled_solid_body_tide(deforming_bodies: list[str], love_numbers: dict[tuple[int, int], dict[tuple[int, int], float]]) -> GravityFieldVariationSettings:
    """Function for creating solid body tides with coupling between different forcing degree and orders.
    
    Function for creating solid body tides with coupling between different forcing degree and orders, similar to
    :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k`, but
    with multiple tide-raising bodies and the option to have forcing and response at different degrees and orders
    
    .. math::
        \\Delta \\bar{C}_{l'm'}-i\\Delta \\bar{S}_{l'm'}=\\sum_{j}\\frac{1}{2l+1}k_{lm}^{l'm'}\\frac{\\mu_{j}}{\\mu}\\left(\\frac{R}{r_{j}}\\right)^{l+1}\\bar{P}_{lm}(\\sin\\phi_{j})\\left(\\cos m\\theta_{j}-i\\sin m\\theta_{j}\\right)
    
    where quantities without subscripts represent properties of the body :math:`B` for which the gravity field's time variations are computed.
    Here, :math:`\\mu_{j}` is the gravitational parameter of the tide-raising body :math:`j`, and :math:`r_{j}`, :math:`\\phi_{j}` and :math:`\\theta_{j}` represent the spherical position
    of body :math:`j` in a frame fixed to body :math:`B`, and :math:`\\bar{P}_{lm}` are the fully normalized associated Legendre polynomials.
    
    The quantities :math:`l` and :math:`m` denote the degree and order of the forcing, while :math:`l'` and :math:`m'` denote the degree and order of the
    response
    
    
    Parameters
    ----------
    tide_raising_body : str
        Name of body raising the tide.
    love_number_per_degree : dict[tuple[int, int], dict[tuple[int,int],float]]
        Double dictionary of Love numbers for each combination for forcing and response degree and order, and the value (``float``) the associated
        mode-coupled Love number :math:`k_{lm}^{l'm'}`. The first tuple is the forcing degree and order :math:`l,m`, the second tuple is the response degree and order :math:`l',m'`
    
    Returns
    -------
    GravityFieldVariationSettings
        Instance of a :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.GravityFieldVariationSettings` derived class containing required settings"""

def periodic(cosine_coefficient_amplitudes_cosine_time: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]], cosine_coefficient_amplitudes_sine_time: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]], sine_coefficient_amplitudes_cosine_time: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]], sine_coefficient_amplitudes_sine_time: list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]], angular_frequencies: list[float], reference_epoch: float, minimum_degree: int=2, minimum_order: int=0) -> GravityFieldVariationSettings:
    """Function for creating time-periodic gravity field variations.
    
    Function for creating gravity field variations that are a superposition of purely sinusoidal variations in the gravity field coefficients.
    The cosine and sine coefficient variations at degree and order :math:`l` and :math:`m` are computed from:
    
    .. math::
         \\Delta{\\bar{C}}_{lm}=\\sum_{i=1}^{N}\\left(A_{i,\\bar{C}_{lm}}\\cos\\left(f_{i}(t-t_{0})\\right) + B_{i,\\bar{C}_{lm}}\\sin\\left(f_{i}(t-t_{0})\\right) \\right)\\\\
         \\Delta{\\bar{S}}_{lm}=\\sum_{i=1}^{N}\\left(A_{i,\\bar{S}_{lm}}\\cos\\left(f_{i}(t-t_{0})\\right) + B_{i,\\bar{S}_{lm}}\\sin\\left(f_{i}(t-t_{0})\\right) \\right)
    
    The summation is over all :math:`N` frequencies that are provided by the user. For each frequency, the user provides a block of coefficients :math:`A_{\\bar{C}_{lm}}`,
    :math:`A_{i,\\bar{S}_{lm}}`, :math:`B_{i,\\bar{S}_{lm}}` and :math:`B_{i,\\bar{S}_{lm}}`. These blocks give the spherical harmonic coefficient variations
    at degree :math:`l_{\\text{min}}` until degree :math:`l_{\\text{min}}+l_{\\text{size}}`, and order :math:`m_{\\text{min}}` until degree :math:`m_{\\text{min}}+m_{\\text{size}}`.
    The :math:`l_{\\text{min}}`  and :math:`m_{\\text{min}}` are defined by the ``minimum_degree`` and ``minimum_order`` inputs. The :math`l_{\\text{size}}`  and :math`m_{\\text{size}}`
    are defined by the size of the matrices in the first four input lists to this function (starting with ``cosine_coefficient_amplitudes_cosine_time``).
    
    For instance, if ``minimum_degree`` is set to 2, and ``minimum_order`` to 0, and each matrix block for coefficient amplitudes provided is a 1x3 matrix, these
    amplitudes make of the variations in :math:`l,m=2,0`, :math:`l,m=2,1` and :math:`l,m=2,2` respectively.
    
    Parameters
    ----------
    cosine_coefficient_amplitudes_cosine_time : list[np.array]
        List of coefficient amplitude blocks :math:`A_{i,\\bar{C}_{lm}}`, with each entry in the list corresponding to the frequency provided by the same entry in ``angular_frequencies``.
        The first entry in each matrix block provides the coefficient variation amplitude at degree equal to ``minimum_degree`` and order equal to ``minimum_order``.
    cosine_coefficient_amplitudes_sine_time : list[np.array]
        List of coefficient amplitude blocks :math:`B_{i,\\bar{C}_{lm}}`, with each entry in the list corresponding to the frequency provided by the same entry in ``angular_frequencies``.
        The first entry in each matrix block provides the coefficient variation amplitude at degree equal to ``minimum_degree`` and order equal to ``minimum_order``.
    sine_coefficient_amplitudes_cosine_time : list[np.array]
        List of coefficient amplitude blocks :math:`A_{i,\\bar{S}_{lm}}`, with each entry in the list corresponding to the frequency provided by the same entry in ``angular_frequencies``.
        The first entry in each matrix block provides the coefficient variation amplitude at degree equal to ``minimum_degree`` and order equal to ``minimum_order``.
        Note that if ``minimum_order`` is equal to 0, the first column of values in each matrix will be unused (since there is no :math:`S_{l0}` coefficient).
    sine_coefficient_amplitudes_sine_time : list[np.array]
        List of coefficient amplitude blocks :math:`B_{i,\\bar{S}_{lm}}`, with each entry in the list corresponding to the frequency provided by the same entry in ``angular_frequencies``.
        The first entry in each matrix block provides the coefficient variation amplitude at degree equal to ``minimum_degree`` and order equal to ``minimum_order``.
        Note that if ``minimum_order`` is equal to 0, the first column of values in each matrix will be unused (since there is no :math:`S_{l0}` coefficient).
    angular_frequencies : list[float]
        List of angular frequencies (in rad/s) at which variations are to be added
    reference_epoch: float
        Reference epoch :math:`t_{0}` for the variations
    minimum_degree: int
        Minimum degree :math:`l_{\\text{min}}` of gravity field variations
    minimum_order: int
        Minimum order :math:`m_{\\text{min}}` of gravity field variations
    
    Returns
    -------
    GravityFieldVariationSettings
        Instance of a :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.GravityFieldVariationSettings` derived class containing required settings"""

def polynomial(cosine_amplitudes_per_power: dict[int, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]], sine_amplitudes_per_power: dict[int, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]], reference_epoch: float, minimum_degree: int=2, minimum_order: int=0) -> GravityFieldVariationSettings:
    """Function for creating time-polynomial gravity field variations.
    
    Function for creating gravity field variations that are a superposition of purely polynomial variations in the gravity field coefficients.
    The cosine and sine coefficient variations at degree and order :math:`l` and :math:`m` are computed from:
    
    .. math::
         \\Delta{\\bar{C}}_{lm}=\\sum_{i=1}^{N}\\left(K_{i,\\bar{C}_{lm}} (t-t_{0})^{p_{i}} \\right)\\\\
         \\Delta{\\bar{S}}_{lm}=\\sum_{i=1}^{N}\\left(K_{i,\\bar{S}_{lm}} (t-t_{0})^{p_{i}} \\right)
    
    The summation is over all :math:`N` polynomial exponents that are provided by the user. For each exponent, the user provides a block of coefficients :math:`K_{i,\\bar{C}_{lm}}`,
    and :math:`K_{i,\\bar{S}_{lm}}`. These blocks give the spherical harmonic coefficient variations
    at degree :math:`l_{\\text{min}}` until degree :math:`l_{\\text{min}}+l_{\\text{size}}`, and order :math:`m_{\\text{min}}` until degree :math:`m_{\\text{min}}+m_{\\text{size}}`.
    The :math:`l_{\\text{min}}`  and :math:`m_{\\text{min}}` are defined by the ``minimum_degree`` and ``minimum_order`` inputs. The :math`l_{\\text{size}}`  and :math`m_{\\text{size}}`
    are defined by the size of the matrices in the first two input lists to this function (starting with ``cosine_amplitudes_per_power``).
    
    For instance, if ``minimum_degree`` is set to 2, and ``minimum_order`` to 0, and each matrix block for coefficient amplitudes provided is a 1x3 matrix, these
    amplitudes make of the variations in :math:`l,m=2,0`, :math:`l,m=2,1` and :math:`l,m=2,2` respectively.
    
    Parameters
    ----------
    cosine_amplitudes_per_power : dict[int, np.array]
        Dictionary of cosine coefficient amplitude blocks, with each key in the list corresponding to the polynomial power :math:`p_{i}`, and
        the dictionary value the corresponding :math:`K_{i,\\bar{C}_{lm}}`. The first entry in each matrix block provides the coefficient variation amplitude at degree equal to ``minimum_degree`` and order equal to ``minimum_order``.
    sine_amplitudes_per_power : list[np.array]
        Dictionary of sine coefficient amplitude blocks, with each key in the list corresponding to the polynomial power :math:`p_{i}`, and
        the dictionary value the corresponding :math:`K_{i,\\bar{S}_{lm}}`. The first entry in each matrix block provides the coefficient variation amplitude at degree equal to ``minimum_degree`` and order equal to ``minimum_order``.
        Note that if ``minimum_order`` is equal to 0, the first column of values in each matrix will be unused (since there is no :math:`S_{l0}` coefficient).
    reference_epoch: float
        Reference epoch :math:`t_{0}` for the variations
    minimum_degree: int
        Minimum degree :math:`l_{\\text{min}}` of gravity field variations
    minimum_order: int
        Minimum order :math:`m_{\\text{min}}` of gravity field variations
    
    Returns
    -------
    GravityFieldVariationSettings
        Instance of a :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.GravityFieldVariationSettings` derived class containing required settings"""

def single_period_periodic(cosine_coefficient_amplitude_cosine_time: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], cosine_coefficient_amplitude_sine_time: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], sine_coefficient_amplitude_cosine_time: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], sine_coefficient_amplitude_sine_time: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], angular_frequency: float, reference_epoch: float, minimum_degree: int=2, minimum_order: int=0) -> GravityFieldVariationSettings:
    """Function for creating time-periodic gravity field variations at a single variation period.
    
    Function for creating time-periodic gravity field variations at a single variation period, same as :math:`~periodic`, but with a single
    frequency :math:`f`.
    
    Parameters
    ----------
    cosine_coefficient_amplitude_cosine_time : np.array
        Coefficient amplitude block :math:`A_{i,\\bar{C}_{lm}}`. The first entry in each matrix block provides the coefficient variation amplitude at degree equal to ``minimum_degree`` and order equal to ``minimum_order``.
    cosine_coefficient_amplitude_sine_time : np.array
        Coefficient amplitude block :math:`B_{i,\\bar{C}_{lm}}`. The first entry in each matrix block provides the coefficient variation amplitude at degree equal to ``minimum_degree`` and order equal to ``minimum_order``.
    sine_coefficient_amplitude_cosine_time : np.array
        Coefficient amplitude block :math:`A_{i,\\bar{S}_{lm}}`. The first entry in each matrix block provides the coefficient variation amplitude at degree equal to ``minimum_degree`` and order equal to ``minimum_order``.
        Note that if ``minimum_order`` is equal to 0, the first column of values in each matrix will be unused (since there is no :math:`S_{l0}` coefficient).
    sine_coefficient_amplitude_sine_time : np.array
        Coefficient amplitude block :math:`B_{i,\\bar{S}_{lm}}`. The first entry in each matrix block provides the coefficient variation amplitude at degree equal to ``minimum_degree`` and order equal to ``minimum_order``.
        Note that if ``minimum_order`` is equal to 0, the first column of values in each matrix will be unused (since there is no :math:`S_{l0}` coefficient).
    angular_frequency : list[float]
        Angular frequencies (in rad/s) at which variations are to be added
    reference_epoch: float
        Reference epoch :math:`t_{0}` for the variations
    minimum_degree: int
        Minimum degree :math:`l_{\\text{min}}` of gravity field variations
    minimum_order: int
        Minimum order :math:`m_{\\text{min}}` of gravity field variations
    
    Returns
    -------
    GravityFieldVariationSettings
        Instance of a :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.GravityFieldVariationSettings` derived class containing required settings"""

def single_power_polynomial(cosine_amplitudes: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], sine_amplitudes: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], polynomial_power: int, reference_epoch: float, minimum_degree: int=2, minimum_order: int=0) -> GravityFieldVariationSettings:
    """Function for creating polynomial gravity field variations at a single variation period.
    
    Function for creating polynomial gravity field variations at a single variation period, same as :math:`~polynomial`, but with a single
    polynmial power :math:`p`.
    
    Parameters
    ----------
    cosine_amplitudes_per_power : np.array
        Cosine coefficient amplitude block :math:`K_{\\bar{C}_{lm}}`. The first entry in each matrix block provides the coefficient variation amplitude at degree equal to ``minimum_degree`` and order equal to ``minimum_order``.
    sine_amplitudes_per_power : list[np.array]
        Sine coefficient amplitude block :math:`K_{\\bar{S}_{lm}}`. The first entry in each matrix block provides the coefficient variation amplitude at degree equal to ``minimum_degree`` and order equal to ``minimum_order``.
        Note that if ``minimum_order`` is equal to 0, the first column of values in each matrix will be unused (since there is no :math:`S_{l0}` coefficient).
    polynomial_power: float
        Polynomial power :math:`p` used in the computation
    reference_epoch: float
        Reference epoch :math:`t_{0}` for the variations
    minimum_degree: int
        Minimum degree :math:`l_{\\text{min}}` of gravity field variations
    minimum_order: int
        Minimum order :math:`m_{\\text{min}}` of gravity field variations
    
    Returns
    -------
    GravityFieldVariationSettings
        Instance of a :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.GravityFieldVariationSettings` derived class containing required settings"""

def solid_body_tide(tide_raising_body: str, love_number: float, degree: int) -> GravityFieldVariationSettings:
    """Function for creating solid body tides.
    
    Function for creating solid body tides, using a single real Love number :math:`k_{l}` at a single degree :math:`l` (e.g. :math:`k_{2}`, :math:`k_{3}`, etc.).
    This function evaluates Eq. (6.6) from the IERS Conventions 2010, with real :math:`k_{l}=k_{lm}`, a single value of :math:`l` and
    a single tide-raising body :math:`j`.
    
     .. math::
         \\Delta \\bar{C}_{lm}-i\\Delta \\bar{S}_{lm}=\\frac{1}{2l+1}k_{l}\\frac{\\mu_{j}}{\\mu}\\left(\\frac{R}{r_{j}}\\right)^{l+1}\\bar{P}_{lm}(\\sin\\phi_{j})\\left(\\cos m\\theta_{j}-i\\sin m\\theta_{j}\\right)
    
    where quantities without subscripts represent properties of the body :math:`B` for which the gravity field's time variations are computed.
    Here, :math:`\\mu_{j}` is the gravitational parameter of the tide-raising body :math:`j`, and :math:`r_{j}`, :math:`\\phi_{j}` and :math:`\\theta_{j}` represent the spherical position
    of body :math:`j` in a frame fixed to body :math:`B`, and :math:`\\bar{P}_{lm}` are the fully normalized associated Legendre polynomials
    
    Parameters
    ----------
    tide_raising_body : str
        Name of body raising the tide.
    love_number : float
        Constant real Love number to use for body undergoing deformation, at the spherical harmonic degree defined by 'degree' input.
    degree : int
        Degree of the spherical harmonic gravity field, and associated Love number, that is to be considered.
    
    Returns
    -------
    BasicSolidBodyGravityFieldVariationSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
    
    Examples
    --------
    In this example, we create gravity field variations of Earth for a tide raised by the Moon, with a single Love number :math:`k_{2}` of 0.301, and add it to the list of gravity field variations
    
    .. code-block:: python
    
        tide_raising_body = "Moon"
        degree = 2
        love_number = 0.301
        gravity_field_variation_list = list()
        gravity_field_variation_list.append( environment_setup.gravity_field_variation.solid_body_tide(
            tide_raising_body, love_number, degree )
        body_settings.get( "Earth" ).gravity_field_variation_settings = gravity_field_variation_list"""

def solid_body_tide_complex_k(tide_raising_body: str, love_number: complex, degree: int) -> GravityFieldVariationSettings:
    """Function for creating solid body tides.
    
    As :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide`, but with complex value for the Love number.
    
    
    Parameters
    ----------
    tide_raising_body : str
        Name of body raising the tide.
    love_number : complex
        Constant real Love number to use for body undergoing deformation, at the spherical harmonic degree defined by 'degree' input.
    degree : int
        Degree of the spherical harmonic gravity field, and associated Love number, that is to be considered.
    
    Returns
    -------
    BasicSolidBodyGravityFieldVariationSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class"""

def solid_body_tide_degree_order_variable_complex_k(tide_raising_body: str, love_number_per_degree_and_order: dict[int, list[complex]]) -> GravityFieldVariationSettings:
    """Function for creating solid body tides.
    
    As :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k`,
    but with complex values for the Love number.
    
    Parameters
    ----------
    tide_raising_body : str
        Name of body raising the tide.
    love_number_per_degree : dict( int, list( complex ) )
        Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the list of Love numbers :math:`k_{lm}` at this degree. Note that, for Love numbers at degree :math:`l`, the associated list can contain up to :math:`l+1` entries, representing the Love numbers (in order) :math:`k_{l0}`, :math:`k_{l1}`...:math:`k_{ll}`.
    
    Returns
    -------
    BasicSolidBodyGravityFieldVariationSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class"""

def solid_body_tide_degree_order_variable_k(tide_raising_body: str, love_number_per_degree_and_order: dict[int, list[float]]) -> GravityFieldVariationSettings:
    """Function for creating solid body tides.
    
    Function for creating solid body tides, using a set of real, separate, Love numbers at any number of degrees and orders
    (e.g. :math:`k_{20}`, :math:`k_{21}`, :math:`k_{22}`, :math:`k_{30}`, etc.).
    This function evaluates Eq. (6.6) from the IERS Conventions 2010, with a set of real values :math:`k_{lm}`,
    at a set of values of :math:`l` and a single tide-raising body :math:`j`.
    
    The mathematical model of this function is effectively nearly equal to that of :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide`
    with the modification that :math:`k_{l}\\rightarrow k_{lm}`.
    
    Parameters
    ----------
    tide_raising_body : str
        Name of body raising the tide.
    love_number_per_degree_and_order : dict( int, list( float ) )
        Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the list of Love numbers :math:`k_{lm}` at this degree. Note that, for Love numbers at degree :math:`l`, the associated list can contain up to :math:`l+1` entries, representing the Love numbers (in order) :math:`k_{l0}`, :math:`k_{l1}`... :math:`k_{ll}`.
    
    Returns
    -------
    BasicSolidBodyGravityFieldVariationSettings
         Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
    
    
    Examples
    --------
    In this example, we create gravity field variations of the Moon, for a tide raised by Earth, with a Love numbers :math:`k_{20}=0.024615`, :math:`k_{21}=0.023915` and :math:`k_{21}=0.024852`, and add it to the list of gravity field variations
    
    .. code-block:: python
    
        tide_raising_body = "Earth"
        love_numbers = dict( )
        love_numbers[ 2 ] = list( )
        love_numbers[ 2 ].append( 0.024615 )
        love_numbers[ 2 ].append( 0     .023915 )
        love_numbers[ 2 ].append( 0.024852 )
        gravity_field_variation_list = list()
        gravity_field_variation_list.append( environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k(
            tide_raising_body, love_numbers )
        body_settings.get( "Earth" ).gravity_field_variation_settings = gravity_field_variation_list"""

def solid_body_tide_degree_variable_complex_k(tide_raising_body: str, love_number_per_degree: dict[int, complex]) -> GravityFieldVariationSettings:
    """Function for creating solid body tides.
    
    As :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide_degree_variable_k`,
    but with complex values for the Love numbers.
    
    
    Parameters
    ----------
    tide_raising_body : str
        Name of body raising the tide.
    love_number_per_degree : dict( int, complex )
        Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the Love number :math:`k_{l}` itself.
    
    Returns
    -------
     BasicSolidBodyGravityFieldVariationSettings
         Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class"""

def solid_body_tide_degree_variable_k(tide_raising_body: str, love_number_per_degree: dict[int, float]) -> GravityFieldVariationSettings:
    """Function for creating solid body tides.
    
    Function for creating solid body tides, using a set of real, separate, Love numbers at any number of degrees
    (e.g. :math:`k_{2}`, :math:`k_{3}`, etc.). This output and mathematical model of this function is effectively identical to a list of outputs to
    :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide`, with differing degrees and associated Love numbers.
    This function evaluates Eq. (6.6) from the IERS Conventions 2010, with a set of real values :math:`k_{l}=k_{lm}`, at a set of values of :math:`l`
    and a single tide-raising body :math:`j`.
    
    Parameters
    ----------
    tide_raising_body : str
        Name of body raising the tide.
    love_number_per_degree : dict( int, float )
        Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the Love number :math:`k_{l}` itself
    
    Returns
    -------
    BasicSolidBodyGravityFieldVariationSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
    
    Examples
    --------
    In this example, we create gravity field variations of Earth for a tide raised by the Moon, with a Love numbers :math:`k_{2}=0.301`, and :math:`k_{3}=0.09`, and add it to the list of gravity field variations
    
    .. code-block:: python
    
        tide_raising_body = "Moon"
        love_numbers = dict( )
        love_numbers[ 2 ] = 0.301
        love_numbers[ 3 ] = 0.09
        gravity_field_variation_list = list()
        gravity_field_variation_list.append( environment_setup.gravity_field_variation.solid_body_tide_degree_variable_k(
            tide_raising_body, love_numbers )
        body_settings.get( "Earth" ).gravity_field_variation_settings = gravity_field_variation_list"""

def solid_multi_body_tide_degree_order_variable_k(tide_raising_bodies: list[str], love_number_per_degree_and_order: dict[int, list[float]]) -> GravityFieldVariationSettings:
    """Function for creating solid body tides raised by multiple bodies.
    
    As :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k`, but with the same
    set of Love numbers used for tides raised by multiple bodies :math:`j`
    
    Parameters
    ----------
    tide_raising_bodies : list[str]
        Names of bodies raising the tide.
    love_number_per_degree_and_order : dict( int, list( float ) )
        Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the list of Love numbers :math:`k_{lm}` at this degree. Note that, for Love numbers at degree :math:`l`, the associated list can contain up to :math:`l+1` entries, representing the Love numbers (in order) :math:`k_{l0}`, :math:`k_{l1}`... :math:`k_{ll}`.
    
    Returns
    -------
    BasicSolidBodyGravityFieldVariationSettings
         Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class"""

def tabulated(cosine_variations_table: dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]], sine_variations_table: dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]], minimum_degree: int, minimum_order: int, interpolation_settings: interpolators.InterpolatorSettings) -> GravityFieldVariationSettings:
    """Function for creating gravity field variations from tabular variation values as a function of time.
    
    Function for creating gravity field variations from tabular variation values as a function of time, which are interpolated to compute the instantaneous
    gravity field variation
    
    The user provides a tables of blocks of coefficient variation at epochs :math:`t_{i}` :math:`\\Delta \\bar{C}_{lm}(t_{i})`. These blocks give the spherical harmonic coefficient variations
    at degree :math:`l_{\\text{min}}` until degree :math:`l_{\\text{min}}+l_{\\text{size}}`, and order :math:`m_{\\text{min}}` until degree :math:`m_{\\text{min}}+m_{\\text{size}}`.
    The :math:`l_{\\text{min}}`  and :math:`m_{\\text{min}}` are defined by the ``minimum_degree`` and ``minimum_order`` inputs. The :math`l_{\\text{size}}`  and :math`m_{\\text{size}}`
    are defined by the size of the matrices in the tabulated inputs
    
    For instance, if ``minimum_degree`` is set to 2, and ``minimum_order`` to 0, and each matrix block for coefficient amplitudes provided is a 1x3 matrix, these
    amplitudes make of the variations in :math:`l,m=2,0`, :math:`l,m=2,1` and :math:`l,m=2,2` respectively.
    
    Parameters
    ----------
    cosine_variations_table : dict[float, np.array]
        Dictionary of cosine coefficient variations, with each key in list corresponding to epoch :math:`t_{i}` and the value to the corresponding variation
        :math:`\\Delta \\bar{C}_{lm}(t_{i})`. The first entry in each matrix block provides the coefficient variation at degree equal to ``minimum_degree`` and order equal to
        ``minimum_order``.
    sine_amplitudes_per_power : dict[float, np.array]
        Dictionary of sine coefficient variations, with each key in list corresponding to epoch :math:`t_{i}` and the value to the corresponding variation
        :math:`\\Delta \\bar{C}_{lm}(t_{i})`. The first entry in each matrix block provides the coefficient variation at degree equal to ``minimum_degree`` and order equal to
        ``minimum_order``.
        Note that if ``minimum_order`` is equal to 0, the first column of values in each matrix will be unused (since there is no :math:`S_{l0}` coefficient).
    minimum_degree: int
        Minimum degree :math:`l_{\\text{min}}` of gravity field variations
    minimum_order: int
        Minimum order :math:`m_{\\text{min}}` of gravity field variations
    interpolation_settings: InterpolatorSettings
         Settings that define the type of interpolator that is to be used
    
    Returns
    -------
    GravityFieldVariationSettings
        Instance of a :class:`~tudatpy.dynamics.environment_setup.gravity_field_variation.GravityFieldVariationSettings` derived class containing required settings"""
basic_solid_body: BodyDeformationTypes
iers_2010_tidal: BodyDeformationTypes
ocean_tide: BodyDeformationTypes
periodic_variation: BodyDeformationTypes
pole_tide: BodyDeformationTypes
polynomial_variation: BodyDeformationTypes
tabulated_deformation: BodyDeformationTypes