import numpy
import pybind11_stubgen.typing_ext
import typing
__all__ = ['legendre_normalization_factor', 'normalize_spherical_harmonic_coefficients', 'spherical_harmonic_coefficients_from_inertia', 'unnormalize_spherical_harmonic_coefficients']

def legendre_normalization_factor(degree: int, order: int) -> float:
    """Function to calculate the normalization factor for spherical harmonics at a given degree and order
    
    Function to calculate the normalization factor for spherical harmonics at a given degree and order.
    Specifically, this function returns the value :math:`\\mathcal{N}_{lm}`, computed from:
    
    .. math::
        \\mathcal{N}_{lm}=\\sqrt{\\frac{(2-\\delta_{0m})(2l+1)(l-m)!)}{(l+m)!}}
    
    with :math:`\\delta_{0m}` is the Kronecker Delta function. The following can be used such that the conversion between unnormalized and fully
    normalized spherical harmonic coefficients and Legendre polynomials can be computed from:
    
    .. math::
        [C,S]_{lm}&=\\mathcal{N}_{lm}[\\bar{C},\\bar{S}]_{lm}\\\\
        \\bar{P}_{lm}&=\\mathcal{N}_{lm}{P}_{lm}
    
    with :math:`[C,S]_{lm}` denoting the unnormalized cosine or sine spherical harmonic coefficients at degree :math:`l` and order :math:`m`,
    and :math:`P_{lm}` and :math:`\\bar{P}_{lm}` representing the unnormalized and normalized associated Legendre polynomials at degree :math:`l` and order :math:`m`.
    
    
    Parameters
    ----------
    degree : int
        Spherical harmonic degree :math:`l`
    order : int
        Spherical harmonic order :math:`m`
    Returns
    -------
    float
        Normalization factor :math:`\\mathcal{N}_{lm}`"""

def normalize_spherical_harmonic_coefficients(unnormalized_cosine_coefficients: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], unnormalized_sine_coefficients: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> tuple[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]:
    """Function to normalize spherical harmonic coefficients
    
    Function to normalize spherical harmonic coefficients, using the equations provided in the :func:`~tudatpy.gravitation.astro.legendre_normalization_factor` function.
    
    Parameters
    ----------
    unnormalized_cosine_coefficients : numpy.ndarray
        Matrix for which entry :math:`(i,j)` contains the unnormalized cosine coefficient :math:`C_{lm}`
    unnormalized_sine_coefficients : numpy.ndarray
        Matrix for which entry :math:`(i,j)` contains the unnormalized sine coefficient :math:`S_{lm}`
    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Tuple of two matrices, containing the normalized coefficients :math:`\\bar{C}_{lm}` (first) and :math:`\\bar{S}_{lm}` (second)"""

def spherical_harmonic_coefficients_from_inertia(inertia_tensor: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)], gravitational_parameter: float, reference_radius: float, output_normalized_coefficients: bool=True) -> tuple[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], float]:
    """Function to compute degree-two spherical harmonic coefficients from an inertia tensor
    
    Function to compute degree-two spherical harmonic coefficients :math:`C_{20}`, :math:`C_{21}`, :math:`C_{22}`, :math:`S_{21}`, :math:`S_{22}` and from an inertia tensor :math:`\\mathbf{I}`, according to the following relation"
    
    .. math::
        \\mathbf{I}=M R^2\\left(\\left(\\begin{array}{ccc} \\frac{C_{20}}{3}-2 C_{22} & -2 S_{22} & -C_{21} \\\\ -2 S_{22} & \\frac{C_{20}}{3}+2 C_{22} & -S_{21} \\\\ -C_{21} & -S_{21} & -\\frac{2 C_{20}}{3} \\end{array}\\right)+\\bar{I} \\mathbf{1}_{3 \\times 3}\\right)
    
    with :math:`M` the mass of the body, and :math:`R` the reference radius of the spherical harmonic coefficients. The term :math:`\\bar{I}` is the mean moment of inertia. The spherical harmonic
    coefficients in the above matrix are unnormalized.
    
    
    Parameters
    ----------
    inertia tensor : numpy.ndarray[numpy.float64[3, 3]]
        Full inertia tensor :math:`\\mathbf{I}` of the body for which spherical harmonic coefficients are to be computed.
    gravitational_parameter : float
        Gravitational parameter :math:`\\mu` of the body for which spherical harmonic coefficients are to be computed.
    reference_radius : float
        Reference radius w.r.t. which spherical harmonic coefficients are to be computed.
    output_normalized_coefficients : bool, default=True
        Boolean denoting whether the coefficients computed are normalized (if true) or unnormalized (if false)
    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Tuple of two matrices, containing the spherical harmonic coefficients :math:`{C}_{lm}` (first) and :math:`{S}_{lm}` (second) up to degree and order 2.
        The degree-two coefficients are computed from the inertia tensor, the degree-one coefficients are set to zero (and :math:`C_{00}=0`)"""

def unnormalize_spherical_harmonic_coefficients(normalized_cosine_coefficients: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], normalized_sine_coefficients: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> tuple[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]]:
    """Function to unnormalize spherical harmonic coefficients
    
    Function to unnormalize spherical harmonic coefficients, using the equations provided in the :func:`~tudatpy.gravitation.astro.legendre_normalization_factor` function.
    
    Parameters
    ----------
    normalized_cosine_coefficients : numpy.ndarray
        Matrix for which entry :math:`(i,j)` contains the normalized cosine coefficient :math:`\\bar{C}_{lm}`
    normalized_sine_coefficients : numpy.ndarray
        Matrix for which entry :math:`(i,j)` contains the normalized sine coefficient :math:`\\bar{S}_{lm}`
    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Tuple of two matrices, containing the unnormalized coefficients :math:`{C}_{lm}` (first) and :math:`{S}_{lm}` (second)"""