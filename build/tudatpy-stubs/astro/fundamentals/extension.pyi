import numpy
import pybind11_stubgen.typing_ext
import typing
__all__ = ['compute_shadow_function']

def compute_shadow_function(occulted_body_position: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)], occulted_body_radius: float, occulting_body_position: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)], occulting_body_radius: float, satellite_position: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]) -> float:
    """Compute the shadow function.
    
    Returns the value of of the shadow function. Returns 0 if the satellite is in umbra, 1 if the
    satellite is fully exposed and a value between 0 and 1 if the satellite is in penumbra or antumbra.
    
    The point of view is from the satellite. The occulting body (for example the Earth) is the body
    that blocks the light from the occulted body (for example the Sun).
    
    Reference: Section 3.4 from ( Montebruck O, Gill E., 2005) and Fig. 5 from (Zhang et al., 2019).
    
    Parameters
    ----------
    occulted_body_position : numpy.ndarray
        Vector containing Cartesian coordinates of the occulted body.
    occulted_body_radius : float
        Mean radius of occulted body.
    occulting_body_position : numpy.ndarray
        Vector containing Cartesian coordinates of the occulting body.
    occulting_body_radius : float
        Mean radius of occulting body.
    satellite_position : numpy.ndarray
        Vector containing Cartesian coordinates of the satellite.
    Returns
    -------
    float
        Shadow function value"""