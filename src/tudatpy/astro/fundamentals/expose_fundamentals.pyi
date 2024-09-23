import typing
import numpy
__all__ = ['compute_shadow_function']

def compute_shadow_function(occulted_body_position: numpy.ndarray, occulted_body_radius: float, occulting_body_position: numpy.ndarray, occulting_body_radius: float, satellite_position: numpy.ndarray) -> float:
    ...