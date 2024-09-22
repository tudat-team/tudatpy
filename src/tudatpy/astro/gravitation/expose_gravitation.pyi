import typing
import numpy
__all__ = ['legendre_normalization_factor', 'normalize_spherical_harmonic_coefficients', 'spherical_harmonic_coefficients_from_inertia', 'unnormalize_spherical_harmonic_coefficients']

def legendre_normalization_factor(degree: int, order: int) -> float:
    ...

def normalize_spherical_harmonic_coefficients(unnormalized_cosine_coefficients: numpy.ndarray, unnormalized_sine_coefficients: numpy.ndarray) -> tuple[numpy.ndarray, numpy.ndarray]:
    ...

def spherical_harmonic_coefficients_from_inertia(inertia_tensor: numpy.ndarray, gravitational_parameter: float, reference_radius: float, output_normalized_coefficients: bool=True) -> tuple[numpy.ndarray, numpy.ndarray, float]:
    ...

def unnormalize_spherical_harmonic_coefficients(normalized_cosine_coefficients: numpy.ndarray, normalized_sine_coefficients: numpy.ndarray) -> tuple[numpy.ndarray, numpy.ndarray]:
    ...