"""
Utility functions for observation corrections
"""

import numpy as np
from numpy.linalg import norm


def _unit(vector: np.ndarray) -> np.ndarray:
    return vector / norm(vector)


def _offset_vector_to_corrections(
    plane_of_sky_offset: np.ndarray,
    right_ascension: float,
    declination: float,
):
    """
    Convert a plane-of-sky offset vector, i.e. the offset that is perpendicular to the line-of-sight, to
    corrections in right ascension and declination.
    """
    observed_direction = np.array(
        [
            np.cos(right_ascension) * np.cos(declination),
            np.sin(right_ascension) * np.cos(declination),
            np.sin(declination),
        ]
    )
    # True direction + plane-of-sky offset = observed direction.
    true_direction = observed_direction - plane_of_sky_offset  # Small-angle approximation
    direction_norm = norm(true_direction)
    if not np.isfinite(direction_norm) or direction_norm == 0:
        raise ValueError("Observation correction must produce a finite, nonzero direction")
    true_direction /= direction_norm

    true_right_ascension = np.arctan2(true_direction[1], true_direction[0])
    true_declination = np.arctan2(
        true_direction[2],
        np.sqrt(true_direction[0] ** 2 + true_direction[1] ** 2),
    )

    right_ascension_correction = ((true_right_ascension - right_ascension) + np.pi) % (
        2 * np.pi
    ) - np.pi
    declination_correction = true_declination - declination

    return right_ascension_correction, declination_correction
