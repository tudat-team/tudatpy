"""
Utility functions for observation corrections
"""
import numpy as np
from numpy.linalg import norm

def _unit(vector:np.ndarray) -> np.ndarray:
    return vector / norm(vector)


def _offset_vector_to_corrections(offset_vec: np.ndarray,
                                       ra: float,
                                       dec: float):
    """
    Calculate corrections in right ascension and declination from a plane-of-sky offset vector.

    Args:
        offset_vec (np.ndarray): Offset vector such that true dir + offset = observed dir.
        ra: Right ascension (rad)
        dec: Declination (rad)

    Returns:
        float: correction in right ascension
        float: correction in declination
    """
    observed_dir = np.array([np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)])
    # true dir. + offset = observed dir.
    true_dir = observed_dir - offset_vec # Small angle approximation
    true_dir = _unit(true_dir)

    ra_true = np.arctan2(true_dir[1], true_dir[0])
    dec_true = np.arctan2(true_dir[2], np.sqrt(true_dir[0] ** 2 + true_dir[1] ** 2))

    ra_corr = ((ra_true - ra) + np.pi) % (2 * np.pi) - np.pi # Wraps to correct range
    dec_corr = dec_true - dec

    return ra_corr, dec_corr
