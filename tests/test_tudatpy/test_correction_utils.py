"""
Unit tests for observation correction utilities in observation_correction_utils.py
"""
from tudatpy.estimation.observable_models_setup.biases.observation_correction_utils import _offset_vector_to_corrections
import numpy as np

def test_corrections_zero_offset_gives_zero():
    """Test that zero offset vector gives zero corrections to RA/DEC"""
    offset_vec = np.zeros(3)
    ra_corr, dec_corr = _offset_vector_to_corrections(offset_vec, ra=0.5, dec=0.3)

    assert np.isclose(ra_corr, 0, rtol=0, atol=1e-15)
    assert np.isclose(dec_corr, 0, rtol=0, atol=1e-15)


def test_corrections_pure_ra_direction():
    """Test correction calculation for offset vector in RA direction"""
    # Point offset vector in direction of purely increasing RA coordinate
    ra, dec = 0.0, 0.0
    correction = 1e-9 # mas level
    offset_vec = np.array([0.0, correction, 0.0])
    ra_corr, dec_corr = _offset_vector_to_corrections(offset_vec, ra, dec)

    # Check that RA is changed and DEC unchanged
    assert np.isclose(dec_corr, 0.0, rtol=0, atol=1e-15)
    # By small angle approximation, correction to RA is equal to offset in y-axis
    # (atol floor accounts for floating point noise of the RA wrapping)
    assert np.isclose(ra_corr, -correction, rtol=1e-10, atol=1e-15)


def test_corrections_pure_dec_direction():
    """Test correction calculation for offset vector in DEC direction"""
    # Point offset vector in direction of purely increasing DEC coordinate
    ra, dec = 0.0, 0.0
    correction = 1e-9  # mas level
    offset_vec = np.array([0.0, 0.0, correction])
    ra_corr, dec_corr = _offset_vector_to_corrections(offset_vec, ra, dec)

    # Check that RA is unchanged and DEC changed
    assert np.isclose(dec_corr, -correction, rtol=1e-10, atol=0)
    # By small angle approximation, correction to Dec is equal to offset in z-axis
    assert np.isclose(ra_corr, 0.0, rtol=0, atol=1e-15)


def test_corrections_ra_and_dec():
    """That that RA/DEC corrections are correctly calculated for some arbitrary offset vector"""
    ra, dec = np.pi / 4, np.pi / 6

    # Get local unit vectors of increasing RA and DEC
    dir_ra = np.array([-np.sin(ra), np.cos(ra), 0.0])
    dir_dec = np.array([-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)])

    # Build an offset with known projections: 2 along dir_ra, 3 along dir_dec
    offset_vec = 2e-9 * dir_ra + 3e-9 * dir_dec
    ra_corr, dec_corr = _offset_vector_to_corrections(offset_vec, ra, dec)

    assert np.isclose(ra_corr, -2e-9/np.cos(dec), rtol=1e-7, atol=1e-15)
    assert np.isclose(dec_corr, -3e-9, rtol=1e-7, atol=0)


def test_corrections_ra_wrapping():
    """Test that RA corrections are wrapped to [-pi, pi) across the RA discontinuity"""
    # Observation just below +pi, offset pushing the true direction across the discontinuity
    ra, dec = np.pi - 1e-9, 0.0
    correction = 4e-9
    dir_ra = np.array([-np.sin(ra), np.cos(ra), 0.0])
    offset_vec = -correction * dir_ra # True direction ends up past +pi
    ra_corr, dec_corr = _offset_vector_to_corrections(offset_vec, ra, dec)

    # Without wrapping the correction would be ~ -2*pi + 4e-9
    assert np.isclose(ra_corr, correction, rtol=1e-7, atol=0)
    assert np.isclose(dec_corr, 0.0, rtol=0, atol=1e-15)
