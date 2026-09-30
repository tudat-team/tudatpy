"""
Unit tests for observation correction utilities in _correction_utils.py
"""

from tudatpy.estimation.observations.observation_corrections._correction_utils import (
    _offset_vector_to_corrections,
)
import numpy as np


def test_corrections_zero_offset_gives_zero():
    """Test that zero offset vector gives zero corrections to RA/DEC"""
    plane_of_sky_offset = np.zeros(3)
    right_ascension_correction, declination_correction = _offset_vector_to_corrections(
        plane_of_sky_offset,
        right_ascension=0.5,
        declination=0.3,
    )

    assert np.isclose(right_ascension_correction, 0, rtol=0, atol=1e-15)
    assert np.isclose(declination_correction, 0, rtol=0, atol=1e-15)


def test_corrections_pure_ra_direction():
    """Test correction calculation for offset vector in RA direction"""
    # Point offset vector in direction of purely increasing RA coordinate
    right_ascension, declination = 0.0, 0.0
    angular_offset = 1e-9  # mas level
    plane_of_sky_offset = np.array([0.0, angular_offset, 0.0])
    right_ascension_correction, declination_correction = _offset_vector_to_corrections(
        plane_of_sky_offset, right_ascension, declination
    )

    # Check that RA is changed and DEC unchanged
    assert np.isclose(declination_correction, 0.0, rtol=0, atol=1e-15)
    # By small angle approximation, correction to RA is equal to offset in y-axis
    # (atol floor accounts for floating point noise of the RA wrapping)
    assert np.isclose(right_ascension_correction, -angular_offset, rtol=1e-10, atol=1e-15)


def test_corrections_pure_dec_direction():
    """Test correction calculation for offset vector in DEC direction"""
    # Point offset vector in direction of purely increasing DEC coordinate
    right_ascension, declination = 0.0, 0.0
    angular_offset = 1e-9  # mas level
    plane_of_sky_offset = np.array([0.0, 0.0, angular_offset])
    right_ascension_correction, declination_correction = _offset_vector_to_corrections(
        plane_of_sky_offset, right_ascension, declination
    )

    # Check that RA is unchanged and DEC changed
    assert np.isclose(declination_correction, -angular_offset, rtol=1e-10, atol=0)
    # By small angle approximation, correction to Dec is equal to offset in z-axis
    assert np.isclose(right_ascension_correction, 0.0, rtol=0, atol=1e-15)


def test_corrections_ra_and_dec():
    """That that RA/DEC corrections are correctly calculated for some arbitrary offset vector"""
    right_ascension, declination = np.pi / 4, np.pi / 6

    # Get local unit vectors of increasing RA and DEC
    increasing_right_ascension_direction = np.array(
        [-np.sin(right_ascension), np.cos(right_ascension), 0.0]
    )
    increasing_declination_direction = np.array(
        [
            -np.sin(declination) * np.cos(right_ascension),
            -np.sin(declination) * np.sin(right_ascension),
            np.cos(declination),
        ]
    )

    # Build an offset with known projections along the two angular-coordinate directions.
    plane_of_sky_offset = (
        2e-9 * increasing_right_ascension_direction + 3e-9 * increasing_declination_direction
    )
    right_ascension_correction, declination_correction = _offset_vector_to_corrections(
        plane_of_sky_offset, right_ascension, declination
    )

    assert np.isclose(
        right_ascension_correction,
        -2e-9 / np.cos(declination),
        rtol=1e-7,
        atol=1e-15,
    )
    assert np.isclose(declination_correction, -3e-9, rtol=1e-7, atol=0)


def test_corrections_ra_wrapping():
    """Test that RA corrections are wrapped to [-pi, pi) across the RA discontinuity"""
    # Observation just below +pi, offset pushing the true direction across the discontinuity
    right_ascension, declination = np.pi - 1e-9, 0.0
    angular_offset = 4e-9
    increasing_right_ascension_direction = np.array(
        [-np.sin(right_ascension), np.cos(right_ascension), 0.0]
    )
    plane_of_sky_offset = (
        -angular_offset * increasing_right_ascension_direction
    )  # True direction ends up past +pi
    right_ascension_correction, declination_correction = _offset_vector_to_corrections(
        plane_of_sky_offset, right_ascension, declination
    )

    # Without wrapping the correction would be ~ -2*pi + 4e-9
    assert np.isclose(right_ascension_correction, angular_offset, rtol=1e-7, atol=0)
    assert np.isclose(declination_correction, 0.0, rtol=0, atol=1e-15)
