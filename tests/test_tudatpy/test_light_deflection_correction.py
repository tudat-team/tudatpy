"""
Unit tests for light deflection correction calculations in light_deflection_correction.py
"""
from unittest.mock import MagicMock

import numpy as np
from numpy.linalg import norm
from tudatpy.estimation.observable_models_setup.biases.light_deflection_correction import (
    _calculate_light_deflection,
    relativistic_light_deflection_from_observations,
)
from tudatpy.constants import SPEED_OF_LIGHT

# Constants
MU_SUN = 1.327e20 # Gravitational parameter of sun to use for tests

def test_light_deflection_zero_when_collinear():
    """Test that light deflection is zero when observer, asteroid and light bending body are in-line"""
    # Define positions to be on 1 line
    observer = np.array([1.0e11, 0.0, 0.0])
    asteroid = np.array([3.0e11, 0.0, 0.0])
    perturber = np.array([0.0,    0.0, 0.0])
    offset_vec = _calculate_light_deflection(observer, asteroid, perturber, MU_SUN)

    np.testing.assert_allclose(offset_vec, np.zeros(3), rtol=0, atol=1e-15)


def test_light_deflection_magnitude_with_analytical_formula():
    """Test that magnitude of calculated deflection vector corresponds to analytical formula for magnitude"""
    # Some arbitrary geometry:
    observer = np.array([1.5e11, 1e10, 0.0])
    asteroid = np.array([3.0e11, 0.0,  0.0])
    perturber = np.array([0.0,    0.0,  0.0])
    observer_wrt_perturber = observer - perturber
    asteroid_wrt_perturber = asteroid - perturber
    phase_angle = np.arccos(np.dot(observer_wrt_perturber, asteroid_wrt_perturber)
                            / (norm(observer_wrt_perturber) * norm(asteroid_wrt_perturber)))

    # Offset magnitude from function:
    offset = norm(_calculate_light_deflection(observer, asteroid, perturber, MU_SUN))
    # Eq 72 from Klioner (2003):
    expected_offset = 2 * MU_SUN / (norm(observer_wrt_perturber) * SPEED_OF_LIGHT ** 2) * np.tan(phase_angle/2)

    assert np.isclose(offset, expected_offset, rtol=1e-10, atol=0)

def test_light_deflection_direction_properties():
    """Test properties of light deflection vector direction"""
    # Some arbitrary geometry:
    observer = np.array([1.5e11, 1e10, 0.0])
    asteroid = np.array([3.0e11, 0.0,  0.0])
    perturber = np.array([0.0,    0.0,  0.0])

    offset_vec = _calculate_light_deflection(observer, asteroid, perturber, MU_SUN)
    line_of_sight = observer - asteroid

    # Check that vector is perpendicular to line of sight
    assert np.isclose(np.dot(line_of_sight, offset_vec), 0, rtol=0, atol=1e-14)

    # Check that vector points away from light deflecting body
    assert np.dot(offset_vec, observer) > 0

    # Check that the vector is in the common plane by checking if triple product is zero.
    triple_prod = np.dot(offset_vec, np.cross(observer, asteroid))
    assert np.allclose(triple_prod, 0.0, rtol=0, atol=1e-15)


def test_singularity_at_opposition():
    """Deflection vector becomes undef. when source and observer are at opposite sides of deflecting body"""
    # Set observer position at +x, asteroid at -x
    observer = np.array([1.0e11, 0, 0.0])
    asteroid = np.array([-1.0e11, 0.0, 0.0])
    perturber = np.array([0.0, 0.0, 0.0])
    with np.errstate(divide='ignore', invalid='ignore'):
        offset_vec = _calculate_light_deflection(observer, asteroid, perturber, MU_SUN)

    assert np.all(np.isnan(offset_vec))


def test_light_deflection_integration():
    """Integration test for light deflection calculation for a simple case"""
    # Create mocks for bodies
    bodies_mock = MagicMock()
    asteroid_mock = MagicMock()
    sun_mock = MagicMock()
    observer_mock = MagicMock()
    # What bodies.get() should return:
    to_return = lambda body_name: {'Sun': sun_mock, 'Observer': observer_mock}.get(body_name, asteroid_mock)
    bodies_mock.get.side_effect = to_return
    bodies_mock.does_body_exist.return_value = True

    # Static geometry: asteroid on +x axis, Sun in origin (positions constant, so light-time iteration has no effect)
    asteroid_mock.state_in_base_frame_from_ephemeris.return_value = np.array([1e11, 0, 0, 0, 0, 0])
    sun_mock.state_in_base_frame_from_ephemeris.return_value = np.zeros(6)
    sun_mock.gravitational_parameter = MU_SUN

    # Observer position is mirrored over the x-axis between the two epochs
    observer_mock.state_in_base_frame_from_ephemeris.side_effect = lambda epoch: (
        np.array([1e11, -1e11, 0, 0, 0, 0]) if epoch < 0.5 else np.array([1e11, 1e11, 0, 0, 0, 0]))

    # Observations [time, RA, DEC] consistent with the mocked geometry
    observations = np.array([[0.0,  np.pi / 2, 0.0],
                             [1.0, -np.pi / 2, 0.0]])

    corrections = relativistic_light_deflection_from_observations(
        observations = observations,
        bodies = bodies_mock,
        body_name = 'Asteroid',
        observer_body_name = 'Observer',
        perturbing_bodies_list = ['Sun'],
    )

    # We check whether hand-calculated corrections to RA/DEC match those returned by the function
    # Observer-Sun-asteroid phase angle is 45 deg at both epochs; eq 72 from Klioner (2003):
    observer_wrt_sun_norm = norm(np.array([1e11, 1e11, 0]))
    expected_offset = 2 * MU_SUN / (observer_wrt_sun_norm * SPEED_OF_LIGHT ** 2) * np.tan(np.pi / 8)
    expected_ra_corr = [expected_offset, -expected_offset]
    expected_dec_corr = [0, 0]

    np.testing.assert_allclose(corrections[:, 0], expected_ra_corr, rtol=1e-8, atol=0) # RA
    np.testing.assert_allclose(corrections[:, 1], expected_dec_corr, rtol=0, atol=1e-15) # DEC
