"""
Unit tests for light deflection correction calculations in light_deflection_correction.py
"""
import numpy as np
from numpy.linalg import norm
from data.light_deflection_correction import _calculate_light_deflection, relativistic_light_deflection
from tudatpy.constants import SPEED_OF_LIGHT
import pandas as pd

# Constants
MU_BODY = 1.327e20 # Gravitational parameter of sun to use for tests

def test_light_deflection_zero_when_collinear():
    """Test that light deflection is zero when observer, asteroid and light bending body are in-line"""
    # Define positions to be on 1 line
    gaia = np.array([1.0e11, 0.0, 0.0])
    asteroid = np.array([3.0e11, 0.0, 0.0])
    body = np.array([0.0,    0.0, 0.0])
    offset_vec = _calculate_light_deflection(gaia, asteroid, body, MU_BODY)

    np.testing.assert_allclose(offset_vec, np.zeros(3), rtol=0, atol=1e-15)


def test_light_deflection_magnitude_with_analytical_formula():
    """Test that magnitude of calculated deflection vector corresponds to analytical formula for magnitude"""
    # Some arbitrary geometry:
    gaia = np.array([1.5e11, 1e10, 0.0])
    asteroid = np.array([3.0e11, 0.0,  0.0])
    body = np.array([0.0,    0.0,  0.0])
    gaia_wrt_body = gaia - body
    asteroid_wrt_body = asteroid - body
    phase_angle = np.arccos(np.dot(gaia_wrt_body, asteroid_wrt_body) / (norm(gaia_wrt_body) * norm(asteroid_wrt_body)))

    # Offset magnitude from function:
    offset = norm(_calculate_light_deflection(gaia, asteroid, body, MU_BODY))
    # Eq 72 from Klioner (2003):
    expected_offset = 2 * MU_BODY / (norm(gaia_wrt_body) * SPEED_OF_LIGHT ** 2) * np.tan(phase_angle/2)

    assert np.isclose(offset, expected_offset, rtol=1e-10, atol=0)

def test_light_deflection_direction_properties():
    """Test properties of light deflection vector direction"""
    # Some arbitrary geometry:
    gaia = np.array([1.5e11, 1e10, 0.0])
    asteroid = np.array([3.0e11, 0.0,  0.0])
    body = np.array([0.0,    0.0,  0.0])

    offset_vec = _calculate_light_deflection(gaia, asteroid, body, MU_BODY)
    line_of_sight = gaia - asteroid

    # Check that vector is perpendicular to line of sight
    assert np.isclose(np.dot(line_of_sight, offset_vec), 0, rtol=0, atol=1e-14)

    # Check that vector points away from light deflecting body
    assert np.dot(offset_vec, gaia) > 0

    # Check that the vector is in the common plane by checking if triple product is zero.
    triple_prod = np.dot(offset_vec, np.cross(gaia, asteroid))
    assert np.allclose(triple_prod, 0.0, rtol=0, atol=1e-15)


def test_singularity_at_opposition():
    """Deflection vector becomes undef. when source and observer are at opposite sides of deflecting body"""
    # Set Gaia position at +x, asteroid at -x
    gaia = np.array([1.0e11, 0, 0.0])
    asteroid = np.array([-1.0e11, 0.0, 0.0])
    body = np.array([0.0, 0.0, 0.0])
    offset_vec = _calculate_light_deflection(gaia, asteroid, body, MU_BODY)

    assert np.all(np.isnan(offset_vec))


def test_light_deflection_integration(mocker):
    """Integration test for light deflection calculation for a simple case"""
    # Create mock for bodies
    bodies_mock = mocker.MagicMock()
    asteroid_mock = mocker.MagicMock()
    sun_mock = mocker.MagicMock()
    # What bodies.get() should return:
    to_return = lambda body_name: sun_mock if body_name == 'Sun' else asteroid_mock
    bodies_mock.get.side_effect = to_return
    asteroid_mock.ephemeris.cartesian_position.return_value = np.array([1e11, 0, 0])
    sun_mock.ephemeris.cartesian_position.return_value = np.zeros(3)
    sun_mock.gravitational_parameter = 1.327e20

    # Create mock observation table
    data = {'epoch': [0,1],
            'number_mp': [0,0],
            'ra': [np.pi/2, -np.pi/2],
            'dec': [0, 0],
            'x_gaia': [1e11, 1e11],
            'y_gaia': [-1e11, 1e11], # For second epoch, Gaia position is flipped over x-axis
            'z_gaia': [0, 0]}
    table = pd.DataFrame(data)
    corrections = relativistic_light_deflection(0, table, bodies_mock, ['Sun'])

    # We check whether hand-calculated corrections to RA/DEC match those returned by the function
    expected_offset = 8.6490584e-9 # Calculated using eq 72 from Klioner
    expected_ra_corr = [expected_offset, -expected_offset]
    expected_dec_error = [0, 0]
    expected_corrections = np.column_stack((expected_ra_corr, expected_dec_error))

    np.testing.assert_allclose(corrections[:,0], expected_corrections[:,0], rtol=1e-8, atol=0) # RA
    np.testing.assert_allclose(corrections[:,1], expected_corrections[:,1], rtol=0, atol=1e-15) # DEC



