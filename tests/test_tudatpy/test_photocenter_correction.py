"""
Unit tests for photocenter correction calculations in photocenter_correction.py
"""
import pytest
from data.photocenter_correction import (_solar_phase_angle,
                                         _offset_direction,
                                         _offset_magnitude,
                                         _offset_vector_to_corrections,
                                         retrieve_asteroid_diameter,
                                         photocenter_offset_spherical)
import numpy as np
import pandas as pd

def _unit(vec):
    return vec / np.linalg.norm(vec)

@pytest.mark.parametrize("solar_phase_angle", [0.2, 0.5, 1.0, 1.5, 2.5])
def test_photocenter_offset_norm_wrt_polynomial(solar_phase_angle):
    """Test consistency of offset magnitude calculation with approximate polynomial from Fuentes-Munoz (2024)"""
    # Set auxiliary parameters (irrelevant for test)
    diameter = 100
    gaia_asteroid_distance = 1e10

    # Offset magnitudes using polynomial approximation and true magnitude
    poly = lambda x: -0.02384*x**3 + 0.05579*x**2 + 0.329*x # From reference paper
    offset_from_polynomial = poly(solar_phase_angle) * (diameter/2) / gaia_asteroid_distance
    offset_from_function = _offset_magnitude(solar_phase_angle, diameter, gaia_asteroid_distance)

    assert np.isclose(offset_from_polynomial, offset_from_function, rtol=1e-2, atol=0) # Within 1%


def test_photocenter_offset_norm_singularity():
    """Test that a zero solar phase angle is a singularity"""
    solar_phase_angle = 0.0

    with pytest.raises(AssertionError):
        _ = _offset_magnitude(solar_phase_angle, 1., 1.)


@pytest.mark.parametrize('mpc_number, diameter_from_lookup',[(433, 16.84),(799, 47.185),(3557, 23.60858589)])
def test_diameter_retrieval(mpc_number, diameter_from_lookup):
    """Test that diameter from function is consistent with that from online lookup"""
    # Note diameter for 3557 manually calculated from its magnitude H
    diameter_from_func = retrieve_asteroid_diameter(mpc_number) / 1000 # Convert to km

    assert np.isclose(diameter_from_lookup, diameter_from_func, rtol=1e-10, atol=0)


def test_solar_phase_angle_calculation():
    """Test that solar phase angle is calculated correctly for arbitrary case"""
    # 45 degree angle between sun and observer directions
    asteroid_wrt_gaia = np.array([1.0e11, 0.0, 0.0])
    asteroid_wrt_sun = np.array([1.0e11, 1.0e11, 0.0])

    solar_phase_angle = _solar_phase_angle(_unit(asteroid_wrt_gaia), _unit(asteroid_wrt_sun))

    assert np.isclose(solar_phase_angle, np.pi/4, rtol=1e-10, atol=0)


def test_offset_direction_properties():
    """Test properties of the direction of the photocenter offset vector"""
    asteroid_wrt_gaia = np.array([2.0e11, 0.0, 0.0])
    asteroid_wrt_sun = np.array([1.0e11, 1.0e11, 0.0])
    offset_dir = _offset_direction(_unit(asteroid_wrt_gaia), _unit(asteroid_wrt_sun))
    line_of_sight = _unit(asteroid_wrt_gaia)

    # Test that the offset is perpendicular to the line of sight
    assert np.isclose(np.dot(offset_dir, line_of_sight), 0.0, rtol=0, atol=1e-15)

    # Test that the vector is in the common plane spanned by asteroid, Gaia, sun
    triple_prod = np.dot(offset_dir, np.cross(asteroid_wrt_gaia, asteroid_wrt_sun))
    assert np.isclose(triple_prod, 0.0, rtol=0, atol=1e-15)

    # Test that the vector points in the direction of the sun
    gaia_wrt_sun = - asteroid_wrt_gaia + asteroid_wrt_sun
    assert np.dot(gaia_wrt_sun, offset_dir) < 0


def test_photocenter_offset_integration(mocker):
    """Integration test for photocenter offset calculation"""
    # Create mock for bodies object
    bodies_mock = mocker.MagicMock()
    sun_mock = mocker.MagicMock()
    asteroid_mock = mocker.MagicMock()
    func = lambda body_name: sun_mock if body_name == "Sun" else asteroid_mock
    bodies_mock.get.side_effect = func
    asteroid_mock.ephemeris.cartesian_position.side_effect = [np.array([1e11, 0, 0])] * 2
    sun_mock.ephemeris.cartesian_position.side_effect = [np.zeros(3)] * 2

    # Patch diameter
    mocker.patch('Code.data.photocenter_correction._retrieve_asteroid_diameter', return_value=1000.0)

    # Create mock table
    # We create 3 observations with 2 transits. The correction for observation 1 and 2 should be identical
    # even though their ra/dec values are different
    data = {'epoch': [0, 1, 2],
            'number_mp': [0, 0, 0],
            'transit_id': [0, 0, 1],
            'x_gaia': [1e11, 1e11, 1e11],
            'y_gaia': [-1e11, -1e11, 1e11],
            'z_gaia': [0, 0, 0],
            'ra': [np.pi/2, np.pi/2, -np.pi/2],
            'dec': [0, np.pi, 0]}
    table_mock = pd.DataFrame(data)

    ra_corr, dec_corr = photocenter_offset_spherical(0, table_mock, bodies_mock).T
    # Compare against results calculated by hand
    true_ra_corr = 2.8160935e-9
    np.testing.assert_allclose(dec_corr, np.zeros(3), rtol=0, atol=1e-15)
    np.testing.assert_allclose(ra_corr, np.array((-true_ra_corr, -true_ra_corr, true_ra_corr)), rtol=1e-7, atol=0)





