"""
Unit tests for photocenter correction calculations in photocenter_correction.py
"""
from unittest.mock import MagicMock, patch

import pytest
import numpy as np
from numpy.linalg import norm
from tudatpy.estimation import observations
from tudatpy.estimation.observable_models_setup import model_settings, links
from tudatpy.estimation.observable_models_setup.biases.photocenter_correction import (
    _photocenter_offset,
    photocenter_corrections_from_observations,
    apply_photocenter_correction_to_observation_collection,
)

_MODULE = 'tudatpy.estimation.observable_models_setup.biases.photocenter_correction'


def _build_angular_observation_collection(observation_pairs, times, body_name, observer_body_name):
    """Build a real ObservationCollection holding angular observations of one body by one observer"""
    link_ends = {
        links.transmitter: links.body_origin_link_end_id(body_name),
        links.receiver: links.body_origin_link_end_id(observer_body_name),
    }
    observation_set = observations.create_single_observation_set(
        model_settings.angular_position_type, link_ends, list(observation_pairs), list(times), links.receiver)
    return observations.ObservationCollection([observation_set])

def _unit(vec):
    return vec / norm(vec)

@pytest.mark.parametrize("solar_phase_angle", [0.2, 0.5, 1.0, 1.5, 2.5])
def test_photocenter_offset_norm_wrt_polynomial(solar_phase_angle):
    """Test consistency of offset magnitude calculation with approximate polynomial from Fuentes-Munoz (2024)"""
    # Set auxiliary parameters (irrelevant for test)
    diameter = 100
    observer_asteroid_distance = 1e10

    # Geometry with the requested solar phase angle
    asteroid_wrt_observer = observer_asteroid_distance * np.array([1.0, 0.0, 0.0])
    asteroid_wrt_sun = 2e11 * np.array([np.cos(solar_phase_angle), np.sin(solar_phase_angle), 0.0])

    # Offset magnitudes using polynomial approximation and true magnitude
    poly = lambda x: -0.02384*x**3 + 0.05579*x**2 + 0.329*x # From reference paper
    offset_from_polynomial = poly(solar_phase_angle) * (diameter/2) / observer_asteroid_distance
    offset_from_function = norm(_photocenter_offset(diameter, asteroid_wrt_observer, asteroid_wrt_sun))

    assert np.isclose(offset_from_polynomial, offset_from_function, rtol=1e-2, atol=0) # Within 1%


def test_photocenter_offset_singularity_at_zero_phase_angle():
    """Test that a zero solar phase angle is a singularity"""
    # Sun and observer directions aligned -> zero solar phase angle
    asteroid_wrt_observer = np.array([1.0e10, 0.0, 0.0])
    asteroid_wrt_sun = np.array([2.0e11, 0.0, 0.0])

    with np.errstate(divide='ignore', invalid='ignore'):
        offset_vec = _photocenter_offset(100.0, asteroid_wrt_observer, asteroid_wrt_sun)

    assert np.all(np.isnan(offset_vec))


def test_offset_direction_properties():
    """Test properties of the direction of the photocenter offset vector"""
    asteroid_wrt_observer = np.array([2.0e11, 0.0, 0.0])
    asteroid_wrt_sun = np.array([1.0e11, 1.0e11, 0.0])
    offset_dir = _unit(_photocenter_offset(100.0, asteroid_wrt_observer, asteroid_wrt_sun))
    line_of_sight = _unit(asteroid_wrt_observer)

    # Test that the offset is perpendicular to the line of sight
    assert np.isclose(np.dot(offset_dir, line_of_sight), 0.0, rtol=0, atol=1e-15)

    # Test that the vector is in the common plane spanned by asteroid, observer, sun
    triple_prod = np.dot(offset_dir, _unit(np.cross(asteroid_wrt_observer, asteroid_wrt_sun)))
    assert np.isclose(triple_prod, 0.0, rtol=0, atol=1e-15)

    # Test that the vector points in the direction of the sun
    observer_wrt_sun = - asteroid_wrt_observer + asteroid_wrt_sun
    assert np.dot(observer_wrt_sun, offset_dir) < 0


def test_photocenter_corrections_integration():
    """Integration test for photocenter correction calculation"""
    # Create mocks for bodies object
    bodies_mock = MagicMock()
    sun_mock = MagicMock()
    asteroid_mock = MagicMock()
    observer_mock = MagicMock()
    to_return = lambda body_name: {'Sun': sun_mock, 'Observer': observer_mock}.get(body_name, asteroid_mock)
    bodies_mock.get.side_effect = to_return
    bodies_mock.does_body_exist.return_value = True

    # Static geometry: asteroid on +x axis, Sun in origin
    asteroid_mock.state_in_base_frame_from_ephemeris.return_value = np.array([1e11, 0, 0, 0, 0, 0])
    sun_mock.state_in_base_frame_from_ephemeris.return_value = np.zeros(6)

    # Observer position is mirrored over the x-axis between the two epochs
    observer_mock.state_in_base_frame_from_ephemeris.side_effect = lambda epoch: (
        np.array([1e11, -1e11, 0, 0, 0, 0]) if epoch < 0.5 else np.array([1e11, 1e11, 0, 0, 0, 0]))

    # Observations [time, RA, DEC] consistent with the mocked geometry
    observations = np.array([[0.0,  np.pi / 2, 0.0],
                             [1.0, -np.pi / 2, 0.0]])

    corrections = photocenter_corrections_from_observations(
        observations = observations,
        diameter = 1000.0,
        bodies = bodies_mock,
        body_name = 'Asteroid',
        observer_body_name = 'Observer',
    )
    ra_corr, dec_corr = corrections.T

    # Compare against results calculated by hand (solar phase angle of 90 deg at both epochs)
    true_ra_corr = 2.8160935e-9
    np.testing.assert_allclose(dec_corr, np.zeros(2), rtol=0, atol=1e-15)
    np.testing.assert_allclose(ra_corr, np.array((-true_ra_corr, true_ra_corr)), rtol=1e-7, atol=0)


def test_input_validation():
    """Test that invalid observation shapes and diameters are rejected"""
    bodies_mock = MagicMock()
    bodies_mock.does_body_exist.return_value = True
    observations = np.array([[0.0, np.pi / 2, 0.0]])

    # Observations not in N x 3 shape
    with pytest.raises(ValueError, match='shape N x 3'):
        photocenter_corrections_from_observations(observations[:, :2], 1000.0, bodies_mock, 'Asteroid', 'Observer')

    # Non-positive diameter
    with pytest.raises(ValueError, match='diameter'):
        photocenter_corrections_from_observations(observations, -1.0, bodies_mock, 'Asteroid', 'Observer')


@pytest.mark.parametrize('in_place', [True, False])
def test_apply_photocenter_correction_to_observation_collection(in_place):
    """Test that the wrapper adds computed corrections to a real collection's observations and wraps RA"""
    # First observation RA sits just below +pi so its correction pushes it over the boundary (tests RA wrapping)
    observation_pairs = [np.array([np.pi - 1e-9, 0.2]), np.array([0.3, 0.4])] # [RA, DEC]
    corrections = np.array([[2e-9, 2e-9], [3e-9, 4e-9]])
    collection = _build_angular_observation_collection(observation_pairs, [0.0, 1.0], 'Asteroid', 'Observer')

    with patch(f'{_MODULE}.photocenter_corrections_from_observations', return_value=corrections):
        result = apply_photocenter_correction_to_observation_collection(
            observation_collection=collection,
            diameter=1000.0,
            bodies=MagicMock(),
            body_name='Asteroid',
            observer_body_name='Observer',
            in_place=in_place,
        )

    # Corrected observations, with RA (column 0) wrapped to (-pi, pi]
    expected = np.array(observation_pairs) + corrections
    expected[:, 0] = (expected[:, 0] + np.pi) % (2 * np.pi) - np.pi
    assert expected[0, 0] < 0 # Wrapping mapped the RA to the negative side of the boundary

    if in_place:
        assert result is None
        target = collection
    else:
        assert result is not None
        target = result

    corrected, _ = target.get_concatenated_observations_and_times()
    np.testing.assert_allclose(np.array(corrected).reshape(-1, 2), expected)
