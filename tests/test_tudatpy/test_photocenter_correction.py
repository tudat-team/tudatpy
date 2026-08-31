"""
Unit tests for photocenter correction calculations in photocenter_correction.py
"""
from unittest.mock import MagicMock, patch

import pytest
import numpy as np
from numpy.linalg import norm
from tudatpy.constants import SPEED_OF_LIGHT
from tudatpy.estimation import observations
from tudatpy.estimation.observable_models_setup import model_settings, links
from tudatpy.estimation.observations.observation_corrections.photocenter_correction import (
    _photocenter_correction_ellipsoidal,
    photocenter_correction_angular_observations,
    apply_photocenter_correction_to_observation_collection,
)

_MODULE = 'tudatpy.estimation.observations.observation_corrections.photocenter_correction'


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


def _spherical_offset_magnitude(radius, solar_phase_angle):
    """Closed form spherical solution of the photocenter-barycenter offset to check against."""
    cot = lambda x: np.cos(x) / np.sin(x)
    numerator = 2 * (np.sin(solar_phase_angle) + (np.pi - solar_phase_angle) * np.cos(solar_phase_angle))
    denominator = 3 * np.pi * (cot(solar_phase_angle / 2)
                               - np.sin(solar_phase_angle / 2) * np.log(cot(solar_phase_angle / 4)))
    return numerator / denominator * radius


@pytest.mark.parametrize("solar_phase_angle", [0.2, 0.5, 1.0, 1.5, 2.5])
def test_offset_norm_wrt_polynomial(solar_phase_angle):
    """Test consistency of offset magnitude calculation with approximate polynomial from Fuentes-Munoz (2024)"""
    radius = 50.0

    # Sun and observer directions as seen from the asteroid, spanning the solar phase angle
    e_observer = np.array([1.0, 0.0, 0.0])
    e_sun = np.array([np.cos(solar_phase_angle), np.sin(solar_phase_angle), 0.0])

    # A sphere is an ellipsoid with three equal semi-axes
    offset = _photocenter_correction_ellipsoidal([radius] * 3, e_sun, e_observer)
    offset_from_function = norm(offset - np.dot(offset, e_observer) * e_observer)

    poly = lambda x: -0.02384*x**3 + 0.05579*x**2 + 0.329*x # From reference paper
    offset_from_polynomial = poly(solar_phase_angle) * radius

    assert np.isclose(offset_from_polynomial, offset_from_function, rtol=1e-2, atol=0) # Within 1%


@pytest.mark.parametrize('in_place', [True, False])
def test_apply_photocenter_correction_to_observation_collection(in_place):
    """Test that the wrapper adds computed corrections to a real collection's observations and wraps RA"""
    # First observation RA sits just below +pi so its correction pushes it over the boundary (tests RA wrapping)
    observation_pairs = [np.array([np.pi - 1e-9, 0.2]), np.array([0.3, 0.4])] # [RA, DEC]
    corrections = np.array([[2e-9, 2e-9], [3e-9, 4e-9]])
    collection = _build_angular_observation_collection(observation_pairs, [0.0, 1.0], 'Asteroid', 'Observer')

    with patch(f'{_MODULE}.photocenter_correction_angular_observations', return_value=corrections):
        result = apply_photocenter_correction_to_observation_collection(
            observation_collection=collection,
            body_dimensions=1000.0,
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


def _mock_bodies(inertial_to_body_fixed=np.eye(3)):
    """Mocked SystemOfBodies with a static geometry: Sun in the origin, asteroid on the +x axis, observer mirrored
    over the x-axis between the two epochs (t < 0.5 and t >= 0.5). Both epochs have a 90 deg solar phase angle.

    Passing None for inertial_to_body_fixed gives an asteroid without a rotation model."""
    bodies_mock = MagicMock()
    sun_mock, asteroid_mock, observer_mock = MagicMock(), MagicMock(), MagicMock()
    bodies_mock.get.side_effect = lambda name: (
        {'Sun': sun_mock, 'Observer': observer_mock}.get(name, asteroid_mock))
    bodies_mock.does_body_exist.return_value = True

    sun_mock.state_in_base_frame_from_ephemeris.return_value = np.zeros(6)
    asteroid_mock.state_in_base_frame_from_ephemeris.return_value = np.array([1e11, 0, 0, 0, 0, 0])
    observer_mock.state_in_base_frame_from_ephemeris.side_effect = lambda epoch: (
        np.array([1e11, -1e11, 0, 0, 0, 0]) if epoch < 0.5 else np.array([1e11, 1e11, 0, 0, 0, 0]))

    if inertial_to_body_fixed is None:
        asteroid_mock.rotation_model = None
    else:
        # Constant rotation model, so that the body-fixed frame does not depend on the epoch
        asteroid_mock.rotation_model.inertial_to_body_fixed_rotation.side_effect = lambda t: inertial_to_body_fixed
    return bodies_mock


@pytest.mark.parametrize("solar_phase_angle", [0.2, 0.5, 1.0, 1.5, 2.5])
def test_offset_reduces_to_spherical_equation_for_equal_semi_axes(solar_phase_angle):
    """Test that the offset matches the closed-form spherical equation when all semi-axes equal the radius"""
    radius = 500.0

    # Sun and observer directions as seen from the asteroid, spanning the solar phase angle
    e_observer = np.array([1.0, 0.0, 0.0])
    e_sun = np.array([np.cos(solar_phase_angle), np.sin(solar_phase_angle), 0.0])

    # The offset is a 3D position in meter, of which only the plane-of-sky part shifts the astrometry
    offset = _photocenter_correction_ellipsoidal([radius] * 3, e_sun, e_observer)
    offset_sky = offset - np.dot(offset, e_observer) * e_observer

    np.testing.assert_allclose(norm(offset_sky), _spherical_offset_magnitude(radius, solar_phase_angle),
                               rtol=1e-12, atol=0)


def test_ellipsoidal_offset_direction_properties():
    """Test properties of the direction of the ellipsoidal photocenter offset vector"""
    semi_axes = [300.0, 200.0, 150.0]
    e_sun, e_observer = _unit(np.array([1.0, 2.0, 3.0])), _unit(np.array([2.0, -1.0, 0.0]))
    offset = _photocenter_correction_ellipsoidal(semi_axes, e_sun, e_observer)

    # The offset lies in the plane spanned by the Sun and observer directions
    plane_normal = _unit(np.cross(e_sun, e_observer))
    assert np.isclose(np.dot(_unit(offset), plane_normal), 0.0, rtol=0, atol=1e-14)

    # Its plane-of-sky component points towards the Sun
    offset_sky = offset - np.dot(offset, e_observer) * e_observer
    sun_sky = e_sun - np.dot(e_sun, e_observer) * e_observer
    assert np.dot(_unit(offset_sky), _unit(sun_sky)) > 0

    # The photocenter lies inside the ellipsoid, so the offset cannot exceed the largest semi-axis
    assert norm(offset) < max(semi_axes)


def test_ellipsoidal_offset_scales_linearly_with_body_size():
    """Test that scaling up the ellipsoid scales the offset by the same factor (the geometry is size-independent)"""
    semi_axes = np.array([300.0, 200.0, 150.0])
    e_sun, e_observer = _unit(np.array([1.0, 2.0, 3.0])), _unit(np.array([2.0, -1.0, 0.0]))
    scale = 7.0

    offset = _photocenter_correction_ellipsoidal(semi_axes, e_sun, e_observer)
    offset_scaled = _photocenter_correction_ellipsoidal(scale * semi_axes, e_sun, e_observer)

    # Absolute tolerance in meter, needed for components that are (near) zero
    np.testing.assert_allclose(offset_scaled, scale * offset, rtol=1e-12, atol=1e-12)


def test_ellipsoidal_offset_symmetries():
    """Test that the offset transforms consistently under relabelling and mirroring of the ellipsoid axes"""
    semi_axes = np.array([300.0, 200.0, 150.0])
    e_sun, e_observer = _unit(np.array([1.0, 2.0, 3.0])), _unit(np.array([2.0, -1.0, 0.0]))
    offset = _photocenter_correction_ellipsoidal(semi_axes, e_sun, e_observer)

    # Swapping the a and b semi-axes (and the x/y components of the input directions) swaps them in the offset too
    swap = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    offset_swapped = _photocenter_correction_ellipsoidal(swap @ semi_axes, swap @ e_sun, swap @ e_observer)
    np.testing.assert_allclose(offset_swapped, swap @ offset, rtol=1e-12, atol=1e-14)

    # The ellipsoid is symmetric about the ab plane, so mirroring both directions in z mirrors the offset
    mirror = np.diag([1.0, 1.0, -1.0])
    offset_mirrored = _photocenter_correction_ellipsoidal(semi_axes, mirror @ e_sun, mirror @ e_observer)
    np.testing.assert_allclose(offset_mirrored, mirror @ offset, rtol=1e-12, atol=1e-14)


def test_ellipsoidal_offset_difference_line_of_sight_direction():
    """Test that the photocenter offset is larger when the long axis is oriented in the direction of the sun,
    than when it is oriented in the direction of the observer line-of-sight.
    """
    semi_axes = [1000.0, 100.0, 100.0] # Strongly stretched ellipsoid, with the long axis along the body-fixed x axis

    # Project 3D offset vector along plane-of-sky
    def plane_of_sky_offset(e_sun, e_observer):
        offset = _photocenter_correction_ellipsoidal(semi_axes, e_sun, e_observer)
        return norm(offset - np.dot(offset, e_observer) * e_observer)

    # Long axis perpendicular to line of sight. Expected offset large
    offset_broadside = plane_of_sky_offset(
        e_sun=np.array([1, 0.0, 0.0]),
        e_observer=np.array([0.0, 0.0, 1.0]))

    # Long axis along the line of sight: expected offset small
    offset_end_on = plane_of_sky_offset(
        e_sun=np.array([0.0, 1.0, 0.0]),
        e_observer=np.array([1.0, 0.0, 0.0]))

    assert offset_broadside > offset_end_on


def test_ellipsoidal_offset_singularity_at_zero_phase_angle():
    """Test that a zero solar phase angle is a singularity, as for the spherical approximation"""
    e_common = _unit(np.array([1.0, 0.5, -0.3]))

    with np.errstate(divide='ignore', invalid='ignore'):
        offset = _photocenter_correction_ellipsoidal([300.0, 200.0, 150.0], e_common, e_common)

    assert np.all(np.isnan(offset))


def test_spherical_corrections_integration():
    """Integration test for a spherical body, whose result can be verified by hand. A sphere has no orientation, so
    this also tests that the correction is computed without the body having a rotation model."""
    radius = 500.0

    # Observations [time, RA, DEC] consistent with the mocked geometry
    observations = np.array([[0.0,  np.pi / 2, 0.0],
                             [1.0, -np.pi / 2, 0.0]])

    corrections = photocenter_correction_angular_observations(
        observations=observations,
        body_dimensions=radius, # A single number, so the body is a sphere
        bodies=_mock_bodies(None), # No rotation model
        body_name='Asteroid',
        observer_body_name='Observer',
    )
    ra_corr, dec_corr = corrections.T

    # Hand-calculated result (solar phase angle of 90 deg at both epochs)
    true_ra_corr = 2.8160935e-9
    np.testing.assert_allclose(dec_corr, np.zeros(2), rtol=0, atol=1e-15)
    np.testing.assert_allclose(ra_corr, np.array((-true_ra_corr, true_ra_corr)), rtol=1e-7, atol=0)


def test_ellipsoidal_corrections_use_rotation_at_emission_time():
    """Test that the body orientation is evaluated at the time of emission, i.e. one light time before reception"""
    bodies_mock = _mock_bodies()
    observations = np.array([[0.0, np.pi / 2, 0.0]])

    photocenter_correction_angular_observations(
        observations=observations, body_dimensions=[300.0, 200.0, 150.0], bodies=bodies_mock,
        body_name='Asteroid', observer_body_name='Observer')

    # Observer and asteroid are 1e11 m apart in the mocked geometry
    expected_epoch = 0.0 - 1e11 / SPEED_OF_LIGHT
    rotation_model = bodies_mock.get('Asteroid').rotation_model
    rotation_model.inertial_to_body_fixed_rotation.assert_called_once_with(expected_epoch)


def test_ellipsoidal_input_validation():
    """Test that invalid observation shapes and semi-axis counts are rejected"""
    bodies_mock = _mock_bodies()
    observations = np.array([[0.0, np.pi / 2, 0.0]])

    # Observations not in N x 3 shape
    with pytest.raises(ValueError, match='N x 3'):
        photocenter_correction_angular_observations(
            observations[:, :2], [300.0, 200.0, 150.0], bodies_mock, 'Asteroid', 'Observer')

    # Semi-axes not of length 3
    with pytest.raises(ValueError, match='length 3'):
        photocenter_correction_angular_observations(
            observations, [300.0, 200.0], bodies_mock, 'Asteroid', 'Observer')


def test_ellipsoid_requires_rotation_model():
    """Test that a genuine ellipsoid is rejected without a rotation model, unlike a sphere"""
    observations = np.array([[0.0, np.pi / 2, 0.0]])

    with pytest.raises(ValueError, match='rotation model'):
        photocenter_correction_angular_observations(
            observations, [300.0, 200.0, 150.0], _mock_bodies(None), 'Asteroid', 'Observer')
