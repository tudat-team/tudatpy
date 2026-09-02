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

_MODULE = "tudatpy.estimation.observations.observation_corrections.photocenter_correction"


def _build_angular_observation_collection(observation_pairs, times, body_name, observer_body_name):
    """Build a real ObservationCollection holding angular observations of one body by one observer"""
    link_ends = {
        links.transmitter: links.body_origin_link_end_id(body_name),
        links.receiver: links.body_origin_link_end_id(observer_body_name),
    }
    observation_set = observations.create_single_observation_set(
        model_settings.angular_position_type,
        link_ends,
        list(observation_pairs),
        list(times),
        links.receiver,
    )
    return observations.ObservationCollection([observation_set])


def _unit(vector):
    return vector / norm(vector)


def _spherical_offset_magnitude(radius, solar_phase_angle):
    """Closed form spherical solution of the photocenter-barycenter offset to check against."""
    cotangent = lambda angle: np.cos(angle) / np.sin(angle)
    numerator = 2 * (
        np.sin(solar_phase_angle) + (np.pi - solar_phase_angle) * np.cos(solar_phase_angle)
    )
    denominator = (
        3
        * np.pi
        * (
            cotangent(solar_phase_angle / 2)
            - np.sin(solar_phase_angle / 2) * np.log(cotangent(solar_phase_angle / 4))
        )
    )
    return numerator / denominator * radius


@pytest.mark.parametrize("solar_phase_angle", [0.2, 0.5, 1.0, 1.5, 2.5])
def test_offset_norm_with_respect_to_polynomial(solar_phase_angle):
    """Test consistency of offset magnitude calculation with approximate polynomial from Fuentes-Munoz (2024)"""
    radius = 50.0

    # Sun and observer directions as seen from the asteroid, spanning the solar phase angle
    unit_vector_to_observer = np.array([1.0, 0.0, 0.0])
    unit_vector_to_sun = np.array([np.cos(solar_phase_angle), np.sin(solar_phase_angle), 0.0])

    # A sphere is an ellipsoid with three equal semi-axes
    photocenter_offset = _photocenter_correction_ellipsoidal(
        [radius] * 3,
        unit_vector_to_sun,
        unit_vector_to_observer,
    )
    plane_of_sky_offset_from_function = norm(
        photocenter_offset
        - np.dot(photocenter_offset, unit_vector_to_observer) * unit_vector_to_observer
    )

    approximate_offset_polynomial = (
        lambda phase_angle: -0.02384 * phase_angle**3
        + 0.05579 * phase_angle**2
        + 0.329 * phase_angle
    )  # From reference paper
    plane_of_sky_offset_from_polynomial = approximate_offset_polynomial(solar_phase_angle) * radius

    assert np.isclose(
        plane_of_sky_offset_from_polynomial,
        plane_of_sky_offset_from_function,
        rtol=1e-2,
        atol=0,
    )  # Within 1%


@pytest.mark.parametrize("in_place", [True, False])
def test_apply_photocenter_correction_to_observation_collection(in_place):
    """Test that the wrapper adds computed corrections to a real collection's observations and wraps RA"""
    # First observation RA sits just below +pi so its correction pushes it over the boundary (tests RA wrapping)
    observation_pairs = [np.array([np.pi - 1e-9, 0.2]), np.array([0.3, 0.4])]  # [RA, DEC]
    angular_corrections = np.array([[2e-9, 2e-9], [3e-9, 4e-9]])
    collection = _build_angular_observation_collection(
        observation_pairs, [0.0, 1.0], "Asteroid", "Observer"
    )

    with patch(
        f"{_MODULE}.photocenter_correction_angular_observations", return_value=angular_corrections
    ):
        result = apply_photocenter_correction_to_observation_collection(
            observation_collection=collection,
            body_dimensions=1000.0,
            bodies=MagicMock(),
            body_name="Asteroid",
            observer_body_name="Observer",
            in_place=in_place,
        )

    # Corrected observations, with RA (column 0) wrapped to (-pi, pi]
    expected_observations = np.array(observation_pairs) + angular_corrections
    expected_observations[:, 0] = (expected_observations[:, 0] + np.pi) % (2 * np.pi) - np.pi
    assert (
        expected_observations[0, 0] < 0
    )  # Wrapping mapped the RA to the negative side of the boundary

    if in_place:
        assert result is None
        target_collection = collection
    else:
        assert result is not None
        target_collection = result

    corrected_observations, _ = target_collection.get_concatenated_observations_and_times()
    np.testing.assert_allclose(
        np.array(corrected_observations).reshape(-1, 2), expected_observations
    )


def _mock_bodies(inertial_to_body_fixed=np.eye(3)):
    """Mocked SystemOfBodies with a static geometry: Sun in the origin, asteroid on the +x axis, observer mirrored
    over the x-axis between the two epochs (t < 0.5 and t >= 0.5). Both epochs have a 90 deg solar phase angle.

    Passing None for inertial_to_body_fixed gives an asteroid without a rotation model."""
    bodies_mock = MagicMock()
    sun_mock, asteroid_mock, observer_mock = MagicMock(), MagicMock(), MagicMock()
    bodies_mock.get.side_effect = lambda name: (
        {"Sun": sun_mock, "Observer": observer_mock}.get(name, asteroid_mock)
    )
    bodies_mock.does_body_exist.return_value = True

    sun_mock.state_in_base_frame_from_ephemeris.return_value = np.zeros(6)
    asteroid_mock.state_in_base_frame_from_ephemeris.return_value = np.array([1e11, 0, 0, 0, 0, 0])
    observer_mock.state_in_base_frame_from_ephemeris.side_effect = lambda epoch: (
        np.array([1e11, -1e11, 0, 0, 0, 0]) if epoch < 0.5 else np.array([1e11, 1e11, 0, 0, 0, 0])
    )

    if inertial_to_body_fixed is None:
        asteroid_mock.rotation_model = None
    else:
        # Constant rotation model, so that the body-fixed frame does not depend on the epoch
        asteroid_mock.rotation_model.inertial_to_body_fixed_rotation.side_effect = (
            lambda t: inertial_to_body_fixed
        )
    return bodies_mock


@pytest.mark.parametrize("solar_phase_angle", [0.2, 0.5, 1.0, 1.5, 2.5])
def test_offset_reduces_to_spherical_equation_for_equal_semi_axes(solar_phase_angle):
    """Test that the offset matches the closed-form spherical equation when all semi-axes equal the radius"""
    radius = 500.0

    # Sun and observer directions as seen from the asteroid, spanning the solar phase angle
    unit_vector_to_observer = np.array([1.0, 0.0, 0.0])
    unit_vector_to_sun = np.array([np.cos(solar_phase_angle), np.sin(solar_phase_angle), 0.0])

    # The offset is a 3D position in meter, of which only the plane-of-sky part shifts the astrometry
    photocenter_offset = _photocenter_correction_ellipsoidal(
        [radius] * 3,
        unit_vector_to_sun,
        unit_vector_to_observer,
    )
    photocenter_offset_plane_of_sky = (
        photocenter_offset
        - np.dot(photocenter_offset, unit_vector_to_observer) * unit_vector_to_observer
    )

    np.testing.assert_allclose(
        norm(photocenter_offset_plane_of_sky),
        _spherical_offset_magnitude(radius, solar_phase_angle),
        rtol=1e-12,
        atol=0,
    )


def test_ellipsoidal_offset_direction_properties():
    """Test properties of the direction of the ellipsoidal photocenter offset vector"""
    semi_axes = [300.0, 200.0, 150.0]
    unit_vector_to_sun = _unit(np.array([1.0, 2.0, 3.0]))
    unit_vector_to_observer = _unit(np.array([2.0, -1.0, 0.0]))
    photocenter_offset = _photocenter_correction_ellipsoidal(
        semi_axes, unit_vector_to_sun, unit_vector_to_observer
    )

    # The offset lies in the plane spanned by the Sun and observer directions
    sun_observer_plane_normal = _unit(np.cross(unit_vector_to_sun, unit_vector_to_observer))
    assert np.isclose(
        np.dot(_unit(photocenter_offset), sun_observer_plane_normal),
        0.0,
        rtol=0,
        atol=1e-14,
    )

    # Its plane-of-sky component points towards the Sun
    photocenter_offset_plane_of_sky = (
        photocenter_offset
        - np.dot(photocenter_offset, unit_vector_to_observer) * unit_vector_to_observer
    )
    sun_direction_plane_of_sky = (
        unit_vector_to_sun
        - np.dot(unit_vector_to_sun, unit_vector_to_observer) * unit_vector_to_observer
    )
    assert (
        np.dot(
            _unit(photocenter_offset_plane_of_sky),
            _unit(sun_direction_plane_of_sky),
        )
        > 0
    )

    # The photocenter lies inside the ellipsoid, so the offset cannot exceed the largest semi-axis
    assert norm(photocenter_offset) < max(semi_axes)


def test_ellipsoidal_offset_scales_linearly_with_body_size():
    """Test that scaling up the ellipsoid scales the offset by the same factor (the geometry is size-independent)"""
    semi_axes = np.array([300.0, 200.0, 150.0])
    unit_vector_to_sun = _unit(np.array([1.0, 2.0, 3.0]))
    unit_vector_to_observer = _unit(np.array([2.0, -1.0, 0.0]))
    scale = 7.0

    photocenter_offset = _photocenter_correction_ellipsoidal(
        semi_axes, unit_vector_to_sun, unit_vector_to_observer
    )
    scaled_photocenter_offset = _photocenter_correction_ellipsoidal(
        scale * semi_axes,
        unit_vector_to_sun,
        unit_vector_to_observer,
    )

    # Absolute tolerance in meter, needed for components that are (near) zero
    np.testing.assert_allclose(
        scaled_photocenter_offset,
        scale * photocenter_offset,
        rtol=1e-12,
        atol=1e-12,
    )


def test_ellipsoidal_offset_symmetries():
    """Test that the offset transforms consistently under relabelling and mirroring of the ellipsoid axes"""
    semi_axes = np.array([300.0, 200.0, 150.0])
    unit_vector_to_sun = _unit(np.array([1.0, 2.0, 3.0]))
    unit_vector_to_observer = _unit(np.array([2.0, -1.0, 0.0]))
    photocenter_offset = _photocenter_correction_ellipsoidal(
        semi_axes, unit_vector_to_sun, unit_vector_to_observer
    )

    # Swapping the a and b semi-axes (and the x/y components of the input directions) swaps them in the offset too
    axis_swap_matrix = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    photocenter_offset_with_swapped_axes = _photocenter_correction_ellipsoidal(
        axis_swap_matrix @ semi_axes,
        axis_swap_matrix @ unit_vector_to_sun,
        axis_swap_matrix @ unit_vector_to_observer,
    )
    np.testing.assert_allclose(
        photocenter_offset_with_swapped_axes,
        axis_swap_matrix @ photocenter_offset,
        rtol=1e-12,
        atol=1e-14,
    )

    # The ellipsoid is symmetric about the ab plane, so mirroring both directions in z mirrors the offset
    axis_mirror_matrix = np.diag([1.0, 1.0, -1.0])
    photocenter_offset_with_mirrored_axes = _photocenter_correction_ellipsoidal(
        semi_axes,
        axis_mirror_matrix @ unit_vector_to_sun,
        axis_mirror_matrix @ unit_vector_to_observer,
    )
    np.testing.assert_allclose(
        photocenter_offset_with_mirrored_axes,
        axis_mirror_matrix @ photocenter_offset,
        rtol=1e-12,
        atol=1e-14,
    )


def test_ellipsoidal_offset_difference_line_of_sight_direction():
    """Test that the photocenter offset is larger when the long axis is oriented in the direction of the sun,
    than when it is oriented in the direction of the observer line-of-sight.
    """
    semi_axes = [
        1000.0,
        100.0,
        100.0,
    ]  # Strongly stretched ellipsoid, with the long axis along the body-fixed x axis

    # Project 3D offset vector along plane-of-sky
    def plane_of_sky_offset(
        unit_vector_to_sun,
        unit_vector_to_observer,
    ):
        photocenter_offset = _photocenter_correction_ellipsoidal(
            semi_axes,
            unit_vector_to_sun,
            unit_vector_to_observer,
        )
        return norm(
            photocenter_offset
            - np.dot(photocenter_offset, unit_vector_to_observer) * unit_vector_to_observer
        )

    # Long axis perpendicular to line of sight. Expected offset large
    offset_broadside = plane_of_sky_offset(
        unit_vector_to_sun=np.array([1, 0.0, 0.0]),
        unit_vector_to_observer=np.array([0.0, 0.0, 1.0]),
    )

    # Long axis along the line of sight: expected offset small
    offset_end_on = plane_of_sky_offset(
        unit_vector_to_sun=np.array([0.0, 1.0, 0.0]),
        unit_vector_to_observer=np.array([1.0, 0.0, 0.0]),
    )

    assert offset_broadside > offset_end_on


def test_ellipsoidal_offset_singularity_at_zero_phase_angle():
    """Test that a zero solar phase angle is a singularity, as for the spherical approximation"""
    common_sun_observer_direction = _unit(np.array([1.0, 0.5, -0.3]))

    with np.errstate(divide="ignore", invalid="ignore"):
        photocenter_offset = _photocenter_correction_ellipsoidal(
            [300.0, 200.0, 150.0],
            common_sun_observer_direction,
            common_sun_observer_direction,
        )

    assert np.all(np.isnan(photocenter_offset))


def test_spherical_corrections_integration():
    """Integration test for a spherical body, whose result can be verified by hand. A sphere has no orientation, so
    this also tests that the correction is computed without the body having a rotation model."""
    radius = 500.0

    # Observations [time, RA, DEC] consistent with the mocked geometry
    observations = np.array([[0.0, np.pi / 2, 0.0], [1.0, -np.pi / 2, 0.0]])

    angular_corrections = photocenter_correction_angular_observations(
        observations=observations,
        body_dimensions=radius,  # A single number, so the body is a sphere
        bodies=_mock_bodies(None),  # No rotation model
        body_name="Asteroid",
        observer_body_name="Observer",
    )
    right_ascension_corrections, declination_corrections = angular_corrections.T

    # Hand-calculated result (solar phase angle of 90 deg at both epochs)
    expected_right_ascension_correction = 2.8160935e-9
    np.testing.assert_allclose(declination_corrections, np.zeros(2), rtol=0, atol=1e-15)
    np.testing.assert_allclose(
        right_ascension_corrections,
        np.array(
            (
                -expected_right_ascension_correction,
                expected_right_ascension_correction,
            )
        ),
        rtol=1e-7,
        atol=0,
    )


def test_ellipsoidal_corrections_use_rotation_at_emission_time():
    """Test that the body orientation is evaluated at the time of emission, i.e. one light time before reception"""
    bodies_mock = _mock_bodies()
    observations = np.array([[0.0, np.pi / 2, 0.0]])

    photocenter_correction_angular_observations(
        observations=observations,
        body_dimensions=[300.0, 200.0, 150.0],
        bodies=bodies_mock,
        body_name="Asteroid",
        observer_body_name="Observer",
    )

    # Observer and asteroid are 1e11 m apart in the mocked geometry
    expected_epoch = 0.0 - 1e11 / SPEED_OF_LIGHT
    rotation_model = bodies_mock.get("Asteroid").rotation_model
    rotation_model.inertial_to_body_fixed_rotation.assert_called_once_with(expected_epoch)


def test_ellipsoidal_input_validation():
    """Test that invalid observation shapes and semi-axis counts are rejected"""
    bodies_mock = _mock_bodies()
    observations = np.array([[0.0, np.pi / 2, 0.0]])

    # Observations not in N x 3 shape
    with pytest.raises(ValueError, match="N x 3"):
        photocenter_correction_angular_observations(
            observations[:, :2], [300.0, 200.0, 150.0], bodies_mock, "Asteroid", "Observer"
        )

    # Semi-axes not of length 3
    with pytest.raises(ValueError, match="length 3"):
        photocenter_correction_angular_observations(
            observations, [300.0, 200.0], bodies_mock, "Asteroid", "Observer"
        )


def test_ellipsoid_requires_rotation_model():
    """Test that a genuine ellipsoid is rejected without a rotation model, unlike a sphere"""
    observations = np.array([[0.0, np.pi / 2, 0.0]])

    with pytest.raises(ValueError, match="rotation model"):
        photocenter_correction_angular_observations(
            observations, [300.0, 200.0, 150.0], _mock_bodies(None), "Asteroid", "Observer"
        )
