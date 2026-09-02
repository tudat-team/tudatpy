"""
Unit tests for light deflection correction calculations in light_deflection_correction.py
"""

from unittest.mock import MagicMock, patch

import pytest
import numpy as np
from numpy.linalg import norm
from tudatpy.estimation import observations
from tudatpy.estimation.observable_models_setup import model_settings, links
from tudatpy.estimation.observations.observation_corrections.light_deflection_correction import (
    _light_deflection_single_contribution,
    light_deflection_correction_angular_observations,
    apply_light_deflection_correction_to_observation_collection,
)
from tudatpy.constants import SPEED_OF_LIGHT

# Constants
MU_SUN = 1.327e20  # Gravitational parameter of sun to use for tests
_MODULE = "tudatpy.estimation.observations.observation_corrections.light_deflection_correction"


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


def test_light_deflection_zero_when_collinear():
    """Test that light deflection is zero when observer, asteroid and light bending body are in-line"""
    # Define positions to be on 1 line
    observer_position = np.array([1.0e11, 0.0, 0.0])
    observed_body_position = np.array([3.0e11, 0.0, 0.0])
    perturbing_body_position = np.array([0.0, 0.0, 0.0])
    light_deflection_offset = _light_deflection_single_contribution(
        observer_position,
        observed_body_position,
        perturbing_body_position,
        MU_SUN,
    )

    np.testing.assert_allclose(light_deflection_offset, np.zeros(3), rtol=0, atol=1e-15)


def test_light_deflection_magnitude_with_analytical_formula():
    """Test that magnitude of calculated deflection vector corresponds to analytical formula for magnitude"""
    # Some arbitrary geometry:
    observer_position = np.array([1.5e11, 1e10, 0.0])
    observed_body_position = np.array([3.0e11, 0.0, 0.0])
    perturbing_body_position = np.array([0.0, 0.0, 0.0])
    position_observer_with_respect_to_perturbing_body = observer_position - perturbing_body_position
    position_observed_body_with_respect_to_perturbing_body = (
        observed_body_position - perturbing_body_position
    )
    phase_angle = np.arccos(
        np.dot(
            position_observer_with_respect_to_perturbing_body,
            position_observed_body_with_respect_to_perturbing_body,
        )
        / (
            norm(position_observer_with_respect_to_perturbing_body)
            * norm(position_observed_body_with_respect_to_perturbing_body)
        )
    )

    # Offset magnitude from function:
    light_deflection_offset_magnitude = norm(
        _light_deflection_single_contribution(
            observer_position,
            observed_body_position,
            perturbing_body_position,
            MU_SUN,
        )
    )
    # Eq 72 from Klioner (2003):
    expected_light_deflection_offset_magnitude = (
        2
        * MU_SUN
        / (norm(position_observer_with_respect_to_perturbing_body) * SPEED_OF_LIGHT**2)
        * np.tan(phase_angle / 2)
    )

    assert np.isclose(
        light_deflection_offset_magnitude,
        expected_light_deflection_offset_magnitude,
        rtol=1e-10,
        atol=0,
    )


def test_light_deflection_direction_properties():
    """Test properties of light deflection vector direction"""
    # Some arbitrary geometry:
    observer_position = np.array([1.5e11, 1e10, 0.0])
    observed_body_position = np.array([3.0e11, 0.0, 0.0])
    perturbing_body_position = np.array([0.0, 0.0, 0.0])

    light_deflection_offset = _light_deflection_single_contribution(
        observer_position,
        observed_body_position,
        perturbing_body_position,
        MU_SUN,
    )
    line_of_sight = observer_position - observed_body_position

    # Check that vector is perpendicular to line of sight
    assert np.isclose(np.dot(line_of_sight, light_deflection_offset), 0, rtol=0, atol=1e-14)

    # Check that vector points away from light deflecting body
    assert np.dot(light_deflection_offset, observer_position) > 0

    # Check that the vector is in the common plane by checking if triple product is zero.
    common_plane_triple_product = np.dot(
        light_deflection_offset,
        np.cross(observer_position, observed_body_position),
    )
    assert np.allclose(common_plane_triple_product, 0.0, rtol=0, atol=1e-15)


def test_singularity_at_opposition():
    """Deflection vector becomes undef. when source and observer are at opposite sides of deflecting body"""
    # Set observer position at +x, asteroid at -x
    observer_position = np.array([1.0e11, 0, 0.0])
    observed_body_position = np.array([-1.0e11, 0.0, 0.0])
    perturbing_body_position = np.array([0.0, 0.0, 0.0])
    with np.errstate(divide="ignore", invalid="ignore"):
        light_deflection_offset = _light_deflection_single_contribution(
            observer_position,
            observed_body_position,
            perturbing_body_position,
            MU_SUN,
        )

    assert np.all(np.isnan(light_deflection_offset))


def test_light_deflection_integration():
    """Integration test for light deflection calculation for a simple case"""
    # Create mocks for bodies
    bodies_mock = MagicMock()
    asteroid_mock = MagicMock()
    sun_mock = MagicMock()
    observer_mock = MagicMock()
    # What bodies.get() should return:
    select_body_mock = lambda body_name: {"Sun": sun_mock, "Observer": observer_mock}.get(
        body_name, asteroid_mock
    )
    bodies_mock.get.side_effect = select_body_mock
    bodies_mock.does_body_exist.return_value = True

    # Static geometry: asteroid on +x axis, Sun in origin (positions constant, so light-time iteration has no effect)
    asteroid_mock.state_in_base_frame_from_ephemeris.return_value = np.array([1e11, 0, 0, 0, 0, 0])
    sun_mock.state_in_base_frame_from_ephemeris.return_value = np.zeros(6)
    sun_mock.gravitational_parameter = MU_SUN

    # Observer position is mirrored over the x-axis between the two epochs
    observer_mock.state_in_base_frame_from_ephemeris.side_effect = lambda epoch: (
        np.array([1e11, -1e11, 0, 0, 0, 0]) if epoch < 0.5 else np.array([1e11, 1e11, 0, 0, 0, 0])
    )

    # Observations [time, RA, DEC] consistent with the mocked geometry
    observations = np.array([[0.0, np.pi / 2, 0.0], [1.0, -np.pi / 2, 0.0]])

    angular_corrections = light_deflection_correction_angular_observations(
        observations=observations,
        bodies=bodies_mock,
        body_name="Asteroid",
        observer_body_name="Observer",
        perturbing_bodies_list=["Sun"],
    )

    # We check whether hand-calculated corrections to RA/DEC match those returned by the function
    # Observer-Sun-asteroid phase angle is 45 deg at both epochs; eq 72 from Klioner (2003):
    observer_with_respect_to_sun_distance = norm(np.array([1e11, 1e11, 0]))
    expected_light_deflection_offset = (
        2 * MU_SUN / (observer_with_respect_to_sun_distance * SPEED_OF_LIGHT**2) * np.tan(np.pi / 8)
    )
    expected_right_ascension_corrections = [
        expected_light_deflection_offset,
        -expected_light_deflection_offset,
    ]
    expected_declination_corrections = [0, 0]

    np.testing.assert_allclose(
        angular_corrections[:, 0],
        expected_right_ascension_corrections,
        rtol=1e-8,
        atol=0,
    )
    np.testing.assert_allclose(
        angular_corrections[:, 1],
        expected_declination_corrections,
        rtol=0,
        atol=1e-15,
    )


@pytest.mark.parametrize("in_place", [True, False])
def test_apply_light_deflection_correction_to_observation_collection(in_place):
    """Test that the wrapper adds computed corrections to a real collection's observations and wraps RA"""
    # First observation RA sits just below +pi so its correction pushes it over the boundary (tests RA wrapping)
    observation_pairs = [np.array([np.pi - 1e-9, 0.2]), np.array([0.3, 0.4])]  # [RA, DEC]
    angular_corrections = np.array([[2e-9, 2e-9], [3e-9, 4e-9]])
    collection = _build_angular_observation_collection(
        observation_pairs, [0.0, 1.0], "Asteroid", "Observer"
    )

    with patch(
        f"{_MODULE}.light_deflection_correction_angular_observations",
        return_value=angular_corrections,
    ):
        result = apply_light_deflection_correction_to_observation_collection(
            observation_collection=collection,
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
