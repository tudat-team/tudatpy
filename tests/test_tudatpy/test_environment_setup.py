import pytest
import numpy as np

from tudatpy.astro import time_representation
from tudatpy.dynamics import environment, environment_setup


def test_deprecated_mass_properties_alias_preserves_behavior(capfd):
    body_settings = environment_setup.BodyListSettings("SSB", "J2000")
    body_settings.add_empty_settings("Canonical")
    body_settings.add_empty_settings("Legacy")
    bodies = environment_setup.create_system_of_bodies(body_settings)
    settings = environment_setup.rigid_body.constant_rigid_body_properties(125.0)

    environment_setup.add_rigid_body_properties(
        bodies=bodies,
        body_name="Canonical",
        rigid_body_property_settings=settings,
    )
    assert capfd.readouterr().err == ""

    environment_setup.add_mass_properties_model(
        bodies=bodies,
        body_name="Legacy",
        mass_property_settings=settings,
    )
    warning = capfd.readouterr().err
    assert "add_mass_properties_model is deprecated" in warning
    assert "tudatpy.dynamics.environment_setup.add_rigid_body_properties" in warning
    assert bodies.get("Legacy").mass == pytest.approx(bodies.get("Canonical").mass)
    assert bodies.get("Legacy").mass == pytest.approx(125.0)


def test_custom_ephemeris_time_object_preserves_sub_double_epochs():
    """Preserve distinct callback and Body states for epochs sharing one double value."""
    left = time_representation.Time(765201600)
    right = left + 20e-9
    velocity = np.array([30000.0, 10.0, -5.0])
    received = []

    def state(epoch):
        """Evaluate linear motion from a precise local time difference."""
        assert isinstance(epoch, time_representation.Time)
        received.append(epoch)

        return np.r_[velocity * float(epoch - left), velocity]

    custom = environment.CustomEphemerisTimeObject(state, "SSB", "J2000")
    settings = environment_setup.BodyListSettings("SSB", "J2000")
    settings.add_empty_settings("Probe")
    bodies = environment_setup.create_system_of_bodies(settings)
    body = bodies.get("Probe")
    body.ephemeris = custom
    bodies.add_body(body, "Probe")
    expected = np.r_[velocity * float(right - left), velocity]

    # Both Python epoch overloads and the body's native extended-Time dispatch reach the callback.
    assert float(left) == float(right)
    np.testing.assert_allclose(custom.cartesian_state(left).reshape(6)[:3], 0.0, atol=0.0)
    np.testing.assert_allclose(
        custom.cartesian_state(right).reshape(6), expected, rtol=0, atol=1e-15
    )
    np.testing.assert_allclose(
        body.state_in_base_frame_from_ephemeris(left).reshape(6)[:3], 0.0, atol=0.0
    )
    np.testing.assert_allclose(
        body.state_in_base_frame_from_ephemeris(right).reshape(6), expected, rtol=0, atol=1e-15
    )

    # The retained float route deliberately receives the already rounded epoch.
    np.testing.assert_allclose(custom.cartesian_state(float(right)).reshape(6)[:3], 0.0, atol=0.0)
    assert any(epoch == right for epoch in received)


def test_earth_rotation_time_object_preserves_sub_double_epochs():
    """Resolve Earth rotation over an interval lost when both epochs become doubles."""

    # Use the native Earth model with epochs separated below the double time resolution.
    settings = environment_setup.BodyListSettings("SSB", "GCRS")
    settings.add_empty_settings("Earth")
    settings.get("Earth").rotation_model_settings = environment_setup.rotation_model.gcrs_to_itrs()
    rotation = environment_setup.create_system_of_bodies(settings).get("Earth").rotation_model
    left = time_representation.Time(765201600)
    right = left + 20e-9
    assert float(left) == float(right)

    # Resolve station motion over the precise interval; the native ERA rounding floor is audited separately.
    station = np.array([6378137.0, 0.0, 0.0])
    left_matrix = rotation.body_fixed_to_inertial_rotation(left)
    right_matrix = rotation.body_fixed_to_inertial_rotation(right)
    derivative = rotation.time_derivative_body_fixed_to_inertial_rotation(left)
    displacement = (right_matrix - left_matrix) @ station
    assert np.linalg.norm(displacement) > 8e-6

    # Retain float behavior and verify both inverse rotation and derivative overloads.
    np.testing.assert_array_equal(
        rotation.body_fixed_to_inertial_rotation(float(left)),
        rotation.body_fixed_to_inertial_rotation(float(right)),
    )
    for epoch, matrix in [(left, left_matrix), (right, right_matrix)]:
        np.testing.assert_allclose(
            rotation.inertial_to_body_fixed_rotation(epoch), matrix.T, rtol=0, atol=1e-15
        )
        np.testing.assert_allclose(
            rotation.time_derivative_inertial_to_body_fixed_rotation(epoch),
            rotation.time_derivative_body_fixed_to_inertial_rotation(epoch).T,
            rtol=0,
            atol=1e-18,
        )
    assert (
        np.linalg.norm(rotation.time_derivative_body_fixed_to_inertial_rotation(right) - derivative)
        > 0.0
    )
