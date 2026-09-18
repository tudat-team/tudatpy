import numpy as np
import pytest

from tudatpy.astro import gravitation
from tudatpy.dynamics import environment, environment_setup, propagation_setup, simulator


def test_gravity_deformation_propagation_and_environment_lifecycle():
    body_settings = environment_setup.BodyListSettings("SSB", "J2000")
    body_settings.add_empty_settings("Earth")
    body_settings.add_empty_settings("Moon")

    body_settings.get("Earth").ephemeris_settings = environment_setup.ephemeris.constant(
        np.array([1.0e5, 0.0, 0.0, 0.0, 0.0, 0.0]), "SSB", "J2000"
    )
    body_settings.get("Earth").gravity_field_settings = environment_setup.gravity_field.central(
        4.0e6
    )

    body_settings.get("Moon").ephemeris_settings = environment_setup.ephemeris.constant(
        np.zeros(6), "SSB", "J2000"
    )
    body_settings.get("Moon").rotation_model_settings = (
        environment_setup.rotation_model.constant_rotation_model(
            "J2000", "MoonFixed", np.identity(3)
        )
    )
    nominal_cosine = np.zeros((3, 3))
    nominal_sine = np.zeros((3, 3))
    nominal_cosine[0, 0] = 1.0
    nominal_cosine[2, 0] = -1.0e-3
    nominal_cosine[2, 2] = 2.0e-4
    moon_gravity_settings = environment_setup.gravity_field.spherical_harmonic(
        4.0e5,
        1.0e3,
        nominal_cosine,
        nominal_sine,
        "MoonFixed",
    )
    body_settings.get("Moon").gravity_field_settings = moon_gravity_settings
    body_settings.get("Moon").rigid_body_settings = environment_setup.rigid_body.from_gravity_field(
        0.4
    )
    bodies = environment_setup.create_system_of_bodies(body_settings)

    deformation_setting = propagation_setup.propagator.maxwell_deformation(
        maxwell_relaxation_time=20.0,
        global_relaxation_time=100.0,
        love_number=0.3,
        maximum_degree=2,
        maximum_order=2,
        perturbing_body="Earth",
    )
    deformation_models = propagation_setup.create_gravity_deformation_models(
        bodies, {"Moon": [deformation_setting]}
    )
    integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(
        1.0, coefficient_set=propagation_setup.integrator.CoefficientSets.rk_4
    )
    propagator_settings = propagation_setup.propagator.gravity_deformation(
        ["Moon"],
        deformation_models,
        np.zeros(5),
        0.0,
        integrator_settings,
        propagation_setup.propagator.time_termination(20.0),
    )
    propagator_settings.processing_settings.set_integrated_result = True

    dynamics_simulator = simulator.create_dynamics_simulator(bodies, propagator_settings)
    state_history = dynamics_simulator.propagation_results.state_history
    assert len(state_history) == 21
    final_variation = state_history[max(state_history)]
    assert np.linalg.norm(final_variation) > 0.0

    correction_cosine = np.zeros((3, 3))
    correction_sine = np.zeros((3, 3))
    correction_cosine[2, :] = final_variation[:3]
    correction_sine[2, 1:] = final_variation[3:]
    normalized_correction = gravitation.normalize_spherical_harmonic_coefficients(
        correction_cosine, correction_sine
    )
    moon = bodies.get("Moon")
    np.testing.assert_allclose(
        moon.gravity_field_model.cosine_coefficients,
        nominal_cosine + normalized_correction[0],
        rtol=2.0e-12,
        atol=1.0e-18,
    )
    np.testing.assert_allclose(
        moon.gravity_field_model.sine_coefficients,
        nominal_sine + normalized_correction[1],
        rtol=2.0e-12,
        atol=1.0e-18,
    )
    assert moon.rigid_body_properties.inertia_tensor_available
    assert moon.rigid_body_properties.inertia_tensor_derivative_available
    assert np.all(np.isfinite(moon.rigid_body_properties.current_inertia_tensor))
    assert np.all(np.isfinite(moon.rigid_body_properties.current_inertia_tensor_derivative))
    assert not bodies.get("Earth").rigid_body_properties.inertia_tensor_available

    # Integration internals are intentionally not part of the public mutation API.
    assert not hasattr(moon, "set_angular_velocity_derivative_in_local_frame")
    assert not hasattr(moon.rigid_body_properties, "update_inertia_tensor_derivative")


def test_scaled_mean_moment_settings_compatibility_and_precedence():
    nominal_cosine = np.zeros((3, 3))
    nominal_sine = np.zeros((3, 3))
    nominal_cosine[0, 0] = 1.0
    nominal_cosine[1, 0] = 1.0e-5
    nominal_cosine[1, 1] = -2.0e-5
    nominal_sine[1, 1] = 3.0e-5
    nominal_cosine[2, 0] = -1.0e-3
    nominal_cosine[2, 2] = 2.0e-4

    def gravity_settings():
        return environment_setup.gravity_field.spherical_harmonic(
            4.0e5, 1.0e3, nominal_cosine, nominal_sine, "J2000"
        )

    body_settings = environment_setup.BodyListSettings("SSB", "J2000")
    for name in ("Legacy", "Canonical", "NoInertia", "Explicit"):
        body_settings.add_empty_settings(name)
        body_settings.get(name).gravity_field_settings = gravity_settings()

    legacy_gravity_settings = body_settings.get("Legacy").gravity_field_settings
    with pytest.warns(DeprecationWarning, match="rigid_body.from_gravity_field"):
        legacy_gravity_settings.scaled_mean_moment_of_inertia = 0.4
    with pytest.warns(DeprecationWarning, match="rigid_body.from_gravity_field"):
        assert legacy_gravity_settings.scaled_mean_moment_of_inertia == 0.4

    body_settings.get("Canonical").rigid_body_settings = (
        environment_setup.rigid_body.from_gravity_field(0.4)
    )
    explicit_center_of_mass = np.array([1.0, 2.0, 3.0])
    explicit_inertia = 7.0 * np.identity(3)
    body_settings.get("Explicit").rigid_body_settings = (
        environment_setup.rigid_body.constant_rigid_body_properties(
            123.0, explicit_center_of_mass, explicit_inertia
        )
    )
    explicit_gravity_settings = body_settings.get("Explicit").gravity_field_settings
    with pytest.warns(DeprecationWarning):
        explicit_gravity_settings.scaled_mean_moment_of_inertia = 0.4

    bodies = environment_setup.create_system_of_bodies(body_settings)
    legacy = bodies.get("Legacy")
    canonical = bodies.get("Canonical")
    no_inertia = bodies.get("NoInertia")
    explicit = bodies.get("Explicit")

    np.testing.assert_allclose(
        legacy.rigid_body_properties.current_center_of_mass,
        canonical.rigid_body_properties.current_center_of_mass,
    )
    np.testing.assert_allclose(
        legacy.rigid_body_properties.current_inertia_tensor,
        canonical.rigid_body_properties.current_inertia_tensor,
    )
    assert legacy.rigid_body_properties.current_mass == pytest.approx(
        canonical.rigid_body_properties.current_mass
    )
    assert not no_inertia.rigid_body_properties.inertia_tensor_available
    assert explicit.rigid_body_properties.current_mass == pytest.approx(123.0)
    np.testing.assert_allclose(
        explicit.rigid_body_properties.current_center_of_mass,
        explicit_center_of_mass,
    )
    np.testing.assert_allclose(
        explicit.rigid_body_properties.current_inertia_tensor, explicit_inertia
    )

    # The compatibility value is an input carrier only; no runtime gravity object owns it.
    assert not hasattr(canonical.gravity_field_model, "scaled_mean_moment_of_inertia")
    assert not hasattr(canonical.gravity_field_model, "inertia_tensor")


def test_deprecated_gravity_field_update_callback_remains_functional():
    callback_calls = []
    with pytest.warns(DeprecationWarning, match="update_inertia_tensor"):
        gravity_field = environment.GravityFieldModel(
            1.0, update_inertia_tensor=lambda: callback_calls.append(True)
        )

    gravity_field.gravitational_parameter = 2.0
    assert callback_calls == [True]
