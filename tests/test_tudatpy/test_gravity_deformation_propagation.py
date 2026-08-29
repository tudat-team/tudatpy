import numpy as np

from tudatpy.astro import gravitation
from tudatpy.dynamics import environment_setup, propagation_setup, simulator


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
    moon_gravity_settings.scaled_mean_moment_of_inertia = 0.4
    body_settings.get("Moon").gravity_field_settings = moon_gravity_settings
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
    solver_settings = propagation_setup.propagator.CoupledStateDerivativeSolverSettings(
        relative_tolerance=2.0e-11,
        absolute_scaled_tolerance=2.0e-13,
        maximum_iterations=30,
    )
    propagator_settings.coupled_state_derivative_solver_settings = solver_settings

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
