"""Regression coverage for multi-body order-invariant k Love number setup (Issue #734).

Uses the locally built kernel bindings so the test exercises the C++ parameter path
without depending on optional pure-Python ephemeris wrappers.
"""

import numpy as np

from tudatpy.kernel.dynamics import environment_setup, parameters_setup, propagation_setup
from tudatpy.kernel.dynamics.parameters_setup import EstimatableParameterTypes
from tudatpy.kernel.estimation import estimation_analysis
from tudatpy.kernel.interface import spice


def _create_multi_body_tide_environment():
    spice.load_standard_kernels()

    bodies_to_create = ["Earth", "Moon", "Sun"]
    body_settings = environment_setup.get_default_body_settings(bodies_to_create, "Earth", "J2000")
    body_settings.get("Earth").gravity_field_variation_settings = [
        environment_setup.gravity_field_variation.solid_body_tide("Moon", 0.301, 2),
        environment_setup.gravity_field_variation.solid_body_tide("Sun", 0.295, 2),
    ]
    bodies = environment_setup.create_system_of_bodies(body_settings)

    bodies.create_empty_body("Delfi")
    acceleration_settings = {
        "Delfi": {"Earth": [propagation_setup.acceleration.spherical_harmonic_gravity(2, 2)]}
    }
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, acceleration_settings, ["Delfi"], ["Earth"]
    )

    initial_state = np.array([7.0e6, 0.0, 0.0, 0.0, 7.5e3, 0.0])
    termination = propagation_setup.propagator.time_termination(3600.0)
    integrator = propagation_setup.integrator.runge_kutta_fixed_step(
        100.0, propagation_setup.integrator.CoefficientSets.rkf_78
    )
    propagator_settings = propagation_setup.propagator.translational(
        ["Earth"],
        acceleration_models,
        ["Delfi"],
        initial_state,
        0.0,
        integrator,
        termination,
    )
    return bodies, propagator_settings


def _assert_love_number_parameter_present(parameter_set):
    identifiers = parameter_set.get_parameter_identifiers()
    love_number_identifiers = [
        identifier
        for identifier in identifiers
        if identifier[0] == EstimatableParameterTypes.full_degree_tidal_love_number_type
    ]
    assert len(love_number_identifiers) == 1
    assert love_number_identifiers[0][1][0] == "Earth"


def _exercise_love_number_setup(deforming_bodies):
    bodies, propagator_settings = _create_multi_body_tide_environment()

    parameter_settings = parameters_setup.initial_states(propagator_settings, bodies)
    parameter_settings.append(
        parameters_setup.order_invariant_k_love_number("Earth", 2, deforming_bodies, False)
    )
    parameter_set = parameters_setup.create_parameter_set(parameter_settings, bodies)

    assert parameter_set.parameter_set_size > 6
    _assert_love_number_parameter_present(parameter_set)

    estimator = estimation_analysis.Estimator(
        bodies,
        parameter_set,
        [],
        propagator_settings,
        integrate_on_creation=False,
    )
    assert estimator is not None
    assert estimator.variational_solver is not None


def test_multi_body_love_number_explicit_deforming_bodies():
    _exercise_love_number_setup(["Moon", "Sun"])


def test_multi_body_love_number_empty_deforming_bodies():
    _exercise_love_number_setup([])
