"""Tests of omitting zero sensitivity matrix columns (EstimatableParameterSet.reduce_sensitivity_matrix_columns)."""

import numpy as np

from tudatpy.dynamics import environment_setup, parameters_setup, propagation_setup, simulator
from tudatpy.interface import spice


def _propagate_variational_equations(reduce_columns):
    spice.load_standard_kernels()
    bodies = environment_setup.create_system_of_bodies(
        environment_setup.get_default_body_settings(["Earth"], "Earth", "J2000")
    )
    bodies.create_empty_body("Delfi")

    acceleration_settings = {
        "Delfi": {
            "Earth": [
                propagation_setup.acceleration.point_mass_gravity(),
                propagation_setup.acceleration.empirical(),
            ]
        }
    }
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, acceleration_settings, ["Delfi"], ["Earth"]
    )
    propagator_settings = propagation_setup.propagator.translational(
        ["Earth"],
        acceleration_models,
        ["Delfi"],
        np.array([7.0e6, 0.0, 0.0, 0.0, 7.5e3, 0.0]),
        0.0,
        propagation_setup.integrator.runge_kutta_fixed_step(
            10.0, propagation_setup.integrator.CoefficientSets.rk_4
        ),
        propagation_setup.propagator.time_termination(4000.0),
    )

    # First empirical acceleration arc ends at, and last arc starts after, the propagation interval [0, 4000]
    parameter_settings = parameters_setup.initial_states(propagator_settings, bodies)
    parameter_settings.append(parameters_setup.gravitational_parameter("Earth"))
    parameter_settings.append(
        parameters_setup.arcwise_constant_empirical_acceleration_terms(
            "Delfi", "Earth", [-1000.0, 0.0, 1000.0, 2000.0, 5000.0]
        )
    )
    parameter_set = parameters_setup.create_parameter_set(
        parameter_settings, bodies, propagator_settings
    )
    assert not parameter_set.reduce_sensitivity_matrix_columns
    parameter_set.reduce_sensitivity_matrix_columns = reduce_columns
    assert parameter_set.reduce_sensitivity_matrix_columns == reduce_columns

    return simulator.create_variational_equations_solver(bodies, propagator_settings, parameter_set)


def _assert_columns_close(matrix, reference_matrix):
    assert matrix.shape == reference_matrix.shape
    column_errors = np.linalg.norm(matrix - reference_matrix, axis=0)
    assert np.all(column_errors <= 1.0e-12 * np.linalg.norm(reference_matrix, axis=0))


def test_reduced_sensitivity_matrix_columns():
    full_solver = _propagate_variational_equations(False)
    reduced_solver = _propagate_variational_equations(True)

    # Sensitivity matrix history is returned with all columns, with the columns of the inactive arcs equal to zero
    full_history = full_solver.sensitivity_matrix_history
    reduced_history = reduced_solver.sensitivity_matrix_history
    assert sorted(full_history.keys()) == sorted(reduced_history.keys())
    for epoch in full_history.keys():
        assert reduced_history[epoch].shape == (6, 16)
        _assert_columns_close(reduced_history[epoch], full_history[epoch])
        assert np.all(reduced_history[epoch][:, 1:4] == 0.0)
        assert np.all(reduced_history[epoch][:, 13:16] == 0.0)

    # Interpolated state transition and sensitivity matrices are returned with all columns
    for epoch in [500.0, 1500.0, 3500.0]:
        _assert_columns_close(
            reduced_solver.state_transition_interface.full_state_transition_sensitivity_at_epoch(
                epoch
            ),
            full_solver.state_transition_interface.full_state_transition_sensitivity_at_epoch(
                epoch
            ),
        )
