"""Binding, ownership and numerical parity checks on short native Mars arcs."""

import gc

import numpy as np
import pytest

from tudatpy.dynamics import propagation, simulator
from tudatpy.data_input.environment_data import spice

EPOCH = 778_118_400.0


@pytest.fixture(scope="module")
def model():
    spice.load_standard_kernels()
    return simulator.MarsGravityArcModel(degree=8)


def states_for(model):
    radius = model.reference_radius + 300_000.0
    speed = np.sqrt(model.gravitational_parameter / radius)
    states = np.array([[radius, 0, 0, 0, 0, speed]] * 2)
    states[1, 1] = 100.0
    return states


def make_simulation(model, **overrides):
    arguments = dict(
        model=model,
        arc_start_times=[EPOCH, EPOCH + 1800],
        initial_states=states_for(model),
        arc_duration=600.0,
        step_size=30.0,
    )
    arguments.update(overrides)
    return simulator.IndependentArcVariationalSimulator(**arguments)


def histories(simulation):
    result = []
    for dynamics, variations in zip(
        simulation.propagation_results.single_arc_results,
        simulation.variational_propagation_results.single_arc_results,
    ):
        for history in (
            dynamics.state_history,
            variations.state_transition_matrix_history,
            variations.sensitivity_matrix_history,
        ):
            times, values = np.array(list(history)), np.array(list(history.values()))
            assert len(times) >= 8
            assert np.isfinite(values).all()
            result.append((times, values))
    return result


def assert_histories_equal(reference, actual):
    # Same equations and per-arc summation order: use the native acceptance
    # tolerance, not an orbit-accuracy claim or a tolerance fitted to this run.
    assert len(reference) == len(actual) == 6
    for (times, values), (other_times, other_values) in zip(reference, actual):
        np.testing.assert_array_equal(times, other_times)
        np.testing.assert_allclose(values, other_values, rtol=1e-11, atol=1e-14)


def test_sequential_parallel_results_and_shared_parameter(model):
    simulation = make_simulation(model)
    assert simulation.number_of_workers == 1
    with pytest.raises(RuntimeError, match="no current"):
        _ = simulation.propagation_results
    simulation.integrate_full_equations()
    assert not simulation.last_integration_used_independent_arcs
    reference = histories(simulation)
    results = simulation.propagation_results
    variations = simulation.variational_propagation_results
    interface = simulation.state_transition_interface
    assert isinstance(results, propagation.MultiArcSimulationResults)
    assert isinstance(variations, propagation.MultiArcVariationalSimulationResults)
    simulation.set_independent_arc_execution(2)
    simulation.integrate_full_equations()
    assert simulation.last_integration_used_independent_arcs
    assert simulation.propagation_results is results
    assert simulation.variational_propagation_results is variations
    assert simulation.state_transition_interface is interface
    assert_histories_equal(reference, histories(simulation))
    matrix = interface.full_state_transition_sensitivity_at_epoch(EPOCH + 300)
    assert matrix.shape == (6, 13)
    np.testing.assert_array_equal(matrix[:, 6:12], 0.0)
    assert np.linalg.norm(matrix[:, -1]) > 0
    parameters = simulation.parameter_vector
    parameters[0] += 5.0
    parameters[-1] *= 1.000001
    simulation.reset_parameter_estimate(parameters)
    parallel_updated = histories(simulation)
    np.testing.assert_allclose(simulation.parameter_vector, parameters)
    simulation.clear_independent_arc_execution()
    simulation.integrate_full_equations()
    assert not simulation.last_integration_used_independent_arcs
    assert_histories_equal(parallel_updated, histories(simulation))


def test_copied_inputs_and_retained_native_results(model):
    states = states_for(model)
    expected = states.ravel().copy()
    epochs = np.array([EPOCH, EPOCH + 1800])
    local_model = simulator.MarsGravityArcModel(degree=8)
    simulation = make_simulation(
        local_model, initial_states=states, arc_start_times=epochs, number_of_workers=2
    )
    states[:] = np.nan
    epochs[:] = np.nan
    del local_model, states, epochs
    gc.collect()
    np.testing.assert_array_equal(simulation.initial_states, expected)
    simulation.integrate_full_equations()
    results = simulation.propagation_results
    variations = simulation.variational_propagation_results
    interface = simulation.state_transition_interface
    del simulation
    gc.collect()
    assert len(results.single_arc_results[0].state_history) >= 20
    assert len(variations.single_arc_results[0].state_transition_matrix_history) >= 20
    assert np.isfinite(interface.full_state_transition_sensitivity_at_epoch(EPOCH + 300)).all()


def test_reset_override_and_result_modes(model):
    simulation = make_simulation(model, number_of_workers=2)
    override = simulation.initial_states
    override[1] += 20.0
    simulation.integrate_full_equations(override)
    np.testing.assert_array_equal(simulation.parameter_vector[:12], override)
    np.testing.assert_array_equal(
        simulation.propagation_results.single_arc_results[0].state_history[EPOCH], override[:6]
    )
    retained_results = simulation.propagation_results
    previous_initial = retained_results.single_arc_results[0].state_history[EPOCH].copy()
    parameters = simulation.parameter_vector
    parameters[0] += 5.0
    simulation.reset_parameter_estimate(parameters, reintegrate=False)
    # A deferred reset must not silently perform a dynamics-only integration.
    np.testing.assert_array_equal(
        retained_results.single_arc_results[0].state_history[EPOCH], previous_initial
    )
    with pytest.raises(RuntimeError, match="no current"):
        _ = simulation.variational_propagation_results
    simulation.integrate_equations_of_motion_only()
    assert simulation.propagation_results.single_arc_results
    with pytest.raises(RuntimeError, match="no current"):
        _ = simulation.state_transition_interface
    simulation.integrate_full_equations()
    assert simulation.variational_propagation_results.single_arc_results


def test_coupled_arcs_remain_sequential_on_reintegration(model):
    sequential = make_simulation(
        model, arc_start_times=[EPOCH, EPOCH + 300], transfer_state_to_next_arc=True
    )
    sequential.integrate_full_equations()
    simulation = make_simulation(
        model,
        arc_start_times=[EPOCH, EPOCH + 300],
        transfer_state_to_next_arc=True,
        number_of_workers=4,
    )
    simulation.integrate_full_equations()
    assert not simulation.last_integration_used_independent_arcs
    assert_histories_equal(histories(sequential), histories(simulation))
    simulation.integrate_full_equations()
    assert not simulation.last_integration_used_independent_arcs


@pytest.mark.parametrize(
    "overrides",
    [
        {"number_of_workers": 0},
        {"number_of_workers": -1},
        {"initial_states": np.zeros((2, 5))},
        {"initial_states": np.full((2, 6), np.nan)},
        {"arc_start_times": [EPOCH, EPOCH]},
        {"arc_start_times": [EPOCH, np.inf]},
        {"arc_start_times": [1e30, 2e30]},
        {"arc_duration": 0.0},
        {"step_size": 0.0},
        {"step_size": 100.0},
        {"arc_duration": 1e30},
        {"arc_duration": 1e-12, "step_size": 1e-14},
    ],
)
def test_invalid_inputs(model, overrides):
    with pytest.raises(ValueError):
        make_simulation(model, **overrides)


@pytest.mark.parametrize("degree,order", [(0, -1), (121, -1), (8, 9), (8, -2)])
def test_invalid_model(degree, order):
    with pytest.raises(ValueError):
        simulator.MarsGravityArcModel(degree, order)


def test_no_python_callback_escape_hatches(model):
    simulation = make_simulation(model)
    with pytest.raises(TypeError):
        simulation.set_independent_arc_execution(2, lambda *args: None)
    with pytest.raises(TypeError):
        simulator.IndependentArcVariationalSimulator(
            lambda: model, [EPOCH, EPOCH + 1800], states_for(model), 600.0
        )
    assert not hasattr(simulation, "bodies")
    assert not hasattr(simulation, "dynamics_simulator")
    parameters = simulation.parameter_vector
    parameters[-1] = -1.0
    with pytest.raises(ValueError):
        simulation.reset_parameter_estimate(parameters)


def test_worker_failure_is_translated_and_recoverable(model):
    simulation = make_simulation(model, number_of_workers=2)
    original = simulation.parameter_vector
    impossible = original.copy()
    impossible[-1] = 1e308
    with pytest.raises(RuntimeError, match="arc"):
        simulation.reset_parameter_estimate(impossible)
    with pytest.raises(RuntimeError, match="no current"):
        _ = simulation.propagation_results
    simulation.reset_parameter_estimate(original)
    assert simulation.last_integration_used_independent_arcs
