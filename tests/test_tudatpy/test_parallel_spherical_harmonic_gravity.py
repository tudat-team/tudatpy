"""Check the spherical harmonic gravity API and compare Mars propagation across thread counts."""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import tudatpy
from tudatpy.astro import element_conversion
from tudatpy.dynamics import environment_setup, parameters_setup, propagation_setup, simulator
from tudatpy.interface import spice

acceleration = propagation_setup.acceleration


def test_backend_information():
    """Verify that is_openmp_available returns a bool, equal to TUDAT_EXPECT_OPENMP if set."""
    available = tudatpy.is_openmp_available()
    assert type(available) is bool
    expected = os.environ.get("TUDAT_EXPECT_OPENMP")
    if expected is not None:
        # The expected backend is a valid flag, and matches the backend of the loaded kernel
        assert expected in {"0", "1"}
        assert available == (expected == "1")


@pytest.mark.parametrize(
    ("spin_count", "wait_policy"),
    [(None, None), ("123", None), (None, "PASSIVE")],
)
def test_libgomp_import_default(spin_count, wait_policy):
    """
    Verify that importing tudatpy, without loading the kernel, sets GOMP_SPINCOUNT=30000 (OpenMP
    builds only) if neither GOMP_SPINCOUNT nor OMP_WAIT_POLICY is set, and leaves user-set values
    unchanged.
    """
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(tudatpy.__file__).absolute().parent.parent)
    for key, value in (("GOMP_SPINCOUNT", spin_count), ("OMP_WAIT_POLICY", wait_policy)):
        if value is None:
            env.pop(key, None)
        else:
            env[key] = value
    # The default spin count is only set for OpenMP builds without user-set wait settings
    expected_spin = spin_count
    if spin_count is None and wait_policy is None and tudatpy.is_openmp_available():
        expected_spin = "30000"
    # In a new interpreter: the tudatpy import does not load the kernel, GOMP_SPINCOUNT has the
    # expected value, and OMP_WAIT_POLICY is unchanged
    subprocess.run(
        [
            sys.executable,
            "-c",
            "import os, sys, tudatpy; "
            "assert 'tudatpy.kernel' not in sys.modules; "
            f"assert os.getenv('GOMP_SPINCOUNT') == {expected_spin!r}; "
            f"assert os.getenv('OMP_WAIT_POLICY') == {wait_policy!r}",
        ],
        env=env,
        check=True,
    )


def test_sequential_default_and_existing_positional_arguments():
    """
    Verify that spherical_harmonic_gravity uses one thread by default with the existing positional
    arguments (with and without remove_point_mass), and accepts one thread as keyword argument.
    """
    # Existing positional calls (without and with remove_point_mass) use one thread
    assert acceleration.spherical_harmonic_gravity(20, 20).number_of_threads == 1
    assert acceleration.spherical_harmonic_gravity(20, 20, True).number_of_threads == 1
    # Explicit keyword call with one thread
    assert (
        acceleration.spherical_harmonic_gravity(
            maximum_degree=20, maximum_order=20, number_of_threads=1
        ).number_of_threads
        == 1
    )


@pytest.mark.parametrize("threads", [0, -1, -20])
def test_reject_invalid_counts(threads):
    """Verify that spherical_harmonic_gravity raises a ValueError for non-positive thread counts."""
    with pytest.raises(ValueError, match="positive"):
        acceleration.spherical_harmonic_gravity(20, 20, number_of_threads=threads)


def test_settings_are_local_and_readonly():
    """
    Verify that more than one thread raises a RuntimeError without OpenMP, and that with OpenMP each
    settings object stores its own number of threads, which cannot be changed after creation.
    """
    if not tudatpy.is_openmp_available():
        # Without OpenMP, more than one thread is rejected
        with pytest.raises(RuntimeError, match="TUDAT_BUILD_WITH_OPENMP"):
            acceleration.spherical_harmonic_gravity(20, 20, number_of_threads=2)
        return
    parallel = acceleration.spherical_harmonic_gravity(20, 20, number_of_threads=4)
    sequential = acceleration.spherical_harmonic_gravity(20, 20)
    # Each settings object keeps its own number of threads
    assert parallel.number_of_threads == 4
    assert sequential.number_of_threads == 1
    # The number of threads is read-only
    with pytest.raises(AttributeError):
        parallel.number_of_threads = 2


START = 778118400.0
DURATION = 300.0


def make_case(threads, mode):
    """Mars orbiter with degree/order 20 JGMRO120D gravity, using the given number of gravity threads."""
    settings = environment_setup.get_default_body_settings(["Mars"], "Mars", "J2000")
    settings.get("Mars").gravity_field_settings = (
        environment_setup.gravity_field.predefined_spherical_harmonic(
            environment_setup.gravity_field.jgmro120d, 20
        )
    )
    bodies = environment_setup.create_system_of_bodies(settings)
    bodies.create_empty_body("Satellite")
    mars = bodies.get("Mars")
    state = element_conversion.keplerian_to_cartesian(
        np.array(
            [
                mars.gravity_field_model.reference_radius + 300000.0,
                0.001,
                np.deg2rad(93.0),
                0.2,
                0.3,
                0.4,
            ]
        ),
        mars.gravitational_parameter,
    )
    accelerations = propagation_setup.create_acceleration_models(
        bodies,
        {
            "Satellite": {
                "Mars": [acceleration.spherical_harmonic_gravity(20, 20, number_of_threads=threads)]
            }
        },
        ["Satellite"],
        ["Mars"],
    )
    propagator = propagation_setup.propagator.translational(
        ["Mars"],
        accelerations,
        ["Satellite"],
        state,
        START,
        propagation_setup.integrator.runge_kutta_fixed_step(
            30.0, propagation_setup.integrator.CoefficientSets.rkdp_87
        ),
        propagation_setup.propagator.time_termination(START + DURATION, True),
    )
    if mode == "variational":
        parameter_settings = parameters_setup.initial_states(propagator, bodies)
        parameter_settings.append(parameters_setup.gravitational_parameter("Mars"))
        parameters = parameters_setup.create_parameter_set(parameter_settings, bodies)
        solver = simulator.SingleArcVariationalSimulator(
            bodies, propagator, parameters, integrate_on_creation=False
        )
    else:
        solver = simulator.create_dynamics_simulator(
            bodies, propagator, simulate_dynamics_on_creation=False
        )
    return bodies, solver, state


def execute(case, mode):
    """Propagate the case; return the state (and for variational propagation STM and sensitivity) histories."""
    _, solver, state = case
    if mode == "variational":
        solver.integrate_full_equations(state)
        results = solver.variational_propagation_results
        histories = [
            results.dynamics_results.state_history,
            results.state_transition_matrix_history,
            results.sensitivity_matrix_history,
        ]
    else:
        solver.integrate_equations_of_motion(state)
        results = solver.propagation_results
        histories = [results.state_history]
    # The propagation completed successfully
    assert (
        results.dynamics_results if mode == "variational" else results
    ).integration_completed_successfully
    return [np.stack([np.asarray(history[t]) for t in sorted(history)]) for history in histories]


@pytest.fixture(scope="module")
def kernels():
    spice.load_standard_kernels()


@pytest.mark.parametrize("mode", ["dynamics", "variational"])
def test_full_solver_history_parity(kernels, mode):
    """
    Verify that the state history (and for variational propagation the state transition and
    sensitivity matrix histories) of a Mars orbiter are bitwise identical with 1 and 4 gravity
    threads (4 only with OpenMP).
    """
    threads = 4 if tudatpy.is_openmp_available() else 1
    reference = execute(make_case(1, mode), mode)
    actual = execute(make_case(threads, mode), mode)
    # Each history is bitwise identical for both numbers of threads
    for reference_history, actual_history in zip(reference, actual):
        np.testing.assert_array_equal(actual_history, reference_history)


def test_python_rotation_callback_stays_on_caller_and_remains_alive(kernels):
    """
    Verify that a Python rotation callback of the body exerting a multi-thread gravity acceleration
    is kept alive by the rotation model after all Python references to it are deleted, and is only
    called on the calling thread.
    """
    import gc
    import threading
    import weakref

    threads = 4 if tudatpy.is_openmp_available() else 1
    case = make_case(threads, "dynamics")
    caller = threading.get_ident()
    calls = []

    class RotationCallback:
        def __call__(self, epoch):
            calls.append(threading.get_ident())
            return np.eye(3)

    callback = RotationCallback()
    reference = weakref.ref(callback)
    rotation = environment_setup.rotation_model.custom_rotation_model(
        "J2000", "IAU_MARS", callback, 1.0
    )
    environment_setup.add_rotation_model(case[0], "Mars", rotation)
    del callback, rotation
    gc.collect()
    # The rotation model keeps the callback alive
    assert reference() is not None
    execute(case, "dynamics")
    # The callback was called during propagation, and only on the calling thread (never by workers)
    assert calls
    assert set(calls) == {caller}
