#!/usr/bin/env python3
"""Validate and time model-local gravity threads through the ordinary Python API.

Requires a consistent build of the current native libraries and Python kernel,
plus the standard kernels and JGMRO120D resource. No replacement objects are used.
"""

import argparse
import hashlib
import os
import json
from pathlib import Path
import platform
import statistics
import sys
import time

import numpy as np
from tudatpy import kernel
from tudatpy.astro import element_conversion
from tudatpy.dynamics import environment_setup, parameters_setup, propagation_setup, simulator
from tudatpy.interface import spice

START = 778118400.0


def make_case(threads, degree, duration, mode):
    """Create isolated environment, acceleration, solver and result storage."""
    settings = environment_setup.get_default_body_settings(["Mars"], "Mars", "J2000")
    settings.get("Mars").gravity_field_settings = (
        environment_setup.gravity_field.predefined_spherical_harmonic(
            environment_setup.gravity_field.jgmro120d, degree
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
    acceleration = propagation_setup.create_acceleration_models(
        bodies,
        {
            "Satellite": {
                "Mars": [
                    propagation_setup.acceleration.spherical_harmonic_gravity(
                        degree, degree, number_of_threads=threads
                    )
                ]
            }
        },
        ["Satellite"],
        ["Mars"],
    )
    propagator = propagation_setup.propagator.translational(
        ["Mars"],
        acceleration,
        ["Satellite"],
        state,
        START,
        propagation_setup.integrator.runge_kutta_fixed_step(
            30.0, propagation_setup.integrator.CoefficientSets.rkdp_87
        ),
        propagation_setup.propagator.time_termination(START + duration, True),
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
    # Explicit owners remain alive until after propagation and collection.
    return bodies, propagator, solver, state


def execute(case, mode):
    _, _, solver, state = case
    if mode == "variational":
        solver.integrate_full_equations(state)
        results = solver.variational_propagation_results
    else:
        solver.integrate_equations_of_motion(state)
        results = solver.propagation_results
    return results


def collect(results, mode):
    dynamics = results.dynamics_results
    if not dynamics.integration_completed_successfully:
        raise RuntimeError("Native integration did not complete successfully")
    histories = {"state": dynamics.state_history}
    if mode == "variational":
        histories.update(
            stm=results.state_transition_matrix_history,
            sensitivity=results.sensitivity_matrix_history,
        )
    arrays = {}
    for name, history in histories.items():
        epochs = np.array(sorted(history), dtype=float)
        values = np.stack([np.asarray(history[t]) for t in epochs])
        if not np.isfinite(epochs).all() or not np.isfinite(values).all():
            raise AssertionError(f"Non-finite {name} history")
        arrays[name] = (epochs, values)
    return arrays


def compare(reference, actual, duration):
    """Preserve the trial's tolerance; check every state/STM/sensitivity entry."""
    expected_shapes = {"state": (6,), "stm": (6, 6), "sensitivity": (6, 1)}
    maximum = 0.0
    if reference.keys() != actual.keys():
        raise AssertionError("Result interfaces differ")
    for name, (epochs, values) in actual.items():
        ref_epochs, ref_values = reference[name]
        np.testing.assert_array_equal(epochs, ref_epochs)
        if epochs[0] != START or epochs[-1] != START + duration:
            raise AssertionError("Unexpected propagation interval")
        if len(epochs) != round(duration / 30.0) + 1:
            raise AssertionError("Unexpected number of saved epochs")
        if values.shape != ref_values.shape or values.shape[1:] != expected_shapes[name]:
            raise AssertionError(f"Unexpected {name} shape: {values.shape}")
        difference = np.abs(values - ref_values)
        if not np.all(difference <= 1e-14 + 1e-11 * np.abs(ref_values)):
            raise AssertionError(f"{name} differs beyond the existing trial tolerance")
        maximum = max(maximum, float(difference.max()))
    return maximum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--degree", type=int, default=100)
    parser.add_argument("--duration", type=float, default=86400.0)
    parser.add_argument("--mode", choices=["dynamics", "variational"], default="dynamics")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.repeats < 1 or args.degree < 2 or args.degree > 120:
        parser.error("Use positive repeats and degree 2..120 for JGMRO120D")
    if args.duration <= 0 or args.duration % 30 != 0:
        parser.error("Duration must be a positive multiple of 30 seconds")
    if len(set(args.workers)) != len(args.workers) or any(t < 1 for t in args.workers):
        parser.error("Workers must be distinct positive integers")
    if max(args.workers) > 1 and not propagation_setup.acceleration.is_openmp_available():
        parser.error("Multiple threads require a Tudat build with OpenMP enabled")
    if args.output.exists() or args.output.with_suffix(".metadata.json").exists():
        parser.error("Choose a new output path; existing results are not overwritten")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    metadata = {
        "python": sys.version,
        "platform": platform.platform(),
        "kernel": kernel.__file__,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "thread_environment": {
            k: os.environ.get(k)
            for k in (
                "OPENBLAS_NUM_THREADS",
                "OMP_NUM_THREADS",
                "OMP_DYNAMIC",
                "OMP_THREAD_LIMIT",
                "OMP_WAIT_POLICY",
            )
        },
        "openmp": propagation_setup.acceleration.is_openmp_available(),
        "arguments": {**vars(args), "output": str(args.output)},
    }
    args.output.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    spice.load_standard_kernels()
    print("Building untimed sequential reference...", flush=True)
    case = make_case(1, args.degree, args.duration, args.mode)
    reference = collect(execute(case, args.mode), args.mode)
    compare(reference, reference, args.duration)
    del case
    rows = []
    with args.output.open("x") as output:
        for repetition in range(args.repeats):
            order = (
                args.workers[repetition % len(args.workers) :]
                + args.workers[: repetition % len(args.workers)]
            )
            for threads in order:
                print(f"Repetition {repetition+1}/{args.repeats}, threads={threads}", flush=True)
                before = time.perf_counter()
                case = make_case(threads, args.degree, args.duration, args.mode)
                ready = time.perf_counter()
                results = execute(case, args.mode)
                done = time.perf_counter()
                actual = collect(results, args.mode)
                difference = compare(reference, actual, args.duration)
                row = dict(
                    threads=threads,
                    repetition=repetition + 1,
                    mode=args.mode,
                    degree=args.degree,
                    duration=args.duration,
                    setup_seconds=ready - before,
                    integration_seconds=done - ready,
                    total_seconds=done - before,
                    max_abs_difference=difference,
                    validation="passed",
                )
                output.write(json.dumps(row) + "\n")
                output.flush()
                rows.append(row)
                del case, results, actual
    baseline = (
        statistics.median(r["integration_seconds"] for r in rows if r["threads"] == 1)
        if 1 in args.workers
        else None
    )
    for threads in args.workers:
        samples = [r["integration_seconds"] for r in rows if r["threads"] == threads]
        totals = [r["total_seconds"] for r in rows if r["threads"] == threads]
        median = statistics.median(samples)
        speedup = f", speedup {baseline/median:.2f}x" if baseline else ""
        print(
            f"threads={threads}: integration median {median:.6f}s (range {min(samples):.6f}..{max(samples):.6f}){speedup}; total median {statistics.median(totals):.6f}s"
        )


if __name__ == "__main__":
    main()
