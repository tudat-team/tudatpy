#!/usr/bin/env python3
"""Benchmark Mars spherical-harmonic gravity, with optional variational equations.

Run each degree/mode in a fresh process. See MARS_GRAVITY.md for profiling.
"""

import argparse
import json
import math
import os
import time


def positive(value):
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--degree",
        type=int,
        choices=range(2, 121),
        default=100,
        metavar="2..120",
        help="maximum degree AND order",
    )
    parser.add_argument("--days", type=positive, default=30.0)
    parser.add_argument("--step", type=positive, default=30.0, help="fixed RKDP87 step, seconds")
    parser.add_argument(
        "--altitude-km",
        type=positive,
        default=300.0,
        help="semi-major axis minus gravity reference radius, km",
    )
    parser.add_argument("--mode", choices=("dynamics", "variational"), default="dynamics")
    args = parser.parse_args()
    if args.days * 86400 < args.step:
        parser.error("duration must be at least one integration step")

    setup_start = time.perf_counter()
    import numpy as np
    import tudatpy
    from tudatpy.interface import spice
    from tudatpy.dynamics import environment_setup, propagation_setup, parameters_setup, simulator
    from tudatpy.astro import element_conversion
    from tudatpy.astro.time_representation import DateTime

    spice.load_standard_kernels()
    start = DateTime(2024, 8, 28).to_epoch()
    end = start + args.days * 86400.0
    settings = environment_setup.get_default_body_settings(["Mars"], "Mars", "J2000")
    # Explicit model selection, independent of any future default changes.
    settings.get("Mars").gravity_field_settings = (
        environment_setup.gravity_field.predefined_spherical_harmonic(
            environment_setup.gravity_field.jgmro120d, maximum_degree=args.degree
        )
    )
    settings.add_empty_settings("Satellite")
    bodies = environment_setup.create_system_of_bodies(settings)
    gravity = bodies.get_body("Mars").gravity_field_model
    coefficient_shape = gravity.cosine_coefficients.shape
    if min(coefficient_shape) < args.degree + 1:
        raise RuntimeError(f"Insufficient gravity coefficients: {coefficient_shape}")
    radius = gravity.reference_radius
    semi_major_axis = radius + args.altitude_km * 1000.0
    eccentricity = 0.001
    if semi_major_axis * (1.0 - eccentricity) <= radius:
        parser.error("initial periapsis must be above the gravity reference radius")
    initial_state = element_conversion.keplerian_to_cartesian_elementwise(
        semi_major_axis=semi_major_axis,
        eccentricity=eccentricity,
        inclination=np.deg2rad(93.0),
        argument_of_periapsis=0.2,
        longitude_of_ascending_node=0.3,
        true_anomaly=0.4,
        gravitational_parameter=gravity.gravitational_parameter,
    )
    models = propagation_setup.create_acceleration_models(
        bodies,
        {
            "Satellite": {
                "Mars": [
                    propagation_setup.acceleration.spherical_harmonic_gravity(
                        args.degree, args.degree
                    )
                ]
            }
        },
        ["Satellite"],
        ["Mars"],
    )
    integrator = propagation_setup.integrator.runge_kutta_fixed_step_size(
        initial_time_step=args.step,
        coefficient_set=propagation_setup.integrator.CoefficientSets.rkdp_87,
    )
    propagation = propagation_setup.propagator.translational(
        ["Mars"],
        models,
        ["Satellite"],
        initial_state,
        start,
        integrator,
        propagation_setup.propagator.time_termination(
            end, terminate_exactly_on_final_condition=True
        ),
    )
    if args.mode == "variational":
        # Six initial-state parameters: propagate the 6x6 STM. No individual
        # gravity coefficients are estimated, to keep the comparison focused.
        parameters = parameters_setup.create_parameter_set(
            parameters_setup.initial_states(propagation, bodies), bodies
        )
        solver = simulator.create_variational_equations_solver(
            bodies, propagation, parameters, simulate_dynamics_on_creation=False
        )
        dynamics = solver.dynamics_simulator
    else:
        dynamics = simulator.create_dynamics_simulator(
            bodies, propagation, simulate_dynamics_on_creation=False
        )

    print(
        json.dumps(
            {
                "phase": "setup",
                "seconds": time.perf_counter() - setup_start,
                "pid": os.getpid(),
                "tudatpy": tudatpy.__file__,
                **vars(args),
                "gravity_model": "jgmro120d",
                "coefficient_shape": coefficient_shape,
                "gravity_reference_radius_m": radius,
                "gravitational_parameter_m3_s2": gravity.gravitational_parameter,
                "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
                "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
            }
        ),
        flush=True,
    )
    print("INTEGRATION_BEGIN", flush=True)
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    if args.mode == "variational":
        solver.integrate_full_equations(initial_state)
    else:
        dynamics.integrate_equations_of_motion(initial_state)
    wall_seconds = time.perf_counter() - wall_start
    cpu_seconds = time.process_time() - cpu_start
    print("INTEGRATION_END", flush=True)
    if not dynamics.integration_completed_successfully:
        raise RuntimeError("Tudat integration did not complete successfully")
    # Conversion to Python is deliberately outside the integration timer.
    history = dynamics.propagation_results.state_history
    final_epoch = max(history)
    final_state = np.asarray(history[final_epoch]).reshape(-1)
    radii = np.array(
        [np.linalg.norm(np.asarray(state).reshape(-1)[:3]) for state in history.values()]
    )
    if (
        not np.all(np.isfinite(final_state))
        or not np.all(np.isfinite(radii))
        or abs(float(final_epoch) - end) > 1e-3
        or np.min(radii) <= radius
    ):
        raise RuntimeError(
            "Invalid state, unexpected final epoch, or sampled radius below reference surface"
        )
    result = {
        "phase": "integration",
        "wall_seconds": wall_seconds,
        "cpu_seconds": cpu_seconds,
        "saved_states": len(history),
        "final_epoch": float(final_epoch),
        "final_state": final_state.tolist(),
        "minimum_sampled_altitude_m": float(np.min(radii) - radius),
    }
    if args.mode == "variational":
        matrix = np.asarray(
            solver.state_transition_interface.full_state_transition_sensitivity_at_epoch(
                final_epoch
            )
        )
        if matrix.shape != (6, 6) or not np.all(np.isfinite(matrix)):
            raise RuntimeError("Invalid final state transition matrix")
        result["final_stm_shape"] = matrix.shape
        result["final_stm_norm"] = float(np.linalg.norm(matrix))
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
