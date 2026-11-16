#!/usr/bin/env python3
"""Profile native Tudat dynamics, optionally including variational equations.

No Python callbacks, plotting, observation simulation, or result-history exports
are performed in the timed integration. See README.md in this directory.
"""

import argparse
import json
import math
import os
import time


def positive(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--days", type=positive, default=120.0)
    parser.add_argument("--step", type=positive, default=30.0, help="fixed step in seconds")
    parser.add_argument("--degree", type=int, default=8, choices=range(2, 101), metavar="2..100")
    parser.add_argument("--mode", choices=("variational", "dynamics"), default="variational")
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
    settings = environment_setup.get_default_body_settings(
        ["Sun", "Earth", "Moon", "Mars"], "Earth", "J2000"
    )
    settings.add_empty_settings("Satellite")
    satellite = settings.get("Satellite")
    satellite.constant_mass = 500.0
    satellite.aerodynamic_coefficient_settings = (
        environment_setup.aerodynamic_coefficients.constant(4.0, [2.2, 0.0, 0.0])
    )
    satellite.radiation_pressure_target_settings = (
        environment_setup.radiation_pressure.cannonball_radiation_target(
            4.0, 1.3, {"Sun": ["Earth"]}
        )
    )
    bodies = environment_setup.create_system_of_bodies(settings)
    acceleration = propagation_setup.acceleration
    models = propagation_setup.create_acceleration_models(
        bodies,
        {
            "Satellite": {
                "Earth": [
                    acceleration.spherical_harmonic_gravity(args.degree, args.degree),
                    acceleration.aerodynamic(),
                ],
                "Sun": [acceleration.point_mass_gravity(), acceleration.radiation_pressure()],
                "Moon": [acceleration.point_mass_gravity()],
                "Mars": [acceleration.point_mass_gravity()],
            }
        },
        ["Satellite"],
        ["Earth"],
    )
    initial_state = element_conversion.keplerian_to_cartesian_elementwise(
        semi_major_axis=7078137.0,
        eccentricity=0.001,
        inclination=np.deg2rad(98.0),
        argument_of_periapsis=0.2,
        longitude_of_ascending_node=0.3,
        true_anomaly=0.4,
        gravitational_parameter=bodies.get_body("Earth").gravitational_parameter,
    )
    integrator = propagation_setup.integrator.runge_kutta_fixed_step_size(
        initial_time_step=args.step,
        coefficient_set=propagation_setup.integrator.CoefficientSets.rkdp_87,
    )
    propagation = propagation_setup.propagator.translational(
        ["Earth"],
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
        parameter_settings = parameters_setup.initial_states(propagation, bodies)
        parameter_settings.extend(
            [
                parameters_setup.gravitational_parameter("Earth"),
                parameters_setup.constant_drag_coefficient("Satellite"),
            ]
        )
        parameters = parameters_setup.create_parameter_set(parameter_settings, bodies)
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
    cpu_seconds = time.process_time() - cpu_start
    wall_seconds = time.perf_counter() - wall_start
    print("INTEGRATION_END", flush=True)
    if not dynamics.integration_completed_successfully:
        raise RuntimeError("Tudat integration did not complete successfully")
    # Export only after the timed region. This is a diagnostic, not an accuracy proof.
    history = dynamics.propagation_results.state_history
    final_epoch = max(history)
    final_state = np.asarray(history[final_epoch]).reshape(-1)
    if not np.all(np.isfinite(final_state)) or abs(float(final_epoch) - end) > 1e-3:
        raise RuntimeError("Non-finite final state or unexpected final epoch")
    print(
        json.dumps(
            {
                "phase": "integration",
                "wall_seconds": wall_seconds,
                "cpu_seconds": cpu_seconds,
                "saved_states": len(history),
                "final_epoch": float(final_epoch),
                "final_state": final_state.tolist(),
                "final_radius_m": float(np.linalg.norm(final_state[:3])),
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
