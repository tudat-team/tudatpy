#!/usr/bin/env python3
"""Native sequential multi-arc Mars variational benchmark (no Python arc scheduler)."""

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
    parser.add_argument("--arcs", type=int, default=8)
    parser.add_argument("--days-per-arc", type=positive, default=2.0)
    parser.add_argument("--degree", type=int, choices=range(2, 121), default=50, metavar="2..120")
    parser.add_argument("--step", type=positive, default=30.0)
    args = parser.parse_args()
    if args.arcs < 1 or args.days_per_arc * 86400 < 8 * args.step:
        parser.error("need at least one arc and eight steps per arc for interpolation")
    setup_start = time.perf_counter()
    import numpy as np
    import tudatpy
    from tudatpy.interface import spice
    from tudatpy.dynamics import environment_setup as es, propagation_setup as ps
    from tudatpy.dynamics import parameters_setup, simulator
    from tudatpy.astro import element_conversion
    from tudatpy.astro.time_representation import DateTime

    spice.load_standard_kernels()
    settings = es.get_default_body_settings(["Mars"], "Mars", "J2000")
    settings.get("Mars").gravity_field_settings = es.gravity_field.predefined_spherical_harmonic(
        es.gravity_field.jgmro120d, maximum_degree=args.degree
    )
    settings.add_empty_settings("Satellite")
    bodies = es.create_system_of_bodies(settings)
    gravity = bodies.get_body("Mars").gravity_field_model
    if min(gravity.cosine_coefficients.shape) < args.degree + 1:
        raise RuntimeError("Insufficient gravity coefficients")
    start = DateTime(2024, 8, 28).to_epoch()
    duration = args.days_per_arc * 86400.0
    starts, ends, arc_settings = [], [], []
    for arc in range(args.arcs):
        # A gap avoids ambiguity at arc boundaries. Every initial state is
        # explicitly prescribed; no state is obtained from the preceding arc.
        epoch = start + arc * (duration + 3600.0)
        state = element_conversion.keplerian_to_cartesian_elementwise(
            semi_major_axis=gravity.reference_radius + 300000.0,
            eccentricity=0.001,
            inclination=np.deg2rad(93.0),
            argument_of_periapsis=0.2,
            longitude_of_ascending_node=0.3,
            true_anomaly=(0.4 + 0.7 * arc) % (2 * np.pi),
            gravitational_parameter=gravity.gravitational_parameter,
        )
        # Each arc gets its own acceleration model and settings. The existing
        # native multi-arc solver still shares SystemOfBodies between arcs.
        models = ps.create_acceleration_models(
            bodies,
            {
                "Satellite": {
                    "Mars": [ps.acceleration.spherical_harmonic_gravity(args.degree, args.degree)]
                }
            },
            ["Satellite"],
            ["Mars"],
        )
        integrator = ps.integrator.runge_kutta_fixed_step_size(
            initial_time_step=args.step, coefficient_set=ps.integrator.CoefficientSets.rkdp_87
        )
        arc_settings.append(
            ps.propagator.translational(
                ["Mars"],
                models,
                ["Satellite"],
                state,
                epoch,
                integrator,
                ps.propagator.time_termination(
                    epoch + duration, terminate_exactly_on_final_condition=True
                ),
            )
        )
        starts.append(epoch)
        ends.append(epoch + duration)
    propagation = ps.propagator.multi_arc(arc_settings, transfer_state_to_next_arc=False)
    parameter_settings = parameters_setup.initial_states(propagation, bodies, starts)
    parameter_settings.append(parameters_setup.gravitational_parameter("Mars"))
    parameters = parameters_setup.create_parameter_set(parameter_settings, bodies)
    print(
        json.dumps(
            {
                "phase": "setup",
                "seconds": time.perf_counter() - setup_start,
                "pid": os.getpid(),
                "tudatpy": tudatpy.__file__,
                **vars(args),
                "parameter_count": 6 * args.arcs + 1,
                "shared_parameter": "Mars gravitational parameter",
                "arc_starts": starts,
                "arc_ends": ends,
                "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
                "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
            }
        ),
        flush=True,
    )
    # Some installed versions expose only the base VariationalSimulator for
    # multi-arc. The factory works there too. Time construction + integration
    # together, consistently for all arc counts and library versions.
    print("NATIVE_SOLVER_BEGIN", flush=True)
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    solver = simulator.create_variational_equations_solver(
        bodies, propagation, parameters, simulate_dynamics_on_creation=True
    )
    wall_seconds, cpu_seconds = time.perf_counter() - wall_start, time.process_time() - cpu_start
    print("NATIVE_SOLVER_END", flush=True)
    interface = solver.state_transition_interface
    summaries = []
    for arc, epoch in enumerate(ends):
        matrix = np.asarray(interface.full_state_transition_sensitivity_at_epoch(epoch))
        if matrix.shape != (6, 6 * args.arcs + 1) or not np.all(np.isfinite(matrix)):
            raise RuntimeError(f"Invalid combined matrix for arc {arc}: {matrix.shape}")
        local = matrix[:, 6 * arc : 6 * (arc + 1)]
        inactive = np.delete(matrix[:, :-1], np.s_[6 * arc : 6 * (arc + 1)], axis=1)
        if np.any(inactive != 0) or not np.any(local) or not np.any(matrix[:, -1]):
            raise RuntimeError(f"Inconsistent local/shared parameter blocks for arc {arc}")
        summaries.append(
            {
                "arc": arc,
                "epoch": epoch,
                "matrix_shape": matrix.shape,
                "stm_norm": float(np.linalg.norm(local)),
                "mu_sensitivity": matrix[:, -1].tolist(),
                "stm": local.tolist(),
            }
        )
    print(
        json.dumps(
            {
                "phase": "native_solver",
                "wall_seconds": wall_seconds,
                "cpu_seconds": cpu_seconds,
                "arc_checks": summaries,
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
