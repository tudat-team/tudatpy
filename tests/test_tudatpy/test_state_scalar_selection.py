"""Runtime selection tests for the double/quad state-scalar Python bindings."""

import json
import math
import subprocess
import sys
import textwrap

import numpy as np
import pytest

_CHILD_PROGRAM = r"""
import json
import math
import numpy as np
import tudatpy

requested_mode = __REQUESTED_MODE__
if requested_mode is not None:
    tudatpy.set_state_scalar_type(requested_mode)

from tudatpy import kernel
from tudatpy.kernel.dynamics import environment_setup, propagation_setup, simulator
from tudatpy.kernel.interface import spice
from tudatpy.kernel.math import interpolators

# Scalar and map conversion at the Python boundary.
scalar_interpolator = interpolators.create_one_dimensional_scalar_interpolator(
    {0.0: 1.25, 1.0: 2.75},
    interpolators.linear_interpolation(),
)
interpolated_value = scalar_interpolator.interpolate(0.25)

# Slimmed-down counterpart of the C++ point-mass high-precision test. All
# force and integration arithmetic remains in C++; Python supplies/receives
# only binary64 values.
spice.load_standard_kernels()
body_settings = environment_setup.get_default_body_settings(
    ["Earth"], "Earth", "J2000"
)
body_settings.add_empty_settings("Vehicle")
bodies = environment_setup.create_system_of_bodies(body_settings)

acceleration_models = propagation_setup.create_acceleration_models(
    bodies,
    {
        "Vehicle": {
            "Earth": [propagation_setup.acceleration.point_mass_gravity()]
        }
    },
    ["Vehicle"],
    ["Earth"],
)

radius = 8.0e6
gravitational_parameter = bodies.get("Earth").gravitational_parameter
initial_state = np.array(
    [radius, 0.0, 0.0, 0.0, math.sqrt(gravitational_parameter / radius), 0.0],
    dtype=np.float64,
)
propagator_settings = propagation_setup.propagator.translational(
    ["Earth"],
    acceleration_models,
    ["Vehicle"],
    initial_state,
    0.0,
    propagation_setup.integrator.runge_kutta_4(60.0),
    propagation_setup.propagator.time_termination(600.0),
)
dynamics_simulator = simulator.create_dynamics_simulator(
    bodies, propagator_settings
)

# This history retains Time keys and the configured C++ state scalar until
# pybind11 converts it to Time/float64 Python objects.
state_history = dynamics_simulator.state_history_time_object
final_state = next(reversed(state_history.values()))
first_epoch = next(iter(state_history))

locked_after_import = False
try:
    tudatpy.set_state_scalar_type("double")
except RuntimeError:
    locked_after_import = True

print(
    "TUDATPY_STATE_SCALAR_RESULT="
    + json.dumps(
        {
            "mode": tudatpy.get_state_scalar_type(),
            "kernel_mode": kernel._state_scalar_type,
            "bits": kernel._state_scalar_bits,
            "locked_after_import": locked_after_import,
            "interpolator_type": type(scalar_interpolator).__name__,
            "interpolated_value": interpolated_value,
            "interpolated_value_is_float": isinstance(interpolated_value, float),
            "epoch_type": type(first_epoch).__name__,
            "state_dtype": str(final_state.dtype),
            "state_shape": list(final_state.shape),
            "state_count": len(state_history),
            "final_state": final_state.tolist(),
        }
    )
)
"""

_ORBIT_CONVERGENCE_CHILD_PROGRAM = r"""
import json
import math
import numpy as np
import tudatpy

tudatpy.set_state_scalar_type(__REQUESTED_MODE__)

from tudatpy import kernel
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel.dynamics import environment_setup, propagation_setup, simulator

gravitational_parameter = 398600441800000.0
semi_major_axis = 8.0e6
propagation_time = 3600.0
step_counts = [256, 512, 1024, 2048, 4096, 8192, 16384, 32768]

body_settings = environment_setup.BodyListSettings("SSB", "J2000")
body_settings.add_empty_settings("Earth")
body_settings.add_empty_settings("Vehicle")
body_settings.get("Earth").ephemeris_settings = (
    environment_setup.ephemeris.constant(
        np.zeros(6), "SSB", "J2000"
    )
)
body_settings.get("Earth").gravity_field_settings = (
    environment_setup.gravity_field.central(gravitational_parameter)
)
bodies = environment_setup.create_system_of_bodies(body_settings)

acceleration_models = propagation_setup.create_acceleration_models(
    bodies,
    {
        "Vehicle": {
            "Earth": [propagation_setup.acceleration.point_mass_gravity()]
        }
    },
    ["Vehicle"],
    ["Earth"],
)

initial_state = element_conversion.keplerian_to_cartesian(
    np.array(
        [
            semi_major_axis,
            0.1,
            math.pi / 6.0,
            math.pi / 5.0,
            math.pi / 7.0,
            math.pi / 9.0,
        ],
        dtype=np.float64,
    ),
    gravitational_parameter,
)

final_states = []
for number_of_steps in step_counts:
    integrator_settings = propagation_setup.integrator.runge_kutta_4(
        propagation_time / number_of_steps
    )
    propagator_settings = propagation_setup.propagator.translational(
        ["Earth"],
        acceleration_models,
        ["Vehicle"],
        initial_state,
        0.0,
        integrator_settings,
        propagation_setup.propagator.time_termination(
            propagation_time, True
        ),
    )
    dynamics_simulator = simulator.create_dynamics_simulator(
        bodies, propagator_settings
    )
    final_states.append(
        next(
            reversed(
                dynamics_simulator.state_history_time_object.values()
            )
        ).tolist()
    )

print(
    "TUDATPY_ORBIT_CONVERGENCE_RESULT="
    + json.dumps(
        {
            "mode": kernel._state_scalar_type,
            "bits": kernel._state_scalar_bits,
            "semi_major_axis": semi_major_axis,
            "gravitational_parameter": gravitational_parameter,
            "step_counts": step_counts,
            "final_states": final_states,
        }
    )
)
"""

_DSN_DOPPLER_CHILD_PROGRAM = r"""
import json
import math
import numpy as np
import tudatpy

tudatpy.set_state_scalar_type(__REQUESTED_MODE__)

from tudatpy import kernel
from tudatpy.kernel.astro import element_conversion, time_representation
from tudatpy.kernel.dynamics import (
    environment,
    environment_setup,
    propagation_setup,
    simulator,
)
from tudatpy.kernel.estimation import observable_models_setup, observations_setup

solar_gravitational_parameter = 1.32712440018e20
body_settings = environment_setup.BodyListSettings("SSB", "J2000")
for body_name in ["Sun", "Earth", "Jupiter"]:
    body_settings.add_empty_settings(body_name)

body_settings.get("Sun").ephemeris_settings = (
    environment_setup.ephemeris.constant(np.zeros(6), "SSB", "J2000")
)
body_settings.get("Sun").gravity_field_settings = (
    environment_setup.gravity_field.central(solar_gravitational_parameter)
)
body_settings.get("Earth").ephemeris_settings = (
    environment_setup.ephemeris.constant(
        np.array([149597870700.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        "SSB",
        "J2000",
    )
)
body_settings.get("Earth").rotation_model_settings = (
    environment_setup.rotation_model.constant_rotation_model(
        "J2000", "IAU_Earth", np.eye(3)
    )
)
bodies = environment_setup.create_system_of_bodies(body_settings)

station_name = "QuadDsnStation"
environment_setup.add_ground_station(
    bodies.get("Earth"),
    station_name,
    np.array([6378137.0, 0.0, 0.0]),
)
frequency_interpolator = environment.PiecewiseLinearFrequencyInterpolator(
    [time_representation.Time(-10000.0)],
    [time_representation.Time(100000.0)],
    [0.25],
    [7.2e9],
)
bodies.get("Earth").get_ground_station(
    station_name
).set_transmitting_frequency_calculator(frequency_interpolator)

bodies.get("Jupiter").system_models = environment.VehicleSystems()
bodies.get(
    "Jupiter"
).system_models.set_default_transponder_turnaround_ratio_function()

initial_state = element_conversion.keplerian_to_cartesian(
    np.array(
        [
            778500000000.0,
            0.0489,
            math.radians(1.3),
            math.pi / 5.0,
            math.pi / 7.0,
            math.pi / 9.0,
        ]
    ),
    solar_gravitational_parameter,
)
acceleration_models = propagation_setup.create_acceleration_models(
    bodies,
    {
        "Jupiter": {
            "Sun": [propagation_setup.acceleration.point_mass_gravity()]
        }
    },
    ["Jupiter"],
    ["Sun"],
)
propagator_settings = propagation_setup.propagator.translational(
    ["Sun"],
    acceleration_models,
    ["Jupiter"],
    initial_state,
    0.0,
    propagation_setup.integrator.runge_kutta_4(60.0),
    propagation_setup.propagator.time_termination(86400.0, True),
)
propagator_settings.processing_settings.set_integrated_result = True
simulator.create_dynamics_simulator(bodies, propagator_settings)

links = observable_models_setup.links
link_ends = {
    links.transmitter: links.body_reference_point_link_end_id(
        "Earth", station_name
    ),
    links.retransmitter: links.body_origin_link_end_id("Jupiter"),
    links.receiver: links.body_reference_point_link_end_id(
        "Earth", station_name
    ),
}
link_definition = links.LinkDefinition(link_ends)
model_settings = observable_models_setup.model_settings
dsn_model_settings = model_settings.dsn_n_way_doppler_averaged(
    link_definition
)
observation_simulators = (
    observations_setup.observations_simulation_settings.create_observation_simulators(
        [dsn_model_settings], bodies
    )
)

ancillary = observations_setup.ancillary_settings
dsn_ancillary_settings = ancillary.dsn_n_way_doppler_ancillary_settings(
    [ancillary.FrequencyBands.x_band, ancillary.FrequencyBands.x_band],
    ancillary.FrequencyBands.x_band,
    7.2e9,
    0.1,
)
epoch_spacing = 1.0e-9
number_of_observations = 101
local_slope_interval = 1.0e-3
observation_times = [
    time_representation.Time(12, 0.125 + index * epoch_spacing)
    for index in range(number_of_observations)
]
observation_times.append(
    time_representation.Time(12, 0.125 + local_slope_interval)
)
simulation_settings = (
    observations_setup.observations_simulation_settings.tabulated_simulation_settings(
        model_settings.dsn_n_way_averaged_doppler_type,
        link_definition,
        observation_times,
        links.receiver,
        ancillary_settings=dsn_ancillary_settings,
    )
)
simulated_observations = (
    observations_setup.observations_wrapper.simulate_observations(
        [simulation_settings], observation_simulators, bodies
    )
)
observable_values = np.asarray(
    simulated_observations.concatenated_observations
).reshape(-1)
closely_spaced_observations = observable_values[:number_of_observations]
local_slope_observation = observable_values[-1]

print(
    "TUDATPY_DSN_DOPPLER_RESULT="
    + json.dumps(
        {
            "mode": kernel._state_scalar_type,
            "bits": kernel._state_scalar_bits,
            "integration_time": 0.1,
            "epoch_spacing": epoch_spacing,
            "observation_count": number_of_observations,
            "local_slope_interval": local_slope_interval,
            "local_slope_observation": local_slope_observation,
            "observations": closely_spaced_observations.tolist(),
        }
    )
)
"""

_INITIAL_STATE_ESTIMATION_CHILD_PROGRAM = r"""
import json
import math
import numpy as np
import tudatpy

tudatpy.set_state_scalar_type(__REQUESTED_MODE__)

from tudatpy import kernel
from tudatpy.kernel.astro import element_conversion, time_representation
from tudatpy.kernel.dynamics import (
    environment_setup,
    parameters_setup,
    propagation_setup,
    simulator,
)
from tudatpy.kernel.estimation import (
    estimation_analysis,
    observable_models_setup,
    observations_setup,
)

gravitational_parameter = 398600441800000.0
semi_major_axis = 8.0e6
propagation_time = 3600.0
step_count = 16384
step_size = propagation_time / step_count
observation_epochs = [
    index * propagation_time / 16 for index in range(17)
]
observation_times = [
    time_representation.Time(epoch) for epoch in observation_epochs
]

body_settings = environment_setup.BodyListSettings("SSB", "J2000")
body_settings.add_empty_settings("Earth")
body_settings.add_empty_settings("Vehicle")
body_settings.get("Earth").ephemeris_settings = (
    environment_setup.ephemeris.constant(np.zeros(6), "SSB", "J2000")
)
body_settings.get("Earth").gravity_field_settings = (
    environment_setup.gravity_field.central(gravitational_parameter)
)
bodies = environment_setup.create_system_of_bodies(body_settings)

acceleration_models = propagation_setup.create_acceleration_models(
    bodies,
    {
        "Vehicle": {
            "Earth": [propagation_setup.acceleration.point_mass_gravity()]
        }
    },
    ["Vehicle"],
    ["Earth"],
)

truth_initial_state = element_conversion.keplerian_to_cartesian(
    np.array(
        [
            semi_major_axis,
            0.1,
            math.pi / 6.0,
            math.pi / 5.0,
            math.pi / 7.0,
            math.pi / 9.0,
        ],
        dtype=np.float64,
    ),
    gravitational_parameter,
)
truth_propagator_settings = propagation_setup.propagator.translational(
    ["Earth"],
    acceleration_models,
    ["Vehicle"],
    truth_initial_state,
    0.0,
    propagation_setup.integrator.runge_kutta_4(step_size),
    propagation_setup.propagator.time_termination(
        propagation_time, True
    ),
)
truth_propagator_settings.processing_settings.set_integrated_result = True
simulator.create_dynamics_simulator(bodies, truth_propagator_settings)

links = observable_models_setup.links
link_ends = {
    links.observed_body: links.body_origin_link_end_id("Vehicle")
}
link_definition = links.LinkDefinition(link_ends)
position_observation_type = (
    observable_models_setup.model_settings.position_observable_type
)
position_model_settings = (
    observable_models_setup.model_settings.cartesian_position(link_definition)
)
observation_simulators = (
    observations_setup.observations_simulation_settings.create_observation_simulators(
        [position_model_settings], bodies
    )
)
simulation_settings = (
    observations_setup.observations_simulation_settings.tabulated_simulation_settings(
        position_observation_type,
        link_definition,
        observation_times,
        links.observed_body,
    )
)
simulated_observations = (
    observations_setup.observations_wrapper.simulate_observations(
        [simulation_settings], observation_simulators, bodies
    )
)
position_observations = np.asarray(
    simulated_observations.concatenated_observations
).reshape(-1, 3)

initial_perturbation = np.array(
    [10.0, -8.0, 6.0, 0.01, -0.008, 0.006],
    dtype=np.float64,
)
initial_guess = truth_initial_state + initial_perturbation
estimation_propagator_settings = propagation_setup.propagator.translational(
    ["Earth"],
    acceleration_models,
    ["Vehicle"],
    initial_guess,
    0.0,
    propagation_setup.integrator.runge_kutta_4(step_size),
    propagation_setup.propagator.time_termination(
        propagation_time, True
    ),
)
estimation_propagator_settings.processing_settings.set_integrated_result = True

parameter_settings = parameters_setup.initial_states(
    estimation_propagator_settings, bodies
)
parameters_to_estimate = parameters_setup.create_parameter_set(
    parameter_settings,
    bodies,
    estimation_propagator_settings,
)
estimator = estimation_analysis.Estimator(
    bodies,
    parameters_to_estimate,
    [position_model_settings],
    estimation_propagator_settings,
)
estimation_input = estimation_analysis.EstimationInput(
    simulated_observations,
    convergence_checker=estimation_analysis.estimation_convergence_checker(
        maximum_iterations=5
    ),
)
# The estimator integrated the dynamics and variational equations on
# construction. Reuse those partials during the differential corrections:
# only the initial state is estimated and its 10 m perturbation is locally
# linear, while avoiding duplicate 16384-step variational integrations.
estimation_input.define_estimation_settings(
    reintegrate_equations_on_first_iteration=False,
    reintegrate_variational_equations=False,
    print_output_to_terminal=False,
    save_state_history_per_iteration=False,
)
estimation_output = estimator.perform_estimation(estimation_input)

final_parameters = np.asarray(estimation_output.final_parameters).reshape(-1)
final_residuals = np.asarray(estimation_output.final_residuals).reshape(-1)
parameter_history = np.asarray(estimation_output.parameter_history)
residual_history = np.asarray(estimation_output.residual_history)
state_error = final_parameters - truth_initial_state

print(
    "TUDATPY_INITIAL_STATE_ESTIMATION_RESULT="
    + json.dumps(
        {
            "mode": kernel._state_scalar_type,
            "bits": kernel._state_scalar_bits,
            "propagation_time": propagation_time,
            "step_count": step_count,
            "step_size": step_size,
            "observation_count": len(observation_times),
            "position_observations": position_observations.tolist(),
            "initial_perturbation": initial_perturbation.tolist(),
            "final_parameters": final_parameters.tolist(),
            "state_error": state_error.tolist(),
            "position_error_norm": float(np.linalg.norm(state_error[:3])),
            "velocity_error_norm": float(np.linalg.norm(state_error[3:])),
            "residual_rms": float(
                np.sqrt(np.mean(np.square(final_residuals)))
            ),
            "maximum_absolute_residual": float(
                np.max(np.abs(final_residuals))
            ),
            "iteration_count": int(parameter_history.shape[1]),
            "residual_rms_history": np.sqrt(
                np.mean(np.square(residual_history), axis=0)
            ).tolist(),
            "exception_during_inversion": (
                estimation_output.exception_during_inversion
            ),
            "exception_during_propagation": (
                estimation_output.exception_during_propagation
            ),
        }
    )
)
"""


def _run_in_fresh_interpreter(requested_mode):
    program = _CHILD_PROGRAM.replace("__REQUESTED_MODE__", repr(requested_mode))
    completed = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(program)],
        check=True,
        capture_output=True,
        text=True,
    )
    result_line = next(
        line
        for line in completed.stdout.splitlines()
        if line.startswith("TUDATPY_STATE_SCALAR_RESULT=")
    )
    return json.loads(result_line.partition("=")[2])


def _run_orbit_convergence_in_fresh_interpreter(requested_mode):
    program = _ORBIT_CONVERGENCE_CHILD_PROGRAM.replace("__REQUESTED_MODE__", repr(requested_mode))
    completed = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(program)],
        check=True,
        capture_output=True,
        text=True,
    )
    result_line = next(
        line
        for line in completed.stdout.splitlines()
        if line.startswith("TUDATPY_ORBIT_CONVERGENCE_RESULT=")
    )
    return json.loads(result_line.partition("=")[2])


def _run_dsn_doppler_in_fresh_interpreter(requested_mode):
    program = _DSN_DOPPLER_CHILD_PROGRAM.replace("__REQUESTED_MODE__", repr(requested_mode))
    completed = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(program)],
        check=True,
        capture_output=True,
        text=True,
    )
    result_line = next(
        line
        for line in completed.stdout.splitlines()
        if line.startswith("TUDATPY_DSN_DOPPLER_RESULT=")
    )
    return json.loads(result_line.partition("=")[2])


def _run_initial_state_estimation_in_fresh_interpreter(requested_mode):
    program = _INITIAL_STATE_ESTIMATION_CHILD_PROGRAM.replace(
        "__REQUESTED_MODE__", repr(requested_mode)
    )
    completed = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(program)],
        check=True,
        capture_output=True,
        text=True,
    )
    result_line = next(
        line
        for line in completed.stdout.splitlines()
        if line.startswith("TUDATPY_INITIAL_STATE_ESTIMATION_RESULT=")
    )
    return json.loads(result_line.partition("=")[2])


def _normalized_step_doubling_errors(result):
    states = np.asarray(result["final_states"])
    state_differences = np.diff(states, axis=0)
    position_errors = np.linalg.norm(state_differences[:, :3], axis=1)
    velocity_errors = np.linalg.norm(state_differences[:, 3:], axis=1)
    characteristic_velocity = math.sqrt(
        result["gravitational_parameter"] / result["semi_major_axis"]
    )
    return np.maximum(
        position_errors / result["semi_major_axis"],
        velocity_errors / characteristic_velocity,
    )


def _assert_common_python_interface(result):
    assert result["mode"] == result["kernel_mode"]
    assert result["locked_after_import"]
    assert result["interpolator_type"] == "OneDimensionalInterpolatorScalar"
    assert result["interpolated_value_is_float"]
    assert result["interpolated_value"] == pytest.approx(1.625)
    assert result["epoch_type"] == "Time"
    assert result["state_dtype"] == "float64"
    assert result["state_shape"] == [6]
    assert result["state_count"] > 1
    assert all(math.isfinite(value) for value in result["final_state"])


def test_default_and_explicit_double_state_scalar_interfaces_match():
    default_result = _run_in_fresh_interpreter(None)
    explicit_result = _run_in_fresh_interpreter("double")

    _assert_common_python_interface(default_result)
    _assert_common_python_interface(explicit_result)
    assert default_result["mode"] == "double"
    assert default_result["bits"] == 53
    assert explicit_result["bits"] == 53
    assert default_result["interpolator_type"] == explicit_result["interpolator_type"]


def test_quad_state_scalar_uses_the_same_python_interface():
    import tudatpy

    if not tudatpy.quad_precision_available():
        pytest.skip("quad-precision Python bindings were disabled at configure time")

    double_result = _run_in_fresh_interpreter("double")
    quad_result = _run_in_fresh_interpreter("quad")

    _assert_common_python_interface(double_result)
    _assert_common_python_interface(quad_result)
    assert quad_result["mode"] == "quad"
    assert quad_result["bits"] == 113
    assert quad_result["interpolator_type"] == double_result["interpolator_type"]
    assert quad_result["state_dtype"] == double_result["state_dtype"]
    assert quad_result["state_shape"] == double_result["state_shape"]


def test_quad_point_mass_orbit_retains_rk4_truncation_convergence():
    import tudatpy

    if not tudatpy.quad_precision_available():
        pytest.skip("quad-precision Python bindings were disabled at configure time")

    double_result = _run_orbit_convergence_in_fresh_interpreter("double")
    quad_result = _run_orbit_convergence_in_fresh_interpreter("quad")

    double_errors = _normalized_step_doubling_errors(double_result)
    quad_errors = _normalized_step_doubling_errors(quad_result)
    double_reductions = double_errors[:-1] / double_errors[1:]
    quad_reductions = quad_errors[:-1] / quad_errors[1:]
    double_orders = np.log2(double_reductions)
    quad_orders = np.log2(quad_reductions)

    print("double step-doubling errors:", double_errors)
    print("double reductions:", double_reductions)
    print("double observed orders:", double_orders)
    print("quad step-doubling errors:", quad_errors)
    print("quad reductions:", quad_reductions)
    print("quad observed orders:", quad_orders)

    assert double_result["mode"] == "double"
    assert quad_result["mode"] == "quad"
    assert double_result["bits"] == 53
    assert quad_result["bits"] == 113
    assert double_result["step_counts"] == quad_result["step_counts"]
    step_counts = np.asarray(quad_result["step_counts"])
    assert np.all(step_counts[1:] == 2 * step_counts[:-1])
    assert np.asarray(double_result["final_states"]).shape == (len(step_counts), 6)
    assert np.asarray(quad_result["final_states"]).shape == (len(step_counts), 6)
    assert np.all(np.isfinite(double_errors))
    assert np.all(np.isfinite(quad_errors))

    # For a fourth-order method, halving the step should reduce the global
    # truncation error by approximately 2**4 = 16. Both modes initially
    # exhibit fourth-order convergence.
    assert np.all((double_orders[:3] > 3.8) & (double_orders[:3] < 4.2))
    assert np.all((quad_orders > 3.8) & (quad_orders < 4.2))
    assert np.all(np.diff(quad_errors) < 0.0)

    # At the finest halving, double precision no longer exhibits fourth-order
    # convergence, whereas the quad propagation remains truncation dominated.
    assert double_orders[-1] < 3.0
    assert quad_errors[-1] * 10.0 < double_errors[-1]


def test_quad_dsn_closed_loop_doppler_resolves_closely_spaced_observations():
    import tudatpy

    if not tudatpy.quad_precision_available():
        pytest.skip("quad-precision Python bindings were disabled at configure time")

    double_result = _run_dsn_doppler_in_fresh_interpreter("double")
    quad_result = _run_dsn_doppler_in_fresh_interpreter("quad")
    double_observations = np.asarray(double_result["observations"])
    quad_observations = np.asarray(quad_result["observations"])
    double_increments = np.diff(double_observations)
    quad_increments = np.diff(quad_observations)
    epoch_spacing = quad_result["epoch_spacing"]
    local_slope_interval = quad_result["local_slope_interval"]
    time_has_extended_long_double_precision = (
        np.finfo(np.longdouble).nmant > np.finfo(np.float64).nmant
    )

    # As in the C++ high-precision test, use a wider 1 ms interval to
    # establish the local physical trend, then check the 1 ns sequence
    # against that trend. The expected 100 ns sweep remains representable
    # after the quad result is converted to a Python float.
    quad_local_slope = (
        quad_result["local_slope_observation"] - quad_observations[0]
    ) / local_slope_interval
    offsets = np.arange(len(quad_observations)) * epoch_spacing
    expected_quad_observations = quad_observations[0] + quad_local_slope * offsets
    expected_increment = quad_local_slope * epoch_spacing
    expected_sweep = quad_local_slope * offsets[-1]
    quad_sweep = quad_observations[-1] - quad_observations[0]
    double_sweep = double_observations[-1] - double_observations[0]
    quad_increment_residuals = quad_increments - expected_increment
    quad_trend_residuals = quad_observations - expected_quad_observations
    double_trend_residuals = (
        double_observations - double_observations[0] - quad_local_slope * offsets
    )

    print(
        "quad local slope [Hz/s], expected/actual 100 ns sweep [Hz]:",
        quad_local_slope,
        expected_sweep,
        quad_sweep,
    )
    print(
        "quad max increment/trend residual [Hz]:",
        np.max(np.abs(quad_increment_residuals)),
        np.max(np.abs(quad_trend_residuals)),
    )
    print(
        "double max increment/trend residual and 100 ns sweep [Hz]:",
        np.max(np.abs(double_increments)),
        np.max(np.abs(double_trend_residuals)),
        double_sweep,
    )
    print(
        "extended long-double Time precision:",
        time_has_extended_long_double_precision,
    )

    assert double_result["mode"] == "double"
    assert quad_result["mode"] == "quad"
    assert double_result["bits"] == 53
    assert quad_result["bits"] == 113
    assert double_result["integration_time"] == pytest.approx(0.1)
    assert quad_result["integration_time"] == pytest.approx(0.1)
    assert double_result["epoch_spacing"] == pytest.approx(1.0e-9)
    assert quad_result["epoch_spacing"] == pytest.approx(1.0e-9)
    assert double_result["local_slope_interval"] == pytest.approx(1.0e-3)
    assert quad_result["local_slope_interval"] == pytest.approx(1.0e-3)
    assert len(double_observations) == double_result["observation_count"]
    assert len(quad_observations) == quad_result["observation_count"]
    assert np.all(np.isfinite(double_observations))
    assert np.all(np.isfinite(quad_observations))

    assert expected_sweep != 0.0
    if time_has_extended_long_double_precision:
        # With an extended-precision Time remainder, quad must retain the
        # small trend: its 100 ns sweep has the sign and magnitude predicted
        # by the independently sampled 1 ms local slope.
        assert quad_sweep * expected_sweep > 0.0
        assert abs(quad_sweep - expected_sweep) < 0.02 * abs(expected_sweep)
        assert np.max(np.abs(quad_increments)) < 1.0e-9
        assert np.max(np.abs(quad_increment_residuals)) < 1.0e-9
        assert np.max(np.abs(quad_trend_residuals)) < 1.0e-9
    else:
        # MSVC and Apple ARM64 use binary64 for long double, including the
        # fractional remainder in Tudat's Time class. Mirror the C++ test's
        # platform-specific limits for the resulting averaged-Doppler floor.
        assert np.max(np.abs(quad_increments)) < 2.0e-6
        assert np.max(np.abs(quad_increment_residuals)) < 2.0e-6
        assert np.max(np.abs(quad_trend_residuals)) < 1.0e-6
        assert abs(quad_sweep - expected_sweep) < 1.0e-6

    # The same model instantiated with double state scalars must not pass
    # the quad behavior accidentally: cancellation in the 0.1 s averaged
    # Doppler calculation overwhelms the nanosecond-scale physical trend.
    assert abs(double_sweep - expected_sweep) > 0.5 * abs(expected_sweep)
    assert np.max(np.abs(double_increments)) > 1.0e-3
    assert np.max(np.abs(double_trend_residuals)) > 1.0e-3


def test_quad_initial_state_estimation_from_noise_free_positions():
    import tudatpy

    if not tudatpy.quad_precision_available():
        pytest.skip("quad-precision Python bindings were disabled at configure time")

    double_result = _run_initial_state_estimation_in_fresh_interpreter("double")
    quad_result = _run_initial_state_estimation_in_fresh_interpreter("quad")
    double_observations = np.asarray(double_result["position_observations"])
    quad_observations = np.asarray(quad_result["position_observations"])
    double_state_error = np.asarray(double_result["state_error"])
    quad_state_error = np.asarray(quad_result["state_error"])
    double_residual_history = np.asarray(double_result["residual_rms_history"])
    quad_residual_history = np.asarray(quad_result["residual_rms_history"])
    expected_initial_perturbation = np.array([10.0, -8.0, 6.0, 0.01, -0.008, 0.006])
    maximum_truth_difference = np.max(np.abs(double_observations - quad_observations))

    print(
        "maximum double-vs-quad truth-position difference [m]:",
        maximum_truth_difference,
    )
    print(
        "double position/velocity error and residual RMS:",
        double_result["position_error_norm"],
        double_result["velocity_error_norm"],
        double_result["residual_rms"],
    )
    print(
        "quad position/velocity error and residual RMS:",
        quad_result["position_error_norm"],
        quad_result["velocity_error_norm"],
        quad_result["residual_rms"],
    )
    print("double residual RMS history:", double_residual_history)
    print("quad residual RMS history:", quad_residual_history)

    assert double_result["mode"] == "double"
    assert double_result["bits"] == 53
    assert quad_result["mode"] == "quad"
    assert quad_result["bits"] == 113
    assert double_result["step_count"] == 16384
    assert quad_result["step_count"] == 16384
    assert double_result["step_size"] == pytest.approx(3600.0 / 16384)
    assert quad_result["step_size"] == pytest.approx(3600.0 / 16384)
    assert double_result["observation_count"] == 17
    assert quad_result["observation_count"] == 17
    assert double_observations.shape == (17, 3)
    assert quad_observations.shape == (17, 3)
    assert double_state_error.shape == (6,)
    assert quad_state_error.shape == (6,)
    np.testing.assert_array_equal(
        double_result["initial_perturbation"],
        expected_initial_perturbation,
    )
    np.testing.assert_array_equal(
        quad_result["initial_perturbation"],
        expected_initial_perturbation,
    )
    assert np.all(np.isfinite(double_observations))
    assert np.all(np.isfinite(quad_observations))
    assert np.all(np.isfinite(double_state_error))
    assert np.all(np.isfinite(quad_state_error))
    assert not double_result["exception_during_inversion"]
    assert not double_result["exception_during_propagation"]
    assert not quad_result["exception_during_inversion"]
    assert not quad_result["exception_during_propagation"]

    # Each mode generated its own noise-free observations from its own
    # propagated trajectory. At this step size, accumulated double and quad
    # arithmetic must therefore produce measurably different numerical truth.
    assert not np.array_equal(double_observations, quad_observations)
    assert maximum_truth_difference > 1.0e-7

    # Both estimations must converge from the same 10 m / 0.01 m/s state
    # perturbation, and estimate exactly the six initial-state components.
    assert len(double_result["final_parameters"]) == 6
    assert len(quad_result["final_parameters"]) == 6
    assert double_residual_history[0] > 10.0
    assert quad_residual_history[0] > 10.0
    assert np.all(np.diff(double_residual_history) < 0.0)
    assert np.all(np.diff(quad_residual_history) < 0.0)
    assert (
        double_result["position_error_norm"]
        < np.linalg.norm(expected_initial_perturbation[:3]) * 1.0e-8
    )
    assert (
        quad_result["position_error_norm"]
        < np.linalg.norm(expected_initial_perturbation[:3]) * 1.0e-9
    )
    assert double_result["position_error_norm"] < 1.0e-7
    assert double_result["velocity_error_norm"] < 1.0e-10
    assert quad_result["position_error_norm"] < 1.0e-8
    assert quad_result["velocity_error_norm"] < 1.0e-11
    assert quad_result["residual_rms"] < 1.0e-10

    # Quad must reach a substantially lower estimation floor. These
    # inequalities also ensure that the double pipeline does not
    # accidentally exhibit the quad-precision behavior.
    assert quad_result["residual_rms"] * 100.0 < double_result["residual_rms"]
    assert quad_result["position_error_norm"] * 10.0 < double_result["position_error_norm"]
