import numpy as np

from tudatpy.astro import element_conversion, two_body_dynamics
from tudatpy.dynamics import environment_setup, propagation_setup, simulator
from tudatpy.interface import spice

SPACECRAFT_NAME = "Spacecraft"
CENTRAL_BODY_NAME = "Earth"
POSITION_TOLERANCE_METERS = 1.0e-4
VELOCITY_TOLERANCE_METERS_PER_SECOND = 1.0e-7
MINIMUM_PERTURBED_ERROR_FACTOR = 1.0e4


class PerturbationRejectingAcceleration:
    def __init__(self, spacecraft, central_body):
        self.spacecraft = spacecraft
        self.central_body = central_body
        self.central_body_gravitational_parameter = central_body.gravitational_parameter
        self.acceleration_models = None

    def set_acceleration_models(self, acceleration_models):
        self.acceleration_models = list(acceleration_models)

    def acceleration(self, time):
        if np.isnan(time):
            return np.zeros(3)
        if self.acceleration_models is None:
            raise RuntimeError("Acceleration models have not been assigned.")

        relative_position = self.spacecraft.position - self.central_body.position
        point_mass_acceleration = (
            -self.central_body_gravitational_parameter
            * relative_position
            / np.linalg.norm(relative_position) ** 3
        )

        configured_acceleration = np.zeros(3)
        for acceleration_model in self.acceleration_models:
            configured_acceleration += acceleration_model.update_and_get_acceleration(time)

        return point_mass_acceleration - configured_acceleration


def create_propagator_settings(acceleration_models, initial_state, final_time):
    return propagation_setup.propagator.translational(
        [CENTRAL_BODY_NAME],
        acceleration_models,
        [SPACECRAFT_NAME],
        initial_state,
        0.0,
        propagation_setup.integrator.runge_kutta_fixed_step(
            30.0,
            coefficient_set=propagation_setup.integrator.CoefficientSets.rkf_78,
        ),
        propagation_setup.propagator.time_termination(
            final_time, terminate_exactly_on_final_condition=True
        ),
    )


def propagate(bodies, acceleration_models, initial_state, final_time):
    dynamics_simulator = simulator.create_dynamics_simulator(
        bodies,
        create_propagator_settings(acceleration_models, initial_state, final_time),
    )
    return dynamics_simulator.state_history[max(dynamics_simulator.state_history)]


def test_custom_acceleration_rejects_all_non_keplerian_accelerations():
    spice.load_standard_kernels()

    body_settings = environment_setup.get_default_body_settings(
        [CENTRAL_BODY_NAME, "Sun", "Moon", "Jupiter"],
        CENTRAL_BODY_NAME,
        "J2000",
    )
    body_settings.add_empty_settings(SPACECRAFT_NAME)
    spacecraft_body_settings = body_settings.get(SPACECRAFT_NAME)
    spacecraft_body_settings.aerodynamic_coefficient_settings = (
        environment_setup.aerodynamic_coefficients.constant(0.035, [1.2, 0.0, 0.0])
    )
    spacecraft_body_settings.radiation_pressure_target_settings = (
        environment_setup.radiation_pressure.cannonball_radiation_target(
            0.035, 1.2, {"Sun": [CENTRAL_BODY_NAME]}
        )
    )
    bodies = environment_setup.create_system_of_bodies(body_settings)
    bodies.get(SPACECRAFT_NAME).mass = 2.2

    earth_gravitational_parameter = bodies.get(CENTRAL_BODY_NAME).gravitational_parameter
    semi_major_axis = 7.0e6
    initial_keplerian_state = np.array([semi_major_axis, 0.05, 0.4, 0.3, 0.2, 0.1])
    initial_state = element_conversion.keplerian_to_cartesian(
        initial_keplerian_state, earth_gravitational_parameter
    )
    orbital_period = 2.0 * np.pi * np.sqrt(semi_major_axis**3 / earth_gravitational_parameter)
    analytical_final_keplerian_state = two_body_dynamics.propagate_kepler_orbit(
        initial_keplerian_state,
        orbital_period,
        earth_gravitational_parameter,
    )
    analytical_final_cartesian_state = element_conversion.keplerian_to_cartesian(
        analytical_final_keplerian_state, earth_gravitational_parameter
    )

    rejecting_acceleration = PerturbationRejectingAcceleration(
        bodies.get(SPACECRAFT_NAME), bodies.get(CENTRAL_BODY_NAME)
    )
    acceleration_settings = {
        CENTRAL_BODY_NAME: [
            propagation_setup.acceleration.spherical_harmonic_gravity(8, 8),
            propagation_setup.acceleration.aerodynamic(),
        ],
        "Sun": [
            propagation_setup.acceleration.radiation_pressure(),
            propagation_setup.acceleration.point_mass_gravity(),
        ],
        "Moon": [propagation_setup.acceleration.point_mass_gravity()],
        "Jupiter": [propagation_setup.acceleration.point_mass_gravity()],
        SPACECRAFT_NAME: [
            propagation_setup.acceleration.custom_acceleration(rejecting_acceleration.acceleration)
        ],
    }
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies,
        {SPACECRAFT_NAME: acceleration_settings},
        [SPACECRAFT_NAME],
        [CENTRAL_BODY_NAME],
    )

    accelerations_to_cancel = [
        acceleration_model
        for body_exerting_acceleration, body_acceleration_models in acceleration_models[
            SPACECRAFT_NAME
        ].items()
        if body_exerting_acceleration != SPACECRAFT_NAME
        for acceleration_model in body_acceleration_models
    ]
    assert len(accelerations_to_cancel) == 6
    rejecting_acceleration.set_acceleration_models(accelerations_to_cancel)

    perturbed_acceleration_models = {
        SPACECRAFT_NAME: {
            body_exerting_acceleration: models
            for body_exerting_acceleration, models in acceleration_models[SPACECRAFT_NAME].items()
            if body_exerting_acceleration != SPACECRAFT_NAME
        }
    }
    perturbed_final_state = propagate(
        bodies, perturbed_acceleration_models, initial_state, orbital_period
    )
    corrected_final_state = propagate(bodies, acceleration_models, initial_state, orbital_period)

    perturbed_position_error = np.linalg.norm(
        perturbed_final_state[:3] - analytical_final_cartesian_state[:3]
    )
    perturbed_velocity_error = np.linalg.norm(
        perturbed_final_state[3:] - analytical_final_cartesian_state[3:]
    )
    corrected_position_error = np.linalg.norm(
        corrected_final_state[:3] - analytical_final_cartesian_state[:3]
    )
    corrected_velocity_error = np.linalg.norm(
        corrected_final_state[3:] - analytical_final_cartesian_state[3:]
    )

    assert corrected_position_error < POSITION_TOLERANCE_METERS
    assert corrected_velocity_error < VELOCITY_TOLERANCE_METERS_PER_SECOND

    assert perturbed_position_error > MINIMUM_PERTURBED_ERROR_FACTOR * POSITION_TOLERANCE_METERS
    assert (
        perturbed_velocity_error
        > MINIMUM_PERTURBED_ERROR_FACTOR * VELOCITY_TOLERANCE_METERS_PER_SECOND
    )
