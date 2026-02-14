# Load required standard modules
import numpy as np
from matplotlib import pyplot as plt

# Load required tudatpy modules
from tudatpy import constants
from tudatpy.interface import spice
from tudatpy import dynamics
from tudatpy.dynamics import environment, environment_setup
from tudatpy.dynamics import propagation_setup, parameters_setup, simulator
from tudatpy import estimation
from tudatpy.estimation import observable_models_setup, observable_models, observations_setup, observations, estimation_analysis
from tudatpy.astro.time_representation import DateTime
from tudatpy.astro import element_conversion
from tudatpy.astro import frame_conversion

# Load spice kernels
spice.load_standard_kernels()

# Set simulation start and end epochs
simulation_start_epoch = DateTime(2024, 8, 28).to_epoch()
simulation_end_epoch   = DateTime(2024, 9, 5).to_epoch()
observation_start_epoch_1 = DateTime(2024, 8, 30).to_epoch()
observation_end_epoch_1   = DateTime(2024, 9, 1).to_epoch()
observation_start_epoch_2   = DateTime(2024, 9, 3).to_epoch()
observation_end_epoch_2   = DateTime(2024, 9, 4).to_epoch()

observation_start_epoch_3   = DateTime(2024, 8, 29).to_epoch()
observation_end_epoch_3 = simulation_end_epoch

# Create default body settings for "Sun", "Earth", "Moon", "Mars", and "Venus"
bodies_to_create = ["Sun", "Earth", "Moon", "Mars", "Venus"]

# Create default body settings for bodies_to_create, with "Earth"/"J2000" as the global frame origin and orientation
global_frame_origin = "Earth"
global_frame_orientation = "J2000"
body_settings = environment_setup.get_default_body_settings(
    bodies_to_create, global_frame_origin, global_frame_orientation)

# Create system of bodies
bodies = environment_setup.create_system_of_bodies(body_settings)

# Create vehicle objects.
body_settings.add_empty_settings("Starlink-32101")
body_settings.get("Starlink-32101").constant_mass = 260

# Create aerodynamic coefficient interface settings
reference_area = 20  # Average projection area of a 3U CubeSat
drag_coefficient = 1.2
aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
    reference_area, [drag_coefficient, 0.0, 0.0]
)
# Add the aerodynamic interface to the environment
body_settings.get("Starlink-32101").aerodynamic_coefficient_settings = aero_coefficient_settings

# Create radiation pressure settings
reference_area_radiation = 20  # Average projection area of a 3U CubeSat
radiation_pressure_coefficient = 1.2
occulting_bodies = dict()
occulting_bodies["Sun"] = ["Earth"]
radiation_pressure_settings = environment_setup.radiation_pressure.cannonball_radiation_target(
    reference_area_radiation, radiation_pressure_coefficient, occulting_bodies)

# Add the radiation pressure interface to the environment
body_settings.get("Starlink-32101").radiation_pressure_target_settings = radiation_pressure_settings

# Create system of bodies
bodies = environment_setup.create_system_of_bodies(body_settings)

# Define bodies that are propagated
bodies_to_propagate = ["Starlink-32101"]

# Define central bodies of propagation
central_bodies = ["Earth"]

# Define the accelerations acting on `Starlink-32101`
accelerations_settings_Starlink_32101 = dict(
    Sun=[
        propagation_setup.acceleration.radiation_pressure(),
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Mars=[
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Moon=[
        propagation_setup.acceleration.point_mass_gravity()
    ],
    Earth=[
        propagation_setup.acceleration.spherical_harmonic_gravity(8, 8),
        propagation_setup.acceleration.aerodynamic()
    ])

# Create global accelerations dictionary
acceleration_settings = {"Starlink-32101": accelerations_settings_Starlink_32101}

# Create acceleration models
acceleration_models = propagation_setup.create_acceleration_models(
    bodies,
    acceleration_settings,
    bodies_to_propagate,
    central_bodies)

# Retrieve the initial state of `Starlink-32101` using Two-Line-Elements (TLEs)
Starlink_tle = environment_setup.ephemeris.sgp4(
    "1 60447U 24144Y   24239.91667824 -.00652022  00000-0 -25508-2 0  9990",
    "2 60447  53.1498 303.6008 0000548  88.4809  23.6264 15.87779028  3478",
)
Starlink_ephemeris = environment_setup.create_body_ephemeris(Starlink_tle, "Starlink-32101")
initial_state = Starlink_ephemeris.cartesian_state( simulation_start_epoch )

# Create numerical integrator settings
integrator_settings = propagation_setup.integrator.\
    runge_kutta_fixed_step_size(initial_time_step=60.0,
                                coefficient_set=propagation_setup.integrator.CoefficientSets.rkdp_87)

# Create termination settings
termination_condition = propagation_setup.propagator.time_termination(simulation_end_epoch)

# Create propagation settings
propagator_settings = propagation_setup.propagator.translational(
    central_bodies,
    acceleration_models,
    bodies_to_propagate,
    initial_state,
    simulation_start_epoch,
    integrator_settings,
    termination_condition
)

# Define the position of the ground station on Earth
station_altitude = 0.0
delft_latitude = np.deg2rad(52.00667)
delft_longitude = np.deg2rad(4.35556)

# Add the ground station to the environment
environment_setup.add_ground_station(
    bodies.get_body("Earth"),
    "TrackingStation",
    [station_altitude, delft_latitude, delft_longitude],
    element_conversion.geodetic_position_type)

# Define the uplink link ends for one-way observable
link_ends = dict()
link_ends[observable_models_setup.links.receiver] = observable_models_setup.links.body_reference_point_link_end_id("Earth", "TrackingStation")
link_ends[observable_models_setup.links.transmitter] = observable_models_setup.links.body_origin_link_end_id("Starlink-32101")

# Create observation settings for each link/observable
link_definition = observable_models_setup.links.LinkDefinition(link_ends)
observation_settings_list = [observable_models_setup.model_settings.one_way_doppler_instantaneous(link_definition)]

# 1 - Define Observation Simulation Settings

simulation_times = np.arange(simulation_start_epoch, simulation_end_epoch, 60.0)

observation_times_1 = np.arange(observation_start_epoch_1, observation_end_epoch_1, 60.0)
observation_times_2 = np.arange(observation_start_epoch_2, observation_end_epoch_2, 60.0)
observation_times_3 = np.arange(observation_start_epoch_3, observation_end_epoch_3, 60.0)

output_times = simulation_times
covariance_output_list = []
formal_errors_list = []
covariances_list = []

for n_scenario in [1,2,3]:
    print(f'Scenario n. {n_scenario}\n')
    print(f'Setting up Estimation and Propagating Covariance ...\n')
    if n_scenario == 1:
        observation_times = observation_times_1
    elif n_scenario == 2:
        observation_times = np.concatenate((observation_times_1, observation_times_2))
    else:
        observation_times = observation_times_3

    observation_simulation_settings = observations_setup.observations_simulation_settings.tabulated_simulation_settings(
        observable_models_setup.model_settings.one_way_instantaneous_doppler_type,
        link_definition,
        observation_times
    )

    # Add noise levels of roughly 1.0E-3 [m/s] and add this as Gaussian noise to the observation
    noise_level = 1.0E-3
    observations_setup.random_noise.add_gaussian_noise_to_observable(
        [observation_simulation_settings],
        noise_level,
        observable_models_setup.model_settings.one_way_instantaneous_doppler_type
    )

    # Create viability settings
    viability_setting = observations_setup.viability.elevation_angle_viability(["Earth", "TrackingStation"], np.deg2rad(15))
    observations_setup.viability.add_viability_check_to_all(
        [observation_simulation_settings],
        [viability_setting]
    )


# 2 - Define the Parameters to Estimate

    # Setup parameters settings to propagate the state transition matrix
    parameter_settings = parameters_setup.initial_states(propagator_settings, bodies)

    # Add estimated parameters to the sensitivity matrix that will be propagated
    parameter_settings.append(parameters_setup.gravitational_parameter("Earth"))
    parameter_settings.append(parameters_setup.constant_drag_coefficient("Starlink-32101"))

    # Create the parameters that will be estimated
    parameters_to_estimate = parameters_setup.create_parameter_set(parameter_settings, bodies)

    # Create the estimator
    estimator = estimation_analysis.Estimator(
        bodies,
        parameters_to_estimate,
        observation_settings_list,
        propagator_settings)

# 3 - Perform the observations simulation

    # Simulate required observations
    simulated_observations = observations_setup.observations_wrapper.simulate_observations(
        [observation_simulation_settings],
        estimator.observation_simulators,
        bodies)

# 4 - Define the Input Covariance

    # Define weighting of the observations in the inversion
    weights_per_observable = {observations.observations_processing.observation_parser(
        observable_models_setup.model_settings.one_way_instantaneous_doppler_type ): noise_level ** -2}
    simulated_observations.set_constant_weight_per_observation_parser( weights_per_observable )

    # Create input object for covariance analysis
    covariance_input = estimation_analysis.CovarianceAnalysisInput(
        simulated_observations)

    # Set methodological options
    covariance_input.define_covariance_settings(
        reintegrate_variational_equations=False)

    # Perform the covariance analysis
    covariance_output = estimator.compute_covariance(covariance_input)
    covariance_output_list.append(covariance_output)

    initial_covariance = covariance_output.covariance
    state_transition_interface = estimator.state_transition_interface

# 5 - Propagate the Covariances and the Formal Errors

    #Propagate the covariancees and the formal errors
    propagated_covariances = estimation_analysis.propagate_covariance_split_output(initial_covariance,state_transition_interface,output_times)
    # Propagate formal errors over the course of the orbit
    propagated_formal_errors = estimation_analysis.propagate_formal_errors_split_output(
        initial_covariance=initial_covariance,
        state_transition_interface=state_transition_interface,
        output_times=output_times)
    # Split tuple into epochs and formal errors
    epochs = np.array(propagated_formal_errors[0])
    formal_errors = np.array(propagated_formal_errors[1])
    formal_errors_list.append(formal_errors)

# 6 - Append Results

    covariances = np.array(propagated_covariances[1])
    covariances_list.append(covariances)
    print('... Done.\n')

print('All Done.\n')

