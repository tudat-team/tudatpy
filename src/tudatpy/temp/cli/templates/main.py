################################################################################
# IMPORT STATEMENTS ############################################################
################################################################################
import numpy as np
from tudatpy import elements
from tudatpy.util import result2array
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.simulation import environment_setup
from tudatpy.kernel.simulation import propagation_setup

################################################################################
# GENERAL SIMULATION SETUP #####################################################
################################################################################

# Load spice kernels.
spice_interface.load_standard_kernels()

# Set simulation start epoch.
simulation_start_epoch = 1.0E7

# Set numerical integration fixed step size.
fixed_step_size = 3600.0

# Set simulation end epoch.
simulation_end_epoch = 1.0E7 + 5.0 * constants.JULIAN_YEAR

# Set vehicle mass.
vehicle_mass = 5.0E3

# Set vehicle thrust magnitude.
thrust_magnitude = 25.0

# Set vehicle specific impulse.
specific_impulse = 5.0E3

################################################################################
# SETUP ENVIRONMENT ############################################################
################################################################################

# Define bodies in simulation.
bodies_to_create = ["Sun", "Earth", "Moon"]

# Create bodies in simulation.
body_settings = environment_setup.get_default_body_settings(bodies_to_create)
body_system = environment_setup.create_bodies(body_settings)
environment_setup.set_global_frame_body_ephemerides(body_system, "SSB", "ECLIPJ2000")

################################################################################
# SETUP ENVIRONMENT : CREATE VEHICLE ###########################################
################################################################################

body_system["Vehicle"] = environment_setup.Body()
body_system["Vehicle"].set_constant_body_mass(vehicle_mass)

################################################################################
# SETUP ENVIRONMENT : FINALIZE BODY CREATION ###################################
################################################################################

environment_setup.set_global_frame_body_ephemerides(body_system, "SSB", "ECLIPJ2000")

################################################################################
# SETUP PROPAGATION : DEFINE THRUST GUIDANCE SETTINGS ##########################
################################################################################

thrust_direction_settings = propagation_setup.ThrustDirectionFromStateGuidanceSettings(
    central_body="Earth",
    is_colinear_with_velocity=True,
    direction_is_opposite_to_vector=False
)

thrust_magnitude_settings = propagation_setup.ConstantThrustMagnitudeSettings(
    thrust_magnitude=thrust_magnitude,
    specific_impulse=specific_impulse
)

################################################################################
# SETUP PROPAGATION : CREATE ACCELERATION MODELS ###############################
################################################################################

acceleration_on_vehicle = dict(
    Vehicle=[
        propagation_setup.ThrustAccelerationSettings(
            thrust_direction_settings=thrust_direction_settings,
            thrust_magnitude_settings=thrust_magnitude_settings)
    ],
    Earth=[
        propagation_setup.AccelerationSettings(
            propagation_setup.AvailableAcceleration.point_mass_gravity)
    ],
    Moon=[
        propagation_setup.AccelerationSettings(
            propagation_setup.AvailableAcceleration.point_mass_gravity)
    ],
    Sun=[
        propagation_setup.AccelerationSettings(
            propagation_setup.AvailableAcceleration.point_mass_gravity)
    ]
)

bodies_to_propagate = ["Vehicle"]

central_bodies = ["Earth"]

acceleration_dict = dict(Vehicle=acceleration_on_vehicle)

# Convert acceleration mappings into acceleration models.
acceleration_models = propagation_setup.create_acceleration_models_dict(
    body_system=body_system,
    selected_acceleration_per_body=acceleration_dict,
    bodies_to_propagate=bodies_to_propagate,
    central_bodies=central_bodies
)

################################################################################
# SETUP PROPAGATION : PROPAGATION SETTINGS #####################################
################################################################################

# Get gravitational parameter of Earth for initial state.
# gravitational_parameter = body_system["Earth"].gravity_field_model.get_gravitational_parameter()

# Get system initial state.
system_initial_state = np.array([8.0E6, 0, 0, 0, 7.5E3, 0])
# system_initial_state = elements.keplerian2cartesian(
#     mu=gravitational_parameter,
#     sma=8.0E6,
#     ecc=0.1,
#     inc=np.deg2rad(0.05),
#     raan=np.deg2rad(0.1),
#     argp=np.deg2rad(0.1),
#     theta=np.deg2rad(0.1)
# )

# termination_settings = propagation_setup.

# Create propagation settings.
propagator_settings = propagation_setup.TranslationalStatePropagatorSettings(
    central_bodies,
    acceleration_models,
    bodies_to_propagate,
    system_initial_state,
    simulation_end_epoch
)
# Create numerical integrator settings.
integrator_settings = propagation_setup.IntegratorSettings(
    propagation_setup.AvailableIntegrators.rk4,
    simulation_start_epoch,
    fixed_step_size
)

################################################################################
# PROPAGATE ####################################################################
################################################################################

# Instantiate the dynamics simulator.
dynamics_simulator = propagation_setup.SingleArcDynamicsSimulator(
    body_system, integrator_settings, propagator_settings, True)

# Propagate and store results to outer loop results dictionary.
result = dynamics_simulator.get_equations_of_motion_numerical_solution()

################################################################################
# VISUALISATION / OUTPUT / PRELIMINARY ANALYSIS ################################
################################################################################

import matplotlib.pyplot as plt
from tudatpy.util import result2array

array = result2array(result)

print(array)
