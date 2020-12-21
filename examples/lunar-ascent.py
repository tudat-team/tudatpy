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
simulation_start_epoch = 1.0e7

# Set numerical integration fixed step size.
fixed_step_size = 3600.0

# Set simulation end epoch.
simulation_end_epoch = 1.0e7 + 5.0 * constants.JULIAN_YEAR

# Set vehicle mass.
vehicle_mass = 4.7E3
vehicle_dry_mass = 2.25E3
specific_impulse = 311.0
maximum_duration = 86400.0
termination_altitude = 100.0E3

# Set vehicle thrust magnitude.
thrust_magnitude = 25.0

# Set vehicle specific impulse.


thrust_parameters = np.array(
    [15629.13262285292,
     21.50263026822358,
     -0.03344538412056863,
     -0.06456210720352829,
     0.3943447499535977,
     0.5358478897251189,
     -0.8607350478880107])

################################################################################
# SETUP ENVIRONMENT ############################################################
################################################################################

# Define bodies in simulation.
bodies_to_create = ["Sun", "Earth", "Moon"]

# Create bodies in simulation.
body_settings = environment_setup.get_default_body_settings(bodies_to_create)
system_of_bodies = environment_setup.create_system_of_bodies(body_settings)

################################################################################
# SETUP ENVIRONMENT : CREATE VEHICLE ###########################################
################################################################################

vehicle_name = "Vehicle"
system_of_bodies.create_empty_body(vehicle_name)
system_of_bodies.get_body(vehicle_name).set_constant_mass(vehicle_mass)

################################################################################
# SETUP PROPAGATION : DEFINE THRUST GUIDANCE SETTINGS ##########################
################################################################################

class LunarAscentGuidance:
    def __init__(self, vehicle_body, initial_time, parameter_vector):
        self.vehicle_body = vehicle_body
        self.initial_time = initial_time
        self.parameter_vector = parameter_vector
        self.thrust_magnitude = parameter_vector[0]
        self.time_interval = parameter_vector[1]
        self.thrust_angle_dict = {}

        current_time = initial_time
        for i in range(len(parameter_vector) - 2):
            self.thrust_angle_dict[current_time] = parameter_vector[i + 2]
            current_time = self.time_interval
        self.thrust_angle_interpolator =  #

    def get_current_thrust_direction(self, time):
        angle = self.thrust_angle_interpolator(time)
        thrust_direction_vertical_frame = np.array([0, np.sin(angle), - np.cos(angle)])
        vertical_to_intertial_frame = self.vehicle_body.flight_conditions.aerodynamic_angle_calculator.get_rotation_quaternion_between_frames(
        # propagation_setup. < ENUM? >.vertical_frame,
        # propagation_setup. < ENUM? >.inertial_frame)
        return vertical_to_intertial_frame * thrust_direction_vertical_frame

    def get_current_thrust_magnitude(self, time):
        return lambda time: self.thrust_magnitude
        # return self.thrust_magnitude


def get_thrust_acceleration_model_from_parameters(
        thrust_parameters,
        system_of_bodies,
        initial_time,
        specific_impulse
):
    #include<pybind11/functional.h>
    thrust_guidance = LunarAscentGuidance(system_of_bodies[vehicle_name], initial_time, thrust_parameters)
    thrust_direction_function = thrust_guidance.get_current_thrust_direction
    thrust_magnitude_function = thrust_guidance.get_current_thrust_magnitude
    thrust_direction_settings = propagation_setup.CustomThrustDirectionSettings(thrust_direction_function)
    thrust_magnitude_settings = propagation_setup.FromFunctionThrustMagnitudeSettings(thrust_magnitude_function, lambda time: specific_impulse)
    return propagation_setup.acceleration.ThrustAccelerationSettings(thrust_direction_settings, thrust_magnitude_settings)


################################################################################
# SETUP PROPAGATION : CREATE ACCELERATION MODELS ###############################
################################################################################


acceleration_on_vehicle = dict(
    Vehicle=[
        get_thrust_acceleration_model_from_parameters(thrust_parameters,
                                                      system_of_bodies,
                                                      simulation_start_epoch,
                                                      specific_impulse)
    ],
    Earth=[propagation_setup.acceleration.point_mass_gravity()],
    Moon=[propagation_setup.acceleration.point_mass_gravity()],
    Sun=[propagation_setup.acceleration.point_mass_gravity()],
)

bodies_to_propagate = ["Vehicle"]

central_bodies = ["Earth"]

acceleration_dict = dict(Vehicle=acceleration_on_vehicle)

# Convert acceleration mappings into acceleration models.
acceleration_models = propagation_setup.create_acceleration_models(
    system_of_bodies=system_of_bodies,
    selected_acceleration_per_body=acceleration_dict,
    bodies_to_propagate=bodies_to_propagate,
    central_bodies=central_bodies,
)

################################################################################
# SETUP PROPAGATION : PROPAGATION SETTINGS #####################################
################################################################################

# Get gravitational parameter of Earth for initial state.
# gravitational_parameter = system_of_bodies["Earth"].gravity_field_model.get_gravitational_parameter()

# Get system initial state.
system_initial_state = np.array([8.0e6, 0, 0, 0, 7.5e3, 0])
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
propagator_settings = propagation_setup.propagator.translational(
    central_bodies,
    acceleration_models,
    bodies_to_propagate,
    system_initial_state,
    simulation_end_epoch,
)
# Create numerical integrator settings.
integrator_settings = propagation_setup.integrator.runge_kutta_4(
    simulation_start_epoch, fixed_step_size
)

################################################################################
# PROPAGATE ####################################################################
################################################################################

# Instantiate the dynamics simulator.
dynamics_simulator = propagation_setup.SingleArcDynamicsSimulator(
    system_of_bodies, integrator_settings, propagator_settings, True
)

# Propagate and store results to outer loop results dictionary.
result = dynamics_simulator.get_equations_of_motion_numerical_solution()

################################################################################
# VISUALISATION / OUTPUT / PRELIMINARY ANALYSIS ################################
################################################################################

array = result2array(result)
print(array)
