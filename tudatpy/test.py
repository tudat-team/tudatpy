import numpy as np

from core import spice_interface

spice_interface.load_standard_spice_kernels()

# 1.1. Set Up the Environment
from core.simulation_setup import get_default_body_settings
from core.simulation_setup import create_bodies
from core.simulation_setup import ConstantEphemerisSettings

bodies_to_create = ["Earth"]
body_settings = get_default_body_settings(bodies_to_create)
body_settings["Earth"].ephemeris_settings = ConstantEphemerisSettings(
    np.zeros(6))
body_dict = create_bodies(body_settings)

# 1.2. Create the Vehicle
from core.simulation_setup import Body
from core.simulation_setup import set_global_frame_body_ephemerides

body_dict["Asterix"] = Body()
set_global_frame_body_ephemerides(body_dict, "SSB", "ECLIPJ2000")
print(body_dict)

# 1.3. Set Up the Acceleration models
from core.basic_astrodynamics import AvailabileAcceleration
from core.simulation_setup import AccelerationSettings
from core.simulation_setup import create_acceleration_models_dict
from core.orbital_element_conversions import KeplerianElementIndices
from core.orbital_element_conversions import convert_keplerian_to_cartesian_elements
from core.propagators import TranslationalStatePropagatorSettings

bodies_to_propagate = ["Asterix"]
central_bodies = ["Earth"]
accelerations_of_asterix = {
    "Earth": [
        AccelerationSettings(AvailabileAcceleration.point_mass_gravity)
    ]}
acceleration_dict = {"Asterix": accelerations_of_asterix}

acceleration_model_dict = create_acceleration_models_dict(
    body_dict,
    acceleration_dict,
    bodies_to_propagate,
    central_bodies)

# 1.4. Set Up the Propagation Settings
KEI = KeplerianElementIndices
asterix_initial_state_in_keplerian_elements = np.zeros(6)
kep_state = asterix_initial_state_in_keplerian_elements
kep_state[int(KEI.semi_major_axis_index)] = 7500.0E3
kep_state[int(KEI.eccentricity_index)] = 0.1
kep_state[int(KEI.inclination_index)] = np.deg2rad(85.3)
kep_state[int(KEI.argument_of_periapsis_index)] = np.deg2rad(235.7)
kep_state[int(KEI.longitude_of_ascending_node_index)] = np.deg2rad(23.4)
kep_state[int(KEI.true_anomaly_index)] = np.deg2rad(139.87)

earth_gravitational_parameter = body_dict[
    "Earth"].gravity_field_model.get_gravitational_parameter()
system_initial_state = convert_keplerian_to_cartesian_elements(
    kep_state, earth_gravitational_parameter)

from core.constants import JULIAN_DAY as simulation_end_epoch

propagator_settings = TranslationalStatePropagatorSettings(
    central_bodies, acceleration_model_dict, bodies_to_propagate,
    system_initial_state,
    simulation_end_epoch
)

# Create numerical integrator settings.
from core.numerical_integrators import IntegratorSettings, AvailableIntegrators

simulation_start_epoch = 0.0
fixed_step_size = 10.0
integrator_settings = IntegratorSettings(AvailableIntegrators.rungeKutta4,
                                         simulation_start_epoch, fixed_step_size);

# 1.5. Perform the Orbit Propagation
from core.propagators import SingleArcDynamicsSimulator

dynamics_simulator = SingleArcDynamicsSimulator(body_dict, integrator_settings, propagator_settings)
integration_result = dynamics_simulator.get_equations_of_motion_numerical_solution()

print(integration_result)
# NOTE: In original tutorial, julian day is not defined. The code snippets
# should be self contained in nature.

# print(system_initial_state)
