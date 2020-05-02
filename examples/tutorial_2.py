# Initial imports and kernel setup
import numpy as np
import tudatpy as tpy
from tudatpy import spice_interface

spice_interface.load_standard_spice_kernels()

# 1.1. Set Up the Environment
from tudatpy.simulation_setup import get_default_body_settings
from tudatpy.simulation_setup import set_global_frame_body_ephemerides
from tudatpy.simulation_setup import create_bodies
from tudatpy.simulation_setup import ConstantAerodynamicCoefficientSettings
from tudatpy.simulation_setup import CannonBallRadiationPressureInterfaceSettings
from tudatpy.simulation_setup import create_aerodynamic_coefficient_interface
from tudatpy.simulation_setup import create_radiation_pressure_interface
from tudatpy.constants import JULIAN_DAY

simulation_end_epoch = JULIAN_DAY

bodies_to_create = ["Sun", "Earth", "Moon", "Mars", "Venus"]

# Create body objects.
simulation_start_epoch = 0.0
time_step = 10.0

body_settings = get_default_body_settings(bodies_to_create,
                                          simulation_start_epoch - 300.0,
                                          simulation_end_epoch + 300.0, time_step)

# TODO: Simplify
for body in bodies_to_create:
    body_settings[body].ephemeris_settings.reset_frame_orientation("J2000")
    body_settings[body].rotation_model_settings.reset_original_frame("J2000")

body_dict = create_bodies(body_settings)

# 2.2. Create the Vehicle
from tudatpy.simulation_setup import Body

body_dict["Asterix"] = Body()
body_dict["Asterix"].set_constant_body_mass(400.0)
set_global_frame_body_ephemerides(body_dict, "SSB", "J2000")

# Create aerodynamic coefficient interface settings
reference_area = 4.0
aerodynamic_coefficient = 1.2
aerodynamic_coefficient_settings = ConstantAerodynamicCoefficientSettings(
    reference_area,
    aerodynamic_coefficient * np.ones(3), 1, 1
)

# Create and set aerodynamic coefficients object
body_dict["Asterix"].set_aerodynamic_coefficient_interface(
    create_aerodynamic_coefficient_interface(aerodynamic_coefficient_settings,
                                             "Asterix")
)

# Create radiation pressure settings
reference_area_radiation = 4.0
radiation_pressure_coefficient = 1.2

occulting_bodies = ["Earth"]
asterix_radiation_pressure_settings = CannonBallRadiationPressureInterfaceSettings(
    "Sun", reference_area_radiation, radiation_pressure_coefficient, occulting_bodies)

# Create and set radiation pressure settings
body_dict["Asterix"].set_radiation_pressure_interface(
    "Sun", create_radiation_pressure_interface(
        asterix_radiation_pressure_settings, "Asterix", body_dict))

# TODO: Simplify (post 1.0.0 work)

# 2.3. Set Up the Acceleration Models
# from tudatpy.basic_astrodynamics import AvailableAcceleration
# from tudatpy.simulation_setup import AccelerationSettings
from tudatpy.simulation_setup import Acceleration
# from tudatpy.simulation_setup import SphericalHarmonicAccelerationSettings
from tudatpy.simulation_setup import create_acceleration_models_dict
from tudatpy.orbital_element_conversions import KeplerianElementIndices
from tudatpy.orbital_element_conversions import convert_keplerian_to_cartesian_elements
from tudatpy.propagators import TranslationalStatePropagatorSettings

# Define propagation settings.
accelerations_of_asterix = dict(
    Sun=
    [
        Acceleration.canon_ball_radiation_pressure()
        # AccelerationSettings(AvailableAcceleration.cannon_ball_radiation_pressure)
    ],
    Earth=
    [
        Acceleration.spherical_harmonic_gravity(5, 5),
        # SphericalHarmonicAccelerationSettings(5, 5),

        Acceleration.aerodynamic()
        # AccelerationSettings(AvailableAcceleration.aerodynamic)
    ])

for other in set(bodies_to_create).difference({"Sun", "Earth"}):
    accelerations_of_asterix[other] = [Acceleration.point_mass_gravity()]

acceleration_dict = dict(Asterix=accelerations_of_asterix)

bodies_to_propagate = ["Asterix"]
central_bodies = ["Earth"]

acceleration_model_dict = create_acceleration_models_dict(
    body_dict,
    acceleration_dict,
    bodies_to_propagate,
    central_bodies)

# 1.4. Set Up the Propagation Settings
from tudatpy.constants import JULIAN_DAY as simulation_end_epoch

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

propagator_settings = TranslationalStatePropagatorSettings(
    central_bodies, acceleration_model_dict, bodies_to_propagate,
    system_initial_state,
    simulation_end_epoch
)

# Create numerical integrator settings.
from tudatpy.numerical_integrators import IntegratorSettings, AvailableIntegrators

simulation_start_epoch = 0.0
fixed_step_size = 10.0
integrator_settings = IntegratorSettings(AvailableIntegrators.rungeKutta4,
                                         simulation_start_epoch, fixed_step_size)

# 1.5. Perform the Orbit Propagation
from tudatpy.propagators import SingleArcDynamicsSimulator

dynamics_simulator = SingleArcDynamicsSimulator(body_dict, integrator_settings, propagator_settings)
integration_result = dynamics_simulator.get_equations_of_motion_numerical_solution()

import pandas as pd

df = pd.DataFrame(index=integration_result.keys(),
                  data=np.vstack(list(integration_result.values())),
                  columns="x y z Vx Vy Vz".split(" "))

df.index.name = "time"
df.to_csv("results_2.dat")

pd.set_option('display.max_rows', 30)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

print(df)
