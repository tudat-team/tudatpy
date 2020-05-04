###############################################################################
# IMPORT STATEMENTS ###########################################################
###############################################################################
import numpy as np
from tudatpy import constants
from tudatpy import elements
from tudatpy import io
from tudatpy import ephemerides
from tudatpy import interpolators
from tudatpy import numerical_integrators
from tudatpy import spice_interface
from tudatpy import basic_astrodynamics
# from tudatpy import orbital_element_conversions # LEGACY MODULE
from tudatpy import propagators
from tudatpy import aerodynamics
from tudatpy import simulation_setup
from tudatpy.prototype import Environment


def main():
    # Load spice kernels.
    spice_interface.load_standard_spice_kernels()

    # Set simulation start epoch.
    simulation_start_epoch = 0.0

    # Set numerical integration fixed step size.
    fixed_step_size = 1.0

    # Set simulation end epoch.
    simulation_end_epoch = constants.JULIAN_DAY

    ###########################################################################
    # CREATE ENVIRONMENT ######################################################
    ###########################################################################

    bodies_to_create = ["Sun", "Earth", "Moon", "Mars", "Venus"]

    env = Environment(
        bodies_to_create,
        user_defined_bodies=None,
        frame_origin="SSB",
        frame_orientation="ECLIPJ2000",
        custom_body_settings=None,
        alternative_kernels=None,
        start_epoch=simulation_start_epoch,
        end_epoch=simulation_end_epoch,
        fixed_dt=fixed_step_size,
        epoch_margin=300.0)

    ###########################################################################
    # CREATE VEHICLE ##########################################################
    ###########################################################################

    env.vehicles = ["Delfi-C3"]

    env.bodies["Delfi-C3"].set_constant_body_mass(400.0)

    ###########################################################################
    # CREATE VEHICLE - ENVIRONMENT INTERFACE ##################################
    ###########################################################################

    # Create aerodynamic interface settings for Delfi-C3.
    env.set_aerodynamic_coefficient_interface(
        target_body_name="Delfi-C3",
        reference_area=4.0,
        coefficients=1.2 * np.ones(3),
        coefficients_in_aerodynamic_frame=True,
        coefficients_in_negative_axis_direction=True,
        model="constant_coefficients")

    # Create radiation pressure interface settings for Delfi-C3.
    env.set_radiation_pressure_interface(
        target_body_name="Delfi-C3",
        reference_area=4.0,
        coefficient=1.2,
        occulting=["Earth"],
        source="Sun",
        model="cannon_ball")

    ###########################################################################
    # FINALIZE BODIES #########################################################
    ###########################################################################

    # Finalize body creation.
    env.finalize_environment(frame_origin="SSB", frame_orientation="J2000")

    ###########################################################################
    # CREATE ACCELERATIONS ####################################################
    ###########################################################################

    # Define bodies that are propagated.
    bodies_to_propagate = env.vehicles

    # Define central bodies.
    central_bodies = ["Earth"]

    # Define unique (Sun, Earth) accelerations acting on Delfi-C3.
    accelerations_of_delfi_c3 = dict(
        Sun=
        [
            simulation_setup.Acceleration.canon_ball_radiation_pressure()
            # AccelerationSettings(AvailableAcceleration.cannon_ball_radiation_pressure) # LEGACY DESIGN.
        ],
        Earth=
        [
            simulation_setup.Acceleration.spherical_harmonic_gravity(5, 5),
            # SphericalHarmonicAccelerationSettings(5, 5), # LEGACY DESIGN.

            simulation_setup.Acceleration.aerodynamic()
            # AccelerationSettings(AvailableAcceleration.aerodynamic) # LEGACY DESIGN.
        ])

    # Define other point mass accelerations acting on Delfi-C3.
    for other in set(bodies_to_create).difference({"Sun", "Earth"}):
        accelerations_of_delfi_c3[other] = [
            simulation_setup.Acceleration.point_mass_gravity()]

    # Create global accelerations dictionary.
    acceleration_dict = {"Delfi-C3": accelerations_of_delfi_c3}

    # Create acceleration models.
    acceleration_models = simulation_setup.create_acceleration_models_dict(
        bodies_to_create,
        acceleration_dict,
        bodies_to_propagate,
        central_bodies)

    ###########################################################################
    # CREATE PROPAGATION SETTINGS #############################################
    ###########################################################################

    # Set initial conditions for the Asterix satellite that will be
    # propagated in this simulation. The initial conditions are given in
    # Keplerian elements and later on converted to Cartesian elements.

    # Set Keplerian elements for Delfi-C3
    earth_gravitational_parameter = env.bodies[
        "Earth"].gravity_field_model.get_gravitational_parameter()

    # LEGACY DESIGN.
    # KEI = orbital_element_conversions.KeplerianElementIndices
    # asterix_initial_state_in_keplerian_elements = np.zeros(6)
    # kep_state = asterix_initial_state_in_keplerian_elements
    # kep_state[int(KEI.semi_major_axis_index)] = 7500.0E3
    # kep_state[int(KEI.eccentricity_index)] = 0.1
    # kep_state[int(KEI.inclination_index)] = np.deg2rad(85.3)
    # kep_state[int(KEI.argument_of_periapsis_index)] = np.deg2rad(235.7)
    # kep_state[int(KEI.longitude_of_ascending_node_index)] = np.deg2rad(23.4)
    # kep_state[int(KEI.true_anomaly_index)] = np.deg2rad(139.87)
    # system_initial_state = corbital_element_conversions.onvert_keplerian_to_cartesian_elements(
    #     kep_state, earth_gravitational_parameter)

    # REVISED CONTEMPORARY DESIGN.
    system_initial_state = elements.keplerian2cartesian(
        mu=earth_gravitational_parameter,
        a=7500.0E3,
        ecc=0.1,
        inc=np.deg2rad(85.3),
        raan=np.deg2rad(23.4),
        argp=np.deg2rad(235.7),
        nu=np.deg2rad(139.87))

    # Create propagation settings.
    propagator_settings = propagators.TranslationalStatePropagatorSettings(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        system_initial_state,
        simulation_end_epoch
    )
    # Create numerical integrator settings.
    integrator_settings = numerical_integrators.IntegratorSettings(
        numerical_integrators.AvailableIntegrators.rungeKutta4,
        simulation_start_epoch,
        fixed_step_size
    )

    ###########################################################################
    # PROPAGATE ORBIT #########################################################
    ###########################################################################

    # Create simulation object and propagate dynamics.
    dynamics_simulator = propagators.SingleArcDynamicsSimulator(
        bodies, integrator_settings, propagator_settings, True)
    result = dynamics_simulator.get_equations_of_motion_numerical_solution()

    ###########################################################################
    # PRINT INITIAL AND FINAL STATES ##########################################
    ###########################################################################

    print(
        f"""
Single Earth-Orbiting Satellite Example.
The initial position vector of Delfi-C3 is [km]: \n{
        result[simulation_start_epoch][:3] / 1E3}
The initial velocity vector of Delfi-C3 is [km]: \n{
        result[simulation_start_epoch][3:] / 1E3}
After {simulation_end_epoch} seconds the position vector of Delfi-C3 is [km]: \n{
        result[simulation_end_epoch][:3] / 1E3}
And the velocity vector of Delfi-C3 is [km]: \n{
        result[simulation_start_epoch][3:] / 1E3}
        """
    )

    ###########################################################################
    # SAVE RESULTS ############################################################
    ###########################################################################

    io.save2txt(
        solution=result,
        filename="singlePerturbedSatellitePropagationHistory.dat",
        directory="./tutorial_2",
    )

    # Final statement (not required, though good practice in a __main__).
    return 0


if __name__ == "__main__":
    main()
