###############################################################################
# IMPORT STATEMENTS ###########################################################
###############################################################################
import numpy as np
from tudatpy import elements
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.simulation import environment_setup
from tudatpy.kernel.simulation import propagation_setup


def main():
    # Load spice kernels.
    spice_interface.load_standard_kernels()

    # Set simulation start epoch.
    simulation_start_epoch = 0.0

    # Set numerical integration fixed step size.
    fixed_step_size = 10.0

    # Set simulation end epoch.
    simulation_end_epoch = constants.JULIAN_DAY

    ###########################################################################
    # CREATE ENVIRONMENT ######################################################
    ###########################################################################

    # Create body objects.
    bodies_to_create = ["Earth"]

    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, "Earth", "J2000")

    # Create Earth Object.
    bodies = environment_setup.create_bodies(body_settings)

    ###########################################################################
    # CREATE VEHICLE ##########################################################
    ###########################################################################

    # Create vehicle objects.
    bodies.add_new_body( "Delfi-C3", 1)

    ###########################################################################
    # CREATE ACCELERATIONS ####################################################
    ###########################################################################

    # Define bodies that are propagated.
    bodies_to_propagate = ["Delfi-C3"]

    # Define central bodies.
    central_bodies = ["Earth"]

    # Define accelerations acting on Delfi-C3.
    acceleration_settings_delfi_c3 = dict(

        Earth=[propagation_setup.acceleration.point_mass_gravity( ) ]
    )

    # Create global accelerations dictionary.
    acceleration_settings = {"Delfi-C3": acceleration_settings_delfi_c3}

    # Create acceleration models.
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, acceleration_settings, bodies_to_propagate, central_bodies)

    ###########################################################################
    # CREATE PROPAGATION SETTINGS #############################################
    ###########################################################################

    # Set initial conditions for the Asterix satellite that will be
    # propagated in this simulation. The initial conditions are given in
    # Keplerian elements and later on converted to Cartesian elements.

    # Set Keplerian elements for Delfi-C3
    earth_gravitational_parameter = environment_setup.get_body_gravitational_parameter( 
	bodies, "Earth" )

    # REVISED CONTEMPORARY DESIGN.
    system_initial_state = elements.keplerian2cartesian(
        mu=earth_gravitational_parameter,
        sma=7500.0E3,
        ecc=0.1,
        inc=np.deg2rad(85.3),
        raan=np.deg2rad(23.4),
        argp=np.deg2rad(235.7),
        theta=np.deg2rad(139.87)
    )

    # Create propagation settings.
    propagator_settings = propagation_setup.TranslationalStatePropagatorSettings(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        system_initial_state,
        simulation_end_epoch
    )
    # Create numerical integrator settings.
    integrator_settings = propagation_setup.integrator.runge_kutta_4(
        simulation_start_epoch,
        fixed_step_size
    )

    ###########################################################################
    # PROPAGATE ORBIT #########################################################
    ###########################################################################

    # Create simulation object and propagate dynamics.
    dynamics_simulator = propagation_setup.SingleArcDynamicsSimulator(
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
    #
    # io.save2txt(
    #     solution=result,
    #     filename="singleSatellitePropagationHistory.dat",
    #     directory="./tutorial_1",
    # )

    # Final statement (not required, though good practice in a __main__).
    return 0


if __name__ == "__main__":
    main()
