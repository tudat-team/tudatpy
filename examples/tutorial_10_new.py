###############################################################################
# IMPORT STATEMENTS ###########################################################
###############################################################################
import numpy as np
import time
from tudatpy import elements
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.simulation import environment_setup
from tudatpy.kernel.simulation import propagation_setup
from tudatpy.kernel.simulation import estimation_setup

def main():

    # Load spice kernels.
    spice_interface.load_standard_kernels()

    # Set simulation start epoch.
    simulation_start_epoch = 0.0

    # Set simulation end epoch.
    simulation_end_epoch = constants.JULIAN_DAY

    ###########################################################################
    # CREATE ENVIRONMENT ######################################################
    ###########################################################################

    # Create body objects.
    bodies_to_create = ["Sun", "Earth", "Moon", "Mars", "Venus"]

    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create,
        simulation_start_epoch - 300.0,
        simulation_end_epoch + 300.0,
        "Earth","J2000")

    bodies = environment_setup.create_bodies(body_settings)

    ###########################################################################
    # CREATE VEHICLE ##########################################################
    ###########################################################################

    # Create vehicle objects.
    bodies.create_body( "Delfi-C3" )
    bodies.get( "Delfi-C3").set_constant_body_mass(400.0)

    ###########################################################################
    # CREATE VEHICLE - ENVIRONMENT INTERFACE ##################################
    ###########################################################################

    # Create aerodynamic coefficient interface settings
    reference_area = 4.0
    drag_coefficient = 1.2
    aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
        reference_area,[drag_coefficient,0,0],
        are_coefficients_in_aerodynamic_frame=True,
        are_coefficients_in_negative_axis_direction=True
    )
    environment_setup.add_aerodynamic_coefficient_interface(
                bodies, "Delfi-C3", aero_coefficient_settings );


    # Create radiation pressure settings
    reference_area_radiation = 4.0
    radiation_pressure_coefficient = 1.2
    occulting_bodies = ["Earth"]
    radiation_pressure_settings = environment_setup.radiation_pressure.cannonball(
        "Sun", reference_area_radiation, radiation_pressure_coefficient, occulting_bodies
    )
    environment_setup.add_radiation_pressure_interface(
                bodies, "Delfi-C3", radiation_pressure_settings );

    ###########################################################################
    # CREATE ACCELERATIONS ####################################################
    ###########################################################################

    # Define bodies that are propagated.
    bodies_to_propagate = ["Delfi-C3"]

    # Define central bodies.
    central_bodies = ["Earth"]

    # Define unique (Sun, Earth) accelerations acting on Delfi-C3.
    accelerations_settings_delfi_c3 = dict(
        Sun=
        [
            propagation_setup.acceleration.cannon_ball_radiation_pressure_(),
            propagation_setup.acceleration.point_mass_gravity()
        ],
        Earth=
        [
            propagation_setup.acceleration.spherical_harmonic_gravity(5, 5),
            propagation_setup.acceleration.aerodynamic()
        ])

    # Define other point mass accelerations acting on Delfi-C3.
    for other in set(bodies_to_create).difference({"Sun", "Earth"}):
        accelerations_settings_delfi_c3[other] = [
            propagation_setup.acceleration.point_mass_gravity()]

    # Create global accelerations dictionary.
    acceleration_settings = {"Delfi-C3": accelerations_settings_delfi_c3}

    # Create acceleration models.
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies,
        acceleration_settings,
        bodies_to_propagate,
        central_bodies)

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
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        system_initial_state,
        simulation_end_epoch
    )

    parameter_settings = estimation_setup.parameter.initial_states(
	propagator_settings, bodies )

    # Create numerical integrator settings.
    integrator_settings = propagation_setup.integrator.runge_kutta_4(
        simulation_start_epoch,
        10.0
    )

    ###########################################################################
    # PROPAGATE ORBIT #########################################################
    ###########################################################################

    # Create simulation object and propagate dynamics.
    variational_equations_solver = estimation_setup.SingleArcVariationalEquationsSolver(
        bodies, integrator_settings, propagator_settings, estimation_setup.create_parameters_to_estimate( parameter_settings, bodies ),
        integrate_on_creation=1 )

    equations_of_motion_result = variational_equations_solver.get_equations_of_motion_solution()
    state_transition_result = variational_equations_solver.get_state_transition_matrix_solution()
    sensitivity_result = variational_equations_solver.get_sensitivity_matrix_solution()

    ###########################################################################
    # PRINT INITIAL AND FINAL STATES ##########################################
    ###########################################################################

    print(
        f"""

The initial position vector of Delfi-C3 is [km]: \n{
        state_transition_result[simulation_start_epoch]}
The initial velocity vector of Delfi-C3 is [km/s]: \n{
        state_transition_result[simulation_end_epoch] }
        """
    )

    # Final statement (not required, though good practice in a __main__).
    return 0


if __name__ == "__main__":
    main()
