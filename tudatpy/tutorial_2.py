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
    bodies_to_create = ["Sun", "Earth", "Moon", "Mars", "Venus"]

    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create,
        simulation_start_epoch - 300.0,
        simulation_end_epoch + 300.0,
        fixed_step_size)

    for body in bodies_to_create:
        body_settings[body].ephemeris_settings.reset_frame_orientation("J2000")
        body_settings[body].rotation_model_settings.reset_original_frame("J2000")

    bodies = environment_setup.create_bodies(body_settings)

    ###########################################################################
    # CREATE VEHICLE ##########################################################
    ###########################################################################

    # Create vehicle objects.
    bodies["Delfi-C3"] = environment_setup.Body()
    bodies["Delfi-C3"].set_constant_body_mass(400.0)

    ###########################################################################
    # CREATE VEHICLE - ENVIRONMENT INTERFACE ##################################
    ###########################################################################

    # Create aerodynamic coefficient interface settings
    reference_area = 4.0
    aerodynamic_coefficient = 1.2
    aero_c_settings = environment_setup.ConstantAerodynamicCoefficientSettings(
        reference_area,
        aerodynamic_coefficient * np.ones(3),
        are_coefficients_in_aerodynamic_frame=True,
        are_coefficients_in_negative_axis_direction=True
    )
    # Create and set aerodynamic coefficients object
    bodies["Delfi-C3"].set_aerodynamic_coefficient_interface(
        environment_setup.create_aerodynamic_coefficient_interface(
            aero_c_settings,
            "Delfi-C3")
    )
    # TODO: Simplify (post 1.0.0 work)

    # Create radiation pressure settings
    reference_area_radiation = 4.0
    radiation_pressure_coefficient = 1.2
    occulting_bodies = ["Earth"]
    rad_press_settings = environment_setup.CannonBallRadiationPressureInterfaceSettings(
        "Sun", reference_area_radiation, radiation_pressure_coefficient, occulting_bodies)

    # Create and set radiation pressure settings
    bodies["Delfi-C3"].set_radiation_pressure_interface(
        "Sun", environment_setup.create_radiation_pressure_interface(
            rad_press_settings, "Delfi-C3", bodies))

    # TODO: Simplify (post 1.0.0 work)

    ###########################################################################
    # FINALIZE BODIES #########################################################
    ###########################################################################

    # Finalize body creation.
    environment_setup.set_global_frame_body_ephemerides(bodies, "SSB", "J2000")

    ###########################################################################
    # CREATE ACCELERATIONS ####################################################
    ###########################################################################

    # Define bodies that are propagated.
    bodies_to_propagate = ["Delfi-C3"]

    # Define central bodies.
    central_bodies = ["Earth"]

    # Define unique (Sun, Earth) accelerations acting on Delfi-C3.
    accelerations_of_delfi_c3 = dict(
        Sun=
        [
            propagation_setup.Acceleration.canon_ball_radiation_pressure()
            # AccelerationSettings(AvailableAcceleration.cannon_ball_radiation_pressure) # LEGACY DESIGN.
        ],
        Earth=
        [
            propagation_setup.Acceleration.spherical_harmonic_gravity(5, 5),
            # SphericalHarmonicAccelerationSettings(5, 5), # LEGACY DESIGN.

            propagation_setup.Acceleration.aerodynamic()
            # AccelerationSettings(AvailableAcceleration.aerodynamic) # LEGACY DESIGN.
        ])

    # Define other point mass accelerations acting on Delfi-C3.
    for other in set(bodies_to_create).difference({"Sun", "Earth"}):
        accelerations_of_delfi_c3[other] = [
            propagation_setup.Acceleration.point_mass_gravity()]

    # Create global accelerations dictionary.
    acceleration_dict = {"Delfi-C3": accelerations_of_delfi_c3}

    # Create acceleration models.
    acceleration_models = propagation_setup.create_acceleration_models_dict(
        bodies,
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
    earth_gravitational_parameter = bodies[
        "Earth"].gravity_field_model.get_gravitational_parameter()

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
    integrator_settings = propagation_setup.IntegratorSettings(
        propagation_setup.AvailableIntegrators.rungeKutta4,
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
The initial velocity vector of Delfi-C3 is [km/s]: \n{
        result[simulation_start_epoch][3:] / 1E3}
After {simulation_end_epoch} seconds the position vector of Delfi-C3 is [km]: \n{
        result[simulation_end_epoch][:3] / 1E3}
And the velocity vector of Delfi-C3 is [km/s]: \n{
        result[simulation_end_epoch][3:] / 1E3}
        """
    )

    # Final statement (not required, though good practice in a __main__).
    return 0


if __name__ == "__main__":
    main()
