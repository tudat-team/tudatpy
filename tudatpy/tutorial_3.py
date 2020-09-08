###############################################################################
# IMPORT STATEMENTS ###########################################################
###############################################################################
import numpy as np
from tudatpy import elements
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.astro import ephemerides
from tudatpy.kernel.simulation import environment_setup
from tudatpy.kernel.simulation import propagation_setup
from tudatpy.kernel import example


def main():
    # Load spice kernels.
    spice_interface.load_standard_kernels()

    # Set simulation start epoch.
    simulation_start_epoch = 0.0

    # Set numerical integration fixed step size.
    fixed_step_size = 10.0

    # Set simulation end epoch.
    # simulation_end_epoch = 3100.0 # UNUSED VARIABLE

    ###########################################################################
    # CREATE ENVIRONMENT ######################################################
    ###########################################################################

    bodies_to_create = ["Earth"]

    body_settings = environment_setup.get_default_body_settings(bodies_to_create)

    body_settings["Earth"].ephemeris_settings = environment_setup.ConstantEphemerisSettings(
        np.zeros(6), "SSB", "J2000")

    body_settings["Earth"].rotation_model_settings.reset_original_frame("J2000")

    # Create Earth Object.
    bodies = environment_setup.create_bodies(body_settings)

    ###########################################################################
    # CREATE VEHICLE ##########################################################
    ###########################################################################

    # Create vehicle objects.
    bodies["Apollo"] = environment_setup.Body()

    bodies["Apollo"].set_aerodynamic_coefficient_interface(
        example.apollo_aerodynamics_coefficient_interface())

    bodies["Apollo"].set_constant_body_mass(5.0E3)

    ###########################################################################
    # FINALIZE BODIES #########################################################
    ###########################################################################

    # Finalize body creation.
    environment_setup.set_global_frame_body_ephemerides(bodies, "SSB", "J2000")

    ###########################################################################
    # CREATE ACCELERATIONS ####################################################
    ###########################################################################

    # Define bodies that are propagated.
    bodies_to_propagate = ["Apollo"]

    # Define central bodies.
    central_bodies = ["Earth"]

    # Define accelerations acting on Apollo.
    accelerations_of_apollo = dict(Earth=[
        propagation_setup.Acceleration.spherical_harmonic_gravity(4, 0),
        propagation_setup.Acceleration.aerodynamic()
    ])

    # Create global accelerations dictionary.
    acceleration_dict = dict(Apollo=accelerations_of_apollo)

    # Create acceleration models.
    acceleration_models = propagation_setup.create_acceleration_models_dict(
        bodies, acceleration_dict,
        bodies_to_propagate, central_bodies)

    # Define constant 30 degree angle of attack.
    constant_aoa = np.deg2rad(30.0)
    bodies["Apollo"].get_flight_conditions(
    ).get_aerodynamic_angle_calculator(
    ).set_orientation_angle_functions(lambda: constant_aoa)

    ###########################################################################
    # CREATE PROPAGATION SETTINGS #############################################
    ###########################################################################

    # Set spherical elements for Apollo and convert to Cartesian.
    # REVISED CONTEMPORARY DESIGN
    cartesian_initial_state = elements.spherical2cartesian(
        r=spice_interface.get_average_radius("Earth") + 120.0E3,
        lat=np.deg2rad(0.0),
        lon=np.deg2rad(68.75),
        speed=7.7E3,
        fpa=np.deg2rad(-0.9),
        heading=np.deg2rad(34.37))

    # Convert the state to the global (inertial) frame.
    earth_rotational_ephemeris = bodies["Earth"].get_rotational_ephemeris()
    system_initial_state = ephemerides.transform_state_to_global_frame(
        cartesian_initial_state,
        simulation_start_epoch,
        earth_rotational_ephemeris)

    # Define list of dependent variables to save.
    # TODO: Revise design of dependent variable saves with Python class layer.
    dependent_variables_list = [
        propagation_setup.SingleDependentVariableSaveSettings(
            propagation_setup.PropagationDependentVariables.mach_number_dependent_variable,
            "Apollo"
        ),
        propagation_setup.SingleDependentVariableSaveSettings(
            propagation_setup.PropagationDependentVariables.altitude_dependent_variable,
            "Apollo", "Earth"
        ),
        propagation_setup.SingleAccelerationDependentVariableSaveSettings(
            propagation_setup.AvailableAcceleration.aerodynamic,
            "Apollo", "Earth", 1
        ),
        propagation_setup.SingleDependentVariableSaveSettings(
            propagation_setup.PropagationDependentVariables.aerodynamic_force_coefficients_dependent_variable,
            "Apollo"
        )
    ]
    dependent_variables_to_save = propagation_setup.DependentVariableSaveSettings(
        dependent_variables_list
    )
    # Define termination conditions.
    termination_dependent_variable = propagation_setup.SingleDependentVariableSaveSettings(
        propagation_setup.PropagationDependentVariables.altitude_dependent_variable,
        "Apollo", "Earth"
    )
    termination_settings = propagation_setup.PropagationDependentVariableTerminationSettings(
        dependent_variable_settings=termination_dependent_variable,
        limit_value=25.0E3,
        use_as_lower_limit=True,
        terminate_exactly_on_final_condition=False
    )

    # Create propagation settings.
    propagator_settings = propagation_setup.TranslationalStatePropagatorSettings(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        system_initial_state,
        termination_settings,
        propagation_setup.TranslationalPropagatorType.cowell,
        dependent_variables_to_save
    )
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
        bodies, integrator_settings, propagator_settings
    )

    # io.save2txt(
    #     solution=dynamics_simulator.get_equations_of_motion_numerical_solution(),
    #     filename="apolloPropagationHistory.dat",
    #     directory="./tutorial_3",
    # )
    # io.save2txt(
    #     solution=dynamics_simulator.get_dependent_variable_history(),
    #     filename="apolloDependentVariableHistory.dat",
    #     directory="./tutorial_3"
    # )

    # Final statement (not required, though good practice in a __main__).
    return 0


if __name__ == "__main__":
    main()
