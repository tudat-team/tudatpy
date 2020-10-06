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

    ###########################################################################
    # CREATE ENVIRONMENT ######################################################
    ###########################################################################

    bodies_to_create = ["Earth"]

    body_settings = environment_setup.get_default_body_settings(
	bodies_to_create,"Earth","J2000")

    # Create Earth Object.
    bodies = environment_setup.create_bodies(body_settings)

    ###########################################################################
    # CREATE VEHICLE ##########################################################
    ###########################################################################

    # Create vehicle objects.
    bodies.create_body( "Apollo" )
    bodies.get( "Apollo" ).set_aerodynamic_coefficient_interface(
        example.apollo_aerodynamics_coefficient_interface())

    bodies.get( "Apollo" ).set_constant_body_mass(5.0E3)

    ###########################################################################
    # CREATE ACCELERATIONS ####################################################
    ###########################################################################

    # Define bodies that are propagated.
    bodies_to_propagate = ["Apollo"]

    # Define central bodies.
    central_bodies = ["Earth"]

    # Define accelerations acting on Apollo.
    accelerations_settings_apollo = dict(Earth=[
        propagation_setup.acceleration.spherical_harmonic_gravity(4, 0),
        propagation_setup.acceleration.aerodynamic()
    ])

    # Create global accelerations dictionary.
    acceleration_settings = {"Apollo": accelerations_settings_apollo}

    # Create acceleration models.
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, acceleration_settings,
        bodies_to_propagate, central_bodies)

    # Define constant 30 degree angle of attack.
    constant_angle_of_attack = np.deg2rad(30.0)
    bodies.get("Apollo").get_flight_conditions(
    ).get_aerodynamic_angle_calculator(
    ).set_orientation_angle_functions(lambda: constant_angle_of_attack)

    ###########################################################################
    # CREATE PROPAGATION SETTINGS #############################################
    ###########################################################################

    # Set spherical elements for Apollo and convert to Cartesian.
    # REVISED CONTEMPORARY DESIGN
    cartesian_initial_state_earth_fixed = elements.spherical2cartesian(
        r=spice_interface.get_average_radius("Earth") + 120.0E3,
        lat=np.deg2rad(0.0),
        lon=np.deg2rad(68.75),
        speed=7.7E3,
        fpa=np.deg2rad(-0.9),
        heading=np.deg2rad(34.37))

    # Convert the state to the global (inertial) frame.
    earth_rotational_ephemeris = bodies.get("Earth").get_rotational_ephemeris()
    system_initial_state = ephemerides.transform_state_to_global_frame(
        cartesian_initial_state_earth_fixed,
        simulation_start_epoch,
        earth_rotational_ephemeris)

    # Define list of dependent variables to save.
    # TODO: Revise design of dependent variable saves with Python class layer.
    dependent_variables_to_save = propagation_setup.dependent_variables.create(
        [
        propagation_setup.dependent_variables.mach_number(
            "Apollo", "Earth"
        ),
        propagation_setup.dependent_variables.altitude(
            "Apollo", "Earth"
        ),

        propagation_setup.dependent_variables.single_acceleration(
            propagation_setup.acceleration.aerodynamic_type, "Apollo", "Earth"
        )#,

        #propagation_setup.SingleDependentVariableSaveSettings(
        #    propagation_setup.PropagationDependentVariables.aerodynamic_force_coefficients_dependent_variable,
        #    "Apollo"
        #)
        ]
    )
    # Define termination conditions.
    termination_variable = propagation_setup.dependent_variables.altitude(
        "Apollo", "Earth"
    )
    termination_settings = propagation_setup.propagator.termination(
        dependent_variadble_settings=termination_variable,
        limit_value=25.0E3,
        use_as_lower_limit=True,
        terminate_exactly_on_final_condition=False
    )

    # Create propagation settings.
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        system_initial_state,
        termination_settings,
        propagation_setup.propagator.cowell,
        dependent_variables_to_save
    )
    integrator_settings = propagation_setup.integrator.runge_kutta_4(
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
    state_solution=dynamics_simulator.get_equations_of_motion_numerical_solution()
    dependent_variable_solution=dynamics_simulator.get_dependent_variable_history()

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
