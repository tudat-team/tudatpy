###############################################################################
# IMPORT STATEMENTS ###########################################################
###############################################################################
import numpy as np
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.simulation import environment_setup
from tudatpy.kernel.simulation import propagation_setup
from tudatpy.kernel.astro import conversion
from tudatpy.kernel import example


def system_of_bodies_error_catch(body_system):
    get_body_fn = getattr(body_system, "get_body")

    def get_body_fn_wrapper(func):
        def catch(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except IndexError:
                raise IndexError("Body name is not recognized in the system of bodies.")

        return catch

    setattr(body_system, "get_body", get_body_fn_wrapper(get_body_fn))


def system_of_bodies_getattr(body_system):
    get_body_fn = getattr(body_system, "get_body")
    setattr(body_system, "__getattr__", get_body_fn)


system_of_bodies_error_catch(environment_setup.SystemOfBodies)
system_of_bodies_getattr(environment_setup.SystemOfBodies)


def main():
    # Load spice kernels.
    spice_interface.load_standard_kernels()

    # Set simulation start epoch.
    simulation_start_epoch = 0.0

    ###########################################################################
    # CREATE ENVIRONMENT ######################################################
    ###########################################################################

    # Create default body settings for "Earth"
    bodies_to_create = ["Earth"]

    # Create default body settings for bodies_to_create, with "Earth"/"J2000" as
    # global frame origin and orientation
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, "Earth", "J2000"
    )

    # Create Earth Object.
    bodies = environment_setup.create_system_of_bodies(body_settings)

    ###########################################################################
    # CREATE VEHICLE ##########################################################
    ###########################################################################

    # Create vehicle object.
    bodies.create_empty_body("Apollo")
    bodies.Apollo.set_constant_mass(5.0e3)

    # Add predefined aerodynamic coefficient database to the body
    bodies.Apollo.set_aerodynamic_coefficient_interface(
        example.apollo_aerodynamics_coefficient_interface()
    )

    ###########################################################################
    # CREATE ACCELERATIONS ####################################################
    ###########################################################################

    # Define bodies that are propagated.
    bodies_to_propagate = ["Apollo"]

    # Define central bodies.
    central_bodies = ["Earth"]

    # Define accelerations acting on Apollo.
    accelerations_settings_apollo = dict(
        Earth=[
            propagation_setup.acceleration.spherical_harmonic_gravity(4, 0),
            propagation_setup.acceleration.aerodynamic(),
        ]
    )
    acceleration_settings = {"Apollo": accelerations_settings_apollo}

    # Create acceleration models.
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, acceleration_settings, bodies_to_propagate, central_bodies
    )

    # Define constant 30 degree angle of attack.
    constant_angle_of_attack = np.deg2rad(30.0)
    environment_setup.set_aerodynamic_orientation_functions(
        body=bodies.Apollo,
        angle_of_attack_function=lambda: constant_angle_of_attack,
    )

    ###########################################################################
    # CREATE PROPAGATION SETTINGS #############################################
    ###########################################################################

    # Set spherical elements for Apollo and convert to Cartesian.
    initial_radial_distance = (
        bodies.get_body("Earth").shape_model.average_radius + 120.0e3
    )
    initial_earth_fixed_state = conversion.spherical_to_cartesian(
        radial_distance=initial_radial_distance,
        latitude=np.deg2rad(0.0),
        longitude=np.deg2rad(68.75),
        speed=7.7e3,
        flight_path_angle=np.deg2rad(-0.9),
        heading_angle=np.deg2rad(34.37),
    )

    # Convert the state from Earth-fixed to inertial frame
    earth_rotation_model = bodies.get_body("Earth").rotation_model
    initial_state = conversion.transform_to_inertial_orientation(
        initial_earth_fixed_state, simulation_start_epoch, earth_rotation_model
    )

    # Define list of dependent variables to save.
    dependent_variables_to_save = [
        propagation_setup.dependent_variable.mach_number("Apollo", "Earth"),
        propagation_setup.dependent_variable.altitude("Apollo", "Earth"),
        propagation_setup.dependent_variable.single_acceleration(
            propagation_setup.acceleration.aerodynamic_type, "Apollo", "Earth"
        ),
        propagation_setup.dependent_variable.aerodynamic_force_coefficients("Apollo"),
    ]

    # Define termination conditions (once altitude goes below 25 km).
    termination_variable = propagation_setup.dependent_variable.altitude(
        "Apollo", "Earth"
    )
    termination_settings = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings=termination_variable,
        limit_value=25.0e3,
        use_as_lower_limit=True,
        terminate_exactly_on_final_condition=False,
    )

    # Create propagation settings.
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        initial_state,
        termination_settings,
        output_variables=dependent_variables_to_save,
    )

    # Create numerical integrator settings.
    fixed_step_size = 1.0
    integrator_settings = propagation_setup.integrator.runge_kutta_4(
        simulation_start_epoch, fixed_step_size
    )

    ###########################################################################
    # PROPAGATE ORBIT #########################################################
    ###########################################################################

    # Create simulation object and propagate dynamics.
    dynamics_simulator = propagation_setup.SingleArcDynamicsSimulator(
        bodies, integrator_settings, propagator_settings
    )
    states = dynamics_simulator.state_history
    dependent_variables = dynamics_simulator.dependent_variable_history

    # from tudatpy
    # array = result2array(states)

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
