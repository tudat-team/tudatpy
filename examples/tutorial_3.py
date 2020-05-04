###############################################################################
# IMPORT STATEMENTS ###########################################################
###############################################################################
import numpy as np
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
from tudatpy import unit_tests


def main():
    # Load spice kernels.
    spice_interface.load_standard_spice_kernels()

    # Set simulation start epoch.
    simulation_start_epoch = 0.0

    # Set numerical integration fixed step size.
    fixed_step_size = 1.0

    # Set simulation end epoch.
    # simulation_end_epoch = 3100.0 # UNUSED VARIABLE

    ###########################################################################
    # CREATE ENVIRONMENT ######################################################
    ###########################################################################

    bodies_to_create = ["Earth"]

    body_settings = simulation_setup.get_default_body_settings(bodies_to_create)

    body_settings["Earth"].ephemeris_settings = simulation_setup.ConstantEphemerisSettings(
        np.zeros(6), "SSB", "J2000")

    body_settings["Earth"].rotation_model_settings.reset_original_frame("J2000")

    # Create Earth Object.
    bodies = simulation_setup.create_bodies(body_settings)

    ###########################################################################
    # CREATE VEHICLE ##########################################################
    ###########################################################################

    # Create vehicle objects.
    bodies["Apollo"] = simulation_setup.Body()

    bodies["Apollo"].set_aerodynamic_coefficient_interface(unit_tests.get_apollo_coefficient_interface())

    bodies["Apollo"].set_constant_body_mass(5.0E3)

    ###########################################################################
    # FINALIZE BODIES #########################################################
    ###########################################################################

    # Finalize body creation.
    simulation_setup.set_global_frame_body_ephemerides(bodies, "SSB", "J2000")

    ###########################################################################
    # CREATE ACCELERATIONS ####################################################
    ###########################################################################

    # Define bodies that are propagated.
    bodies_to_propagate = ["Apollo"]

    # Define central bodies.
    central_bodies = ["Earth"]

    # Define accelerations acting on Apollo.
    accelerations_of_apollo = dict(Earth=[
        simulation_setup.Acceleration.spherical_harmonic_gravity(4, 0),
        simulation_setup.Acceleration.aerodynamic()
    ])

    # Create global accelerations dictionary.
    acceleration_dict = dict(Apollo=accelerations_of_apollo)

    # Create acceleration models.
    acceleration_models = simulation_setup.create_acceleration_models_dict(
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

    # LEGACY DESIGN
    # spherical_idx = orbital_element_conversions.SphericalOrbitalStateElementIndices
    # apollo_spherical_entry_state = np.zeros(6)
    # apollo_spherical_entry_state[spherical_idx.radiusIndex] = spice_interface.get_average_radius("Earth") + 120.0E3
    # apollo_spherical_entry_state[spherical_idx.latitudeIndex] = np.deg2rad(0.0)
    # apollo_spherical_entry_state[spherical_idx.longitudeIndex] = np.deg2rad(68.75)
    # apollo_spherical_entry_state[spherical_idx.speedIndex] = 7.7E3
    # apollo_spherical_entry_state[spherical_idx.flightPathIndex] = np.deg2rad(-0.9)
    # apollo_spherical_entry_state[spherical_idx.headingAngleIndex] = np.deg2rad(34.37)
    # system_initial_state = orbital_element_conversions.convert_spherical_orbital_to_cartesian_state(
    #     apollo_spherical_entry_state)

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
        propagators.SingleDependentVariableSaveSettings(
            propagators.PropagationDependentVariables.mach_number_dependent_variable,
            "Apollo"
        ),
        propagators.SingleDependentVariableSaveSettings(
            propagators.PropagationDependentVariables.altitude_dependent_variable,
            "Apollo", "Earth"
        ),
        propagators.SingleAccelerationDependentVariableSaveSettings(
            basic_astrodynamics.AvailableAcceleration.aerodynamic,
            "Apollo", "Earth", 1
        ),
        propagators.SingleDependentVariableSaveSettings(
            propagators.PropagationDependentVariables.aerodynamic_force_coefficients_dependent_variable,
            "Apollo"
        )
    ]
    dependent_variables_to_save = propagators.DependentVariableSaveSettings(
        dependent_variables_list)

    # Define termination conditions.
    termination_dependent_variable = propagators.SingleDependentVariableSaveSettings(
        propagators.PropagationDependentVariables.altitude_dependent_variable,
        "Apollo", "Earth"
    )
    termination_settings = propagators.PropagationDependentVariableTerminationSettings(
        dependent_variable_settings=termination_dependent_variable,
        limit_value=25.0E3,
        use_as_lower_limit=True,
        terminate_exactly_on_final_condition=False
    )

    # Create propagation settings.
    propagator_settings = propagators.TranslationalStatePropagatorSettings(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        system_initial_state,
        termination_settings,
        propagators.TranslationalPropagatorType.cowell,
        dependent_variables_to_save
    )
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
        bodies, integrator_settings, propagator_settings)

    io.save2txt(
        solution=dynamics_simulator.get_equations_of_motion_numerical_solution(),
        filename="apolloPropagationHistory.dat",
        directory="./tutorial_3",
    )
    io.save2txt(
        solution=dynamics_simulator.get_dependent_variable_history(),
        filename="apolloDependentVariableHistory.dat",
        directory="./tutorial_3"
    )

    # Final statement (not required, though good practice in a __main__).
    return 0


if __name__ == "__main__":
    main()
