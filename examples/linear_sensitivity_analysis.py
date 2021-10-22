
###############################################################################
# IMPORT STATEMENTS ###########################################################
###############################################################################
import numpy as np
from tudatpy.kernel import constants
from tudatpy.kernel import numerical_simulation
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.numerical_simulation import estimation_setup
from matplotlib import pyplot as plt

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

    global_frame_origin = "Earth"
    global_frame_orientation = "J2000"
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation )

    bodies = environment_setup.create_system_of_bodies(body_settings)

    ###########################################################################
    # CREATE VEHICLE ##########################################################
    ###########################################################################

    # Create vehicle objects.
    bodies.create_empty_body( "Delfi-C3" )
    bodies.get( "Delfi-C3").mass = 400.0

    ###########################################################################
    # CREATE VEHICLE - ENVIRONMENT INTERFACE ##################################
    ###########################################################################

    # Create aerodynamic coefficient interface settings
    reference_area = 4.0
    drag_coefficient = 1.2
    aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
        reference_area,[drag_coefficient,0,0]
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
            propagation_setup.acceleration.cannonball_radiation_pressure(),
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
    earth_gravitational_parameter = bodies.get( "Earth" ).gravitational_parameter
    initial_state = element_conversion.keplerian_to_cartesian_elementwise(
        gravitational_parameter=earth_gravitational_parameter,
        semi_major_axis=7500.0E3,
        eccentricity=0.1,
        inclination=np.deg2rad(85.3),
        argument_of_periapsis=np.deg2rad(235.7),
        longitude_of_ascending_node=np.deg2rad(23.4),
        true_anomaly=np.deg2rad(139.87)
    )
    # Create propagation settings.
    termination_condition = propagation_setup.propagator.time_termination( simulation_end_epoch )
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        initial_state,
        termination_condition
    )

    # Create list of parameters for which the variational equations are to be
    # propagated
    parameter_settings = estimation_setup.parameter.initial_states(
	propagator_settings, bodies )
    parameter_settings.append(estimation_setup.parameter.gravitational_parameter("Earth"))
    parameter_settings.append(estimation_setup.parameter.constant_drag_coefficient("Delfi-C3"))

    # Create numerical integrator settings.
    fixed_step_size = 10.0
    integrator_settings = propagation_setup.integrator.runge_kutta_4(
        simulation_start_epoch,
        fixed_step_size
    )

    ###########################################################################
    # PROPAGATE ORBIT #########################################################
    ###########################################################################

    # Create simulation object and propagate dynamics.
    variational_equations_solver = numerical_simulation.SingleArcVariationalSimulator(
        bodies, integrator_settings, propagator_settings, estimation_setup.create_parameters_to_estimate( parameter_settings, bodies ),
        integrate_on_creation=1 )

    states = variational_equations_solver.state_history
    state_transition_matrices = variational_equations_solver.state_transition_matrix_history
    sensitivity_matrices = variational_equations_solver.sensitivity_matrix_history

    ###########################################################################
    # PERFORM SENSITIVITY ANALYSIS    #########################################
    ###########################################################################

    initial_state_variation = [1, 0, 0, 1.0E-3, 0, 0]
    earth_standard_param_variation = [-2.0E+5, 0.0]
    drag_coeff_variation = [0.0, 0.05]

    delta_initial_state_dict = dict()
    earth_standard_param_dict = dict()
    delta_drag_coeff_dict = dict()

    for epoch in state_transition_matrices:
        delta_initial_state_dict[epoch] = np.dot(state_transition_matrices[epoch], initial_state_variation)
        earth_standard_param_dict[epoch] = np.dot(sensitivity_matrices[epoch], earth_standard_param_variation)
        delta_drag_coeff_dict[epoch] = np.dot(sensitivity_matrices[epoch], drag_coeff_variation)

    ###########################################################################
    # PLOT RESULTS   ##########################################################
    ###########################################################################

    font_size = 20
    plt.rcParams.update({'font.size': font_size})

    time = state_transition_matrices.keys()
    time_hours = [t / 3600 for t in time]

    delta_initial_state = np.vstack(list(delta_initial_state_dict.values()))
    delta_earth_standard_param = np.vstack(list(earth_standard_param_dict.values()))
    delta_drag_coefficient = np.vstack(list(delta_drag_coeff_dict.values()))

    # 1 // due to initial state variation
    delta_r1 = np.linalg.norm( delta_initial_state[:, 0:3], axis = 1 )
    delta_v1 = np.linalg.norm( delta_initial_state[:, 3:7], axis = 1 )

    # 2 // due to gravitational parameter variation
    delta_r2 = np.linalg.norm( delta_earth_standard_param[:, 0:3], axis = 1 )
    delta_v2 = np.linalg.norm( delta_earth_standard_param[:, 3:7], axis = 1 )

    # 3 // due to drag coefficient variation
    delta_r3 = np.linalg.norm( delta_drag_coefficient[:, 0:3], axis = 1 )
    delta_v3 = np.linalg.norm( delta_drag_coefficient[:, 3:7], axis = 1 )

    # Plot deviations of position
    plt.figure(figsize=(17, 5))
    plt.grid()
    plt.plot(time_hours, delta_r1, color='tomato', label='variation initial state')
    plt.plot(time_hours, delta_r2, color='orange', label='variation grav. parameter (Earth)')
    plt.plot(time_hours, delta_r3, color='cyan', label='variation drag coefficient')
    plt.yscale('log')
    plt.xlabel('Time [hr]')
    plt.ylabel('$\Delta r (t_1)$ [m]')
    plt.xlim([min(time_hours), max(time_hours)])
    plt.legend()

    # Plot deviations of speed
    plt.figure(figsize=(17, 5))
    plt.grid()
    plt.plot(time_hours, delta_v1, color='tomato', label='variation initial state')
    plt.plot(time_hours, delta_v2, color='orange', label='variation grav. parameter (Earth)')
    plt.plot(time_hours, delta_v3, color='cyan', label='variation drag coefficient')
    plt.yscale('log')
    plt.xlabel('Time [hr]')
    plt.ylabel('$\Delta v (t_1)$ [m/s]')
    plt.xlim([min(time_hours), max(time_hours)])
    plt.legend()

    plt.show()
    # Final statement (not required, though good practice in a __main__).
    return 0


if __name__ == "__main__":
    main()
