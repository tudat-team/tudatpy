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
from tudatpy.kernel.numerical_simulation import propagation
from matplotlib import pyplot as plt


def main():
    # Load spice kernels.
    spice_interface.load_standard_kernels()

    # Set simulation start and end epochs.
    simulation_start_epoch = 0.0
    simulation_end_epoch = constants.JULIAN_DAY

    ###########################################################################
    # CREATE ENVIRONMENT ######################################################
    ###########################################################################

    # Define string names for bodies to be created from default.
    bodies_to_create = ["Sun", "Earth", "Moon", "Mars", "Venus"]

    # Use "Earth"/"J2000" as global frame origin and orientation.
    global_frame_origin = "Earth"
    global_frame_orientation = "J2000"

    # Create default body settings, usually from `spice`.
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create,
        global_frame_origin,
        global_frame_orientation)

    # Create system of selected celestial bodies
    bodies = environment_setup.create_system_of_bodies(body_settings)

    # Create vehicle objects.
    bodies.create_empty_body("Delfi-C3")

    bodies.get("Delfi-C3").mass = 400.0

    # Create aerodynamic coefficient interface settings, and add to vehicle
    reference_area = 4.0
    drag_coefficient = 1.2
    aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
        reference_area, [drag_coefficient, 0, 0]
    )
    environment_setup.add_aerodynamic_coefficient_interface(
        bodies, "Delfi-C3", aero_coefficient_settings)

    # Create radiation pressure settings, and add to vehicle
    reference_area_radiation = 4.0
    radiation_pressure_coefficient = 1.2
    occulting_bodies = ["Earth"]
    radiation_pressure_settings = environment_setup.radiation_pressure.cannonball(
        "Sun", reference_area_radiation, radiation_pressure_coefficient, occulting_bodies
    )
    environment_setup.add_radiation_pressure_interface(
        bodies, "Delfi-C3", radiation_pressure_settings)

    ###########################################################################
    # CREATE ACCELERATIONS ####################################################
    ###########################################################################

    # Define bodies that are propagated.
    bodies_to_propagate = ["Delfi-C3"]

    # Define central bodies.
    central_bodies = ["Earth"]

    # Define accelerations acting on Delfi-C3 by Sun and Earth.
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
        ],
        Moon=
        [
            propagation_setup.acceleration.point_mass_gravity()
        ],
        Mars=
        [
            propagation_setup.acceleration.point_mass_gravity()
        ],
        Venus=
        [
            propagation_setup.acceleration.point_mass_gravity()
        ]
    )

    # Create global accelerations settings dictionary.
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
    earth_gravitational_parameter = bodies.get("Earth").gravitational_parameter
    initial_state = element_conversion.keplerian_to_cartesian_elementwise(
        gravitational_parameter=earth_gravitational_parameter,
        semi_major_axis=7500.0E3,
        eccentricity=0.1,
        inclination=np.deg2rad(85.3),
        argument_of_periapsis=np.deg2rad(235.7),
        longitude_of_ascending_node=np.deg2rad(23.4),
        true_anomaly=np.deg2rad(139.87)
    )

    # Define list of dependent variables to save.
    dependent_variables_to_save = [
        propagation_setup.dependent_variable.total_acceleration("Delfi-C3"),
        propagation_setup.dependent_variable.keplerian_state("Delfi-C3", "Earth"),
        propagation_setup.dependent_variable.latitude("Delfi-C3", "Earth"),
        propagation_setup.dependent_variable.longitude("Delfi-C3", "Earth"),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.point_mass_gravity_type, "Delfi-C3", "Sun"
        ),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.point_mass_gravity_type, "Delfi-C3", "Moon"
        ),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.point_mass_gravity_type, "Delfi-C3", "Mars"
        ),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.point_mass_gravity_type, "Delfi-C3", "Venus"
        ),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.spherical_harmonic_gravity_type, "Delfi-C3", "Earth"
        ),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.aerodynamic_type, "Delfi-C3", "Earth"
        ),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.cannonball_radiation_pressure_type, "Delfi-C3", "Sun"
        )
    ]

    # Create propagation settings.
    termination_condition = propagation_setup.propagator.time_termination(simulation_end_epoch)
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        initial_state,
        termination_condition,
        output_variables=dependent_variables_to_save
    )
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
    dynamics_simulator = numerical_simulation.SingleArcSimulator(
        bodies, integrator_settings, propagator_settings)
    states = dynamics_simulator.state_history
    dependent_variables = dynamics_simulator.dependent_variable_history

    ###########################################################################
    # PLOT RESULTS    #########################################################
    ###########################################################################

    time = dependent_variables.keys()
    time_hours = [t / 3600 for t in time]

    dependent_variable_list = np.vstack(list(dependent_variables.values()))

    # Plot total acceleration as function of time
    total_acceleration_norm = np.linalg.norm(dependent_variable_list[:, 0:3], axis=1)
    plt.figure(figsize=(17, 5))
    plt.plot(time_hours, total_acceleration_norm)
    plt.xlabel('Time [hr]')
    plt.ylabel('Total Acceleration [m/s$^2$]')
    plt.xlim([min(time_hours), max(time_hours)])
    plt.grid()

    # Plot ground track for a period of 3 hours
    latitude = dependent_variable_list[:, 9]
    longitude = dependent_variable_list[:, 10]
    hours = 3
    subset = int(len(time) / 24 * hours)
    latitude = np.rad2deg(latitude[0: subset])
    longitude = np.rad2deg(longitude[0: subset])
    plt.figure(figsize=(17, 5))
    plt.scatter(longitude, latitude, s=1)
    plt.xlabel('Longitude [deg]')
    plt.ylabel('Latitude [deg]')
    plt.xlim([min(longitude), max(longitude)])
    plt.yticks(np.arange(-90, 91, step=45))
    plt.grid()

    # Plot Kepler elements as a function of time
    kepler_elements = dependent_variable_list[:, 3:9]
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(20, 17))

    # Semi-major Axis
    semi_major_axis = [element / 1000 for element in kepler_elements[:, 0]]
    ax1.plot(time_hours, semi_major_axis)
    ax1.set_ylabel('Semi-major axis [km]')

    # Eccentricity
    eccentricity = kepler_elements[:, 1]
    ax2.plot(time_hours, eccentricity)
    ax2.set_ylabel('Eccentricity [-]')

    # Inclination
    inclination = [np.rad2deg(element) for element in kepler_elements[:, 2]]
    ax3.plot(time_hours, inclination)
    ax3.set_ylabel('Inclination [deg]')

    # Argument of Periapsis
    argument_of_periapsis = [np.rad2deg(element) for element in kepler_elements[:, 3]]
    ax4.plot(time_hours, argument_of_periapsis)
    ax4.set_ylabel('Argument of Periapsis [deg]')

    # Right Ascension of the Ascending Node
    raan = [np.rad2deg(element) for element in kepler_elements[:, 4]]
    ax5.plot(time_hours, raan)
    ax5.set_ylabel('RAAN [deg]')

    # True Anomaly
    true_anomaly = [np.rad2deg(element) for element in kepler_elements[:, 5]]
    ax6.scatter(time_hours, true_anomaly, s=1)
    ax6.set_ylabel('True Anomaly [deg]')
    ax6.set_yticks(np.arange(0, 361, step=60))

    for ax in fig.get_axes():
        ax.set_xlabel('Time [hr]')
        ax.set_xlim([min(time_hours), max(time_hours)])
        ax.grid()
    plt.figure(figsize=(17, 5))

    # Plot accelerations as a function of time

    # Point Mass Gravity Acceleration Sun
    acceleration_norm_pm_sun = dependent_variable_list[:, 11]
    plt.plot(time_hours, acceleration_norm_pm_sun, label='PM Sun')

    # Point Mass Gravity Acceleration Moon
    acceleration_norm_pm_moon = dependent_variable_list[:, 12]
    plt.plot(time_hours, acceleration_norm_pm_moon, label='PM Moon')

    # Point Mass Gravity Acceleration Mars
    acceleration_norm_pm_mars = dependent_variable_list[:, 13]
    plt.plot(time_hours, acceleration_norm_pm_mars, label='PM Mars')

    # Point Mass Gravity Acceleration Venus
    acceleration_norm_pm_venus = dependent_variable_list[:, 14]
    plt.plot(time_hours, acceleration_norm_pm_venus, label='PM Venus')

    # Spherical Harmonic Gravity Acceleration Earth
    acceleration_norm_sh_earth = dependent_variable_list[:, 15]
    plt.plot(time_hours, acceleration_norm_sh_earth, label='SH Earth')

    # Aerodynamic Acceleration Earth
    acceleration_norm_aero_earth = dependent_variable_list[:, 16]
    plt.plot(time_hours, acceleration_norm_aero_earth, label='Aerodynamic Earth')

    # Cannonball Radiation Pressure Acceleration Sun
    acceleration_norm_rp_sun = dependent_variable_list[:, 17]
    plt.plot(time_hours, acceleration_norm_rp_sun, label='Radiation Pressure Sun')

    plt.grid()

    plt.xlim([min(time_hours), max(time_hours)])
    plt.xlabel('Time [hr]')
    plt.ylabel('Acceleration Norm [m/s$^2$]')

    plt.legend(bbox_to_anchor=(1.04, 1))
    plt.yscale('log')

    plt.show()

    # Final statement (not required, though good practice in a __main__).
    return 0


if __name__ == "__main__":
    main()
