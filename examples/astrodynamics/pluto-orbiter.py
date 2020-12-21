"""Example of initial interplanetary Earth-Jupiter-Pluto Patched Conics

This script demonstrates the preliminary mission design for a
Earth-Jupiter-Pluto transfer. The launch date, sequence of gravity assists
and duration of interplanetary were optimised using a simple custom genetic
algorithm subject to the following systems engineering requirements:

- R-MIS-071 The duration of the mission shall be no more than 25 years.
- R-MIS-010 The exploration mission shall be an orbital exploration of the
    Pluto-Charon system.
- R-MIS-020 The vehicle shall orbit the Pluto-Charon system for a minimum
    of one year.

The analysis assumes impulsive thrust, therefore velocity changes occur
instantly. The optimisation modelled the trajectory as a sequence of
conic sections, commonly known as the patched conic approach.

Required files
==============
- de430.bsp
- de431_part-1.bsp
- de431_part-2.bsp

Required packages
=================
- numpy [conda/pip]
- tudatpy >= 0.8.21 [conda]
- python-decouple [pip]

Reference
=========
DSE Group 13 2018: AE3200 Design Synthesis Exercise
"""

###############################################################################
# IMPORT STATEMENTS ###########################################################
###############################################################################
import numpy as np
from decouple import config
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.astro.fundamentals import propagate_kepler_orbit
from tudatpy.kernel.astro.conversion import cartesian_to_keplerian
from tudatpy.kernel.astro.conversion import keplerian_to_cartesian
from tudatpy.kernel.astro.two_body_dynamics import (
    calculate_departure_parameters,
    calculate_arrival_parameters,
    calculate_gravity_assist_parameters,
    LambertTargeterIzzo)

if __name__ == "__main__":

    # retrieve spice kernel paths from .env file
    SPICE_KERNEL_BSP_DE430 = config("SPICE_KERNEL_BSP_DE430")
    SPICE_KERNEL_BSP_DE431_PART1 = config("SPICE_KERNEL_BSP_DE431_PART1")
    SPICE_KERNEL_BSP_DE431_PART2 = config("SPICE_KERNEL_BSP_DE431_PART2")

    # append the additional spice kernels to the standard ones
    spice_kernels = spice_interface.get_standard_kernels() + [
        SPICE_KERNEL_BSP_DE430,
        SPICE_KERNEL_BSP_DE431_PART1,
        SPICE_KERNEL_BSP_DE431_PART2]

    # load all desired spice kernels
    for kernel in spice_kernels:
        spice_interface.load_kernel(kernel)

    # define reference frame
    frame_origin = "SSB"
    frame_orientation = "ECLIPJ2000"

    # mission system parent body
    parent_body = "Sun"
    parent_body_mu = spice_interface.get_body_gravitational_parameter("Sun")

    # mission planetary departure parameter definition
    departure_date = "November 27, 2027"
    departure_body = "Earth"
    departure_body_system = "EARTH BARYCENTER"
    departure_r_limit = 1.048

    # mission gravity assist parameter definition
    gravity_assist_date = "December 19, 2029"
    gravity_assist_body = "Jupiter"
    gravity_assist_body_system = "JUPITER BARYCENTER"
    gravity_assist_r_limit = 1.60

    # mission planetary arrival parameter definition
    arrival_date = "12 November, 2051"
    arrival_body = "Pluto"
    arrival_body_system = "PLUTO BARYCENTER"
    arrival_r_limit = 1.076

    # calculate ephemeris time and time of flight parameters
    departure_ephemeris_time = spice_interface.convert_date_string_to_ephemeris_time(departure_date)
    arrival_ephemeris_time = spice_interface.convert_date_string_to_ephemeris_time(arrival_date)
    gravity_assist_ephemeris_time = spice_interface.convert_date_string_to_ephemeris_time(gravity_assist_date)

    # determine departure geometry
    departure_body_mu = spice_interface.get_body_gravitational_parameter(departure_body)
    departure_body_radius = spice_interface.get_average_radius(departure_body)
    departure_body_state = spice_interface.get_body_cartesian_state_at_epoch(
        target_body_name=departure_body,
        observer_body_name=frame_origin,
        reference_frame_name=frame_orientation,
        aberration_corrections="none",
        ephemeris_time=departure_ephemeris_time)
    departure_body_velocity = departure_body_state[3:]
    departure_body_position = departure_body_state[:3]

    # determine gravity assist geometry
    gravity_assist_body_mu = spice_interface.get_body_gravitational_parameter(gravity_assist_body)
    gravity_assist_body_radius = spice_interface.get_average_radius(gravity_assist_body)
    gravity_assist_body_state = spice_interface.get_body_cartesian_state_at_epoch(
        target_body_name=gravity_assist_body_system,
        observer_body_name=frame_origin,
        reference_frame_name=frame_orientation,
        aberration_corrections="none",
        ephemeris_time=gravity_assist_ephemeris_time)
    gravity_assist_body_velocity = gravity_assist_body_state[3:]
    gravity_assist_body_position = gravity_assist_body_state[:3]

    # determine arrival geometry
    arrival_body_mu = spice_interface.get_body_gravitational_parameter(arrival_body)
    arrival_body_radius = spice_interface.get_average_radius(arrival_body)
    arrival_body_state = spice_interface.get_body_cartesian_state_at_epoch(
        target_body_name=arrival_body_system,
        observer_body_name=frame_origin,
        reference_frame_name=frame_orientation,
        aberration_corrections="none",
        ephemeris_time=arrival_ephemeris_time)
    arrival_body_velocity = arrival_body_state[3:]
    arrival_body_position = arrival_body_state[:3]

    # calculate first interplanetary coast leg
    departure_1_velocity, arrival_1_velocity = LambertTargeterIzzo(
        departure_position=departure_body_position,
        arrival_position=gravity_assist_body_position,
        time_of_flight=gravity_assist_ephemeris_time - departure_ephemeris_time,
        gravitational_parameter=parent_body_mu
    ).get_velocity_vectors()

    # calculate second interplanetary coast leg
    departure_2_velocity, arrival_2_velocity = LambertTargeterIzzo(
        departure_position=gravity_assist_body_position,
        arrival_position=arrival_body_position,
        time_of_flight=arrival_ephemeris_time - gravity_assist_ephemeris_time,
        gravitational_parameter=parent_body_mu
    ).get_velocity_vectors()

    # calculate departure parameters
    departure_parameters = calculate_departure_parameters(
        central_body_gravitational_parameter=departure_body_mu,
        central_body_velocity_on_exit_soi=departure_body_velocity,
        outgoing_velocity=departure_1_velocity,
        departure_periapsis_distance=departure_body_radius * departure_r_limit)

    # calculate gravity assist parameters
    gravity_assist_parameters = calculate_gravity_assist_parameters(
        central_body_gravitational_parameter=departure_body_mu,
        central_body_velocity_on_entry_soi=gravity_assist_body_velocity,
        central_body_velocity_on_exit_soi=gravity_assist_body_velocity,
        incoming_velocity=arrival_1_velocity,
        outgoing_velocity=departure_2_velocity,
        smallest_periapsis_distance=gravity_assist_body_radius * gravity_assist_r_limit,
        speed_tolerance=1e-6)

    # calculate arrival parameters
    arrival_parameters = calculate_arrival_parameters(
        central_body_gravitational_parameter=arrival_body_mu,
        central_body_velocity_on_entry_soi=arrival_body_velocity,
        incoming_velocity=arrival_2_velocity,
        arrival_periapsis_distance=arrival_body_radius * arrival_r_limit)

    # calculate samples of first interplanetary leg for plotting
    departure_1_inertial_state = np.hstack([departure_body_position,
                                            departure_1_velocity])
    trajectory_1_initial_keplerian = cartesian_to_keplerian(
        cartesian_elements=departure_1_inertial_state,
        gravitational_parameter=parent_body_mu)
    trajectory_1_cartesian = np.vstack([keplerian_to_cartesian(
        keplerian_elements=propagate_kepler_orbit(
            initial_keplerian_state=trajectory_1_initial_keplerian,
            gravitational_parameter=parent_body_mu,
            propagation_time=t),
        gravitational_parameter=parent_body_mu) for t in
        np.linspace(0, gravity_assist_ephemeris_time - departure_ephemeris_time, 100)])

    # calculate samples of second interplanetary leg for plotting
    departure_2_inertial_state = np.hstack([gravity_assist_body_position,
                                            departure_2_velocity])
    trajectory_2_initial_keplerian = cartesian_to_keplerian(
        cartesian_elements=departure_2_inertial_state,
        gravitational_parameter=parent_body_mu)
    trajectory_2_cartesian = np.vstack([keplerian_to_cartesian(
        keplerian_elements=propagate_kepler_orbit(
            initial_keplerian_state=trajectory_2_initial_keplerian,
            gravitational_parameter=parent_body_mu,
            propagation_time=t),
        gravitational_parameter=parent_body_mu) for t in
        np.linspace(0, arrival_ephemeris_time - gravity_assist_ephemeris_time, 100)])

    spice_kernel_repr = "\n".join(spice_kernels)
    print(departure_parameters)
    mission_details = f"""
Spice Kernels Loaded
====================
{spice_kernel_repr}

{departure_body} Departure: {departure_date}
{"=" * len(f"{departure_body} Departure: {departure_date}")}
departure_parameters.radius_of_periapsis: {departure_parameters.radius_of_periapsis}
departure_parameters.initial_velocity_at_periapsis: {departure_parameters.initial_velocity_at_periapsis}
departure_parameters.outgoing_semi_major_axis:  {departure_parameters.outgoing_semi_major_axis}
departure_parameters.outgoing_eccentricity:  {departure_parameters.outgoing_eccentricity}
departure_parameters.final_velocity_at_periapsis:  {departure_parameters.final_velocity_at_periapsis}
departure_parameters.outgoing_hyperbolic_excess_velocity: {departure_parameters.outgoing_hyperbolic_excess_velocity}
departure_parameters.total_delta_v: {departure_parameters.total_delta_v}

{gravity_assist_body} Gravity Assist: {gravity_assist_date}
{"=" * len(f"{gravity_assist_body} Gravity Assist: {gravity_assist_date}")}
gravity_assist_parameters.radius_of_periapsis: {gravity_assist_parameters.radius_of_periapsis}
gravity_assist_parameters.incoming_semi_major_axis: {gravity_assist_parameters.incoming_semi_major_axis}
gravity_assist_parameters.incoming_eccentricity: {gravity_assist_parameters.incoming_eccentricity}
gravity_assist_parameters.incoming_velocity_at_periapsis: {gravity_assist_parameters.incoming_velocity_at_periapsis}
gravity_assist_parameters.outgoing_semi_major_axis: {gravity_assist_parameters.outgoing_semi_major_axis}
gravity_assist_parameters.outgoing_eccentricity: {gravity_assist_parameters.outgoing_eccentricity}
gravity_assist_parameters.outgoing_velocity_at_periapsis: {gravity_assist_parameters.outgoing_velocity_at_periapsis}
gravity_assist_parameters.total_delta_v: {gravity_assist_parameters.total_delta_v}
gravity_assist_parameters.incoming_hyperbolic_excess_velocity: {gravity_assist_parameters.incoming_hyperbolic_excess_velocity}
gravity_assist_parameters.outgoing_hyperbolic_excess_velocity: {gravity_assist_parameters.outgoing_hyperbolic_excess_velocity}

{arrival_body} Arrival: {arrival_date}
{"=" * len(f"{arrival_body} Arrival: {arrival_date}")}
arrival_parameters.radius_of_periapsis: {arrival_parameters.radius_of_periapsis}
arrival_parameters.initial_velocity_at_periapsis: {arrival_parameters.initial_velocity_at_periapsis}
arrival_parameters.incoming_semi_major_axis: {arrival_parameters.incoming_semi_major_axis}
arrival_parameters.incoming_eccentricity: {arrival_parameters.incoming_eccentricity}
arrival_parameters.final_velocity_at_periapsis: {arrival_parameters.final_velocity_at_periapsis}
arrival_parameters.incoming_hyperbolic_excess_velocity: {arrival_parameters.incoming_hyperbolic_excess_velocity}
arrival_parameters.total_delta_v: {arrival_parameters.total_delta_v}
"""
    print(mission_details)

    import matplotlib.pyplot as plt
    import numpy as np

    scale = np.linalg.norm(departure_body_position[:2])
    sun_circle = plt.Circle((0, 0), 0.6, color='yellow')
    earth_circle = plt.Circle(departure_body_position[:2] / scale, 0.4, color='green')
    earth_orbit_circle = plt.Circle(
        (0, 0),
        np.linalg.norm(departure_body_position[:2]) / scale,
        color='green')

    jupiter_circle = plt.Circle(gravity_assist_body_position[:2] / scale, 0.4, color='red')
    jupiter_orbit_circle = plt.Circle(
        (0, 0),
        np.linalg.norm(gravity_assist_body_position[:2]) / scale,
        color='orange')

    pluto_circle = plt.Circle(arrival_body_position[:2] / scale, 0.4, color='pink')
    pluto_orbit_circle = plt.Circle(
        (0, 0),
        np.linalg.norm(arrival_body_position[:2] / scale),
        color='pink')

    fig, ax = plt.subplots()  # note we must use plt.subplots, not plt.subplot
    # (or if you have an existing figure)
    # fig = plt.gcf()
    # ax = fig.gca()

    # ax.add_artist(pluto_orbit_circle)
    # ax.add_artist(jupiter_orbit_circle)
    # ax.add_artist(earth_orbit_circle)
    #
    ax.add_artist(sun_circle)
    ax.add_artist(earth_circle)
    ax.add_artist(jupiter_circle)
    ax.add_artist(pluto_circle)

    ax.set_aspect(aspect="equal")


    def plot_2d_trajectory(trajectory):
        ax.plot(trajectory[:, 0] / scale, trajectory[:, 1] / scale)


    ax.plot(trajectory_1_cartesian[:, 0] / scale, trajectory_1_cartesian[:, 1] / scale)
    ax.plot(trajectory_2_cartesian[:, 0] / scale, trajectory_2_cartesian[:, 1] / scale)

    prior_ylim = ax.get_ylim()
    prior_xlim = ax.get_xlim()


    # Plot planet orbits
    def plot_circular(ax, fig, r):
        a = np.linspace(0, np.pi * 2, 200)
        ax.plot(r * np.cos(a), r * np.sin(a))


    plot_circular(ax, fig, np.linalg.norm(departure_body_position) / scale)
    plot_circular(ax, fig, np.linalg.norm(gravity_assist_body_position) / scale)
    plot_circular(ax, fig, np.linalg.norm(arrival_body_position) / scale)

    plt.ylim(prior_ylim)
    plt.xlim(prior_xlim)

    # plt.ylim((-np.abs(arrival_body_position[1] / scale),
    #           +np.abs(arrival_body_position[1] / scale)))
    #
    # plt.xlim((-np.abs(arrival_body_position[0] / scale),
    #           +np.abs(arrival_body_position[0] / scale)))

    plt.show()
    # fig.savefig('plotcircles.png')
