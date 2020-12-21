"""
Case Study: Mars 2020
---------------------
Mars 2020 is a Mars rover mission by NASA's Mars Exploration Program that
includes the Perseverance rover and the Ingenuity helicopter drone. It was
launched on 30 July 2020 at 11:50 UTC, and will touch down in Jezero
crater on Mars on 18 February 2021.
"""

###############################################################################
# IMPORT STATEMENTS ###########################################################
###############################################################################
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.astro.fundamentals import propagate_kepler_orbit
from tudatpy.kernel.astro.two_body_dynamics import (
    calculate_departure_parameters,
    calculate_arrival_parameters,
    calculate_gravity_assist_parameters,
    LambertTargeterIzzo,
    PlanetaryDeparture,
    PlanetaryRendezvous,
    KeplerianState,
    CartesianState
)

import numpy as np
import matplotlib.pyplot as plt
from tudatpy.kernel.astro.conversion import keplerian_to_cartesian


def print_rendezvous(test_planetary_departure):
    print("v_inf: ", test_planetary_departure.hyperbolic_excess_velocity)
    print("v_inf: ", test_planetary_departure.hyperbolic_excess_speed)
    print("ecc: ", test_planetary_departure.eccentricity)
    print("h: ", test_planetary_departure.angular_momentum)
    print("v_p: ", test_planetary_departure.periapsis_speed)
    print("v_c: ", test_planetary_departure.parking_orbit_speed)
    print("dv: ", test_planetary_departure.required_delta_v)
    print("beta_angle", np.rad2deg(test_planetary_departure.beta_angle))
    print("turning_angle", np.rad2deg(test_planetary_departure.turning_angle))
    print("orbital_plane_vector: ", test_planetary_departure.orbital_plane_vector)
    print("periapsis_velocity_unit_vector: ", test_planetary_departure.periapsis_velocity_unit_vector)
    print("outgoing: ", test_planetary_departure.periapsis_outgoing_velocity)
    print("incoming: ", test_planetary_departure.periapsis_incoming_velocity)
    print(test_planetary_departure.periapsis_position)
    print(test_planetary_departure.sphere_of_influence_distance_at_entrance)
    print(test_planetary_departure.semi_major_axis_incoming)
    print(test_planetary_departure.true_anomaly_entrance)
    print(test_planetary_departure.hyperbolic_anomaly_at_entrance)
    print(test_planetary_departure.time_of_flight_from_entrance)


def print_departure(test_planetary_departure):
    print("v_inf: ", test_planetary_departure.hyperbolic_excess_velocity)
    print("v_inf: ", test_planetary_departure.hyperbolic_excess_speed)
    print("ecc: ", test_planetary_departure.eccentricity)
    print("h: ", test_planetary_departure.angular_momentum)
    print("v_p: ", test_planetary_departure.periapsis_speed)
    print("v_c: ", test_planetary_departure.parking_orbit_speed)
    print("dv: ", test_planetary_departure.required_delta_v)
    print("beta_angle", np.rad2deg(test_planetary_departure.beta_angle))
    print("orbital_plane_vector: ", test_planetary_departure.orbital_plane_vector)
    print("periapsis_velocity_unit_vector: ", test_planetary_departure.periapsis_velocity_unit_vector)
    print("outgoing: ", test_planetary_departure.periapsis_outgoing_velocity)
    print("incoming: ", test_planetary_departure.periapsis_incoming_velocity)
    print(test_planetary_departure.periapsis_position)
    print(test_planetary_departure.sphere_of_influence_distance_at_exit)
    print(test_planetary_departure.semi_major_axis_outgoing)
    print(test_planetary_departure.true_anomaly_exit)
    print(test_planetary_departure.hyperbolic_anomaly_at_exit)
    print("tof: ", test_planetary_departure.time_of_flight_to_exit / 3600 / 24)


def plot_planetary_departure(
        planetary_departure: PlanetaryDeparture,
        simulated_trajectory: np.ndarray = None,

):
    outgoing_cartesian_state = planetary_departure.outgoing_cartesian_state.to_cartesian().state_vector
    incoming_cartesian_state = planetary_departure.outgoing_cartesian_state.to_cartesian().state_vector
    outgoing_keplerian_state = planetary_departure.outgoing_cartesian_state.to_keplerian().state_vector
    incoming_keplerian_state = planetary_departure.incoming_cartesian_state.to_keplerian().state_vector
    # incoming_mean_motion = planetary_departure.incoming_cartesian_state.mean_motion
    departure_body_mu = planetary_departure.outgoing_cartesian_state.gravitational_parameter
    radius_soi = planetary_departure.sphere_of_influence_distance_at_exit

    # sample points for outgoing hyperbolic trajectory
    time_array_outgoing = np.linspace(0, planetary_departure.time_of_flight_to_exit, 300, endpoint=True)
    trajectory_2_cartesian = np.vstack([keplerian_to_cartesian(
        keplerian_elements=propagate_kepler_orbit(
            initial_keplerian_state=outgoing_keplerian_state,
            gravitational_parameter=departure_body_mu,
            propagation_time=t),
        gravitational_parameter=departure_body_mu) for t in
        time_array_outgoing])

    # sample points for incoming circular parking orbit
    trajectory_1_cartesian = np.vstack([keplerian_to_cartesian(
        keplerian_elements=propagate_kepler_orbit(
            initial_keplerian_state=incoming_keplerian_state,
            gravitational_parameter=departure_body_mu,
            propagation_time=t),
        gravitational_parameter=departure_body_mu) for t in
        np.linspace(-3600, 0, 300)])

    # scaling
    scale = 1 / radius_soi

    # create body
    earth_circle = plt.Circle((0, 0), 6700 * 1e3 * scale, color='green')
    fig, ax = plt.subplots()  # note we must use plt.subplots, not plt.subplot
    ax.add_artist(earth_circle)
    ax.set_aspect(aspect="equal")

    # plot body velocity vector
    vp = planetary_departure.central_body_state.velocity_vector
    vp = vp / np.linalg.norm(vp)
    ax.arrow(0, 0, vp[0], vp[1], head_width=0.05, head_length=0.1, fc='k', ec='k')

    # plot parent body unit vector
    r_sun = planetary_departure.central_body_state.position_vector
    r_sun = - r_sun / np.linalg.norm(r_sun)
    ax.arrow(0, 0, r_sun[0], r_sun[1], head_width=0.05, head_length=0.1, fc='k', ec='k')

    # plot excess velocity vector
    v_out = planetary_departure.hyperbolic_excess_velocity
    v_out = v_out / np.linalg.norm(v_out)
    ax.arrow(0, 0, v_out[0], v_out[1], head_width=0.05, head_length=0.1, fc='k', ec='k')

    # sample points for plotting sphere of influence
    theta = np.linspace(0, np.pi * 2, 300)
    x, y = radius_soi * np.cos(theta) * scale, radius_soi * np.sin(theta) * scale

    # plot sphere of influence
    ax.plot(x, y, linestyle=":")

    ax.plot(trajectory_1_cartesian[:, 0] * scale, trajectory_1_cartesian[:, 1] * scale)
    ax.plot(trajectory_2_cartesian[:, 0] * scale, trajectory_2_cartesian[:, 1] * scale)

    plt.show()


def plot_planetary_rendezvous(
        planetary_rendezvous: PlanetaryRendezvous,
        simulated_trajectory: np.ndarray = None,

):
    outgoing_cartesian_state = planetary_rendezvous.outgoing_cartesian_state.to_cartesian().state_vector
    incoming_cartesian_state = planetary_rendezvous.outgoing_cartesian_state.to_cartesian().state_vector
    outgoing_keplerian_state = planetary_rendezvous.outgoing_cartesian_state.to_keplerian().state_vector
    incoming_keplerian_state = planetary_rendezvous.incoming_cartesian_state.to_keplerian().state_vector
    # incoming_mean_motion = planetary_departure.incoming_cartesian_state.mean_motion
    departure_body_mu = planetary_rendezvous.outgoing_cartesian_state.gravitational_parameter
    radius_soi = planetary_rendezvous.sphere_of_influence_distance_at_entrance
    tof = planetary_rendezvous.time_of_flight_from_entrance
    print("ROF: ", tof)

    # sample points for outgoing hyperbolic trajectory
    time_array_outgoing = np.linspace(tof, 0, 700, endpoint=True)
    trajectory_2_cartesian = np.vstack([keplerian_to_cartesian(
        keplerian_elements=propagate_kepler_orbit(
            initial_keplerian_state=incoming_keplerian_state,
            gravitational_parameter=departure_body_mu,
            propagation_time=t),
        gravitational_parameter=departure_body_mu) for t in
        time_array_outgoing])

    # sample points for incoming circular parking orbit
    trajectory_1_cartesian = np.vstack([keplerian_to_cartesian(
        keplerian_elements=propagate_kepler_orbit(
            initial_keplerian_state=outgoing_keplerian_state,
            gravitational_parameter=departure_body_mu,
            propagation_time=t),
        gravitational_parameter=departure_body_mu) for t in
        np.linspace(0, 3600, 100)])

    # scaling
    scale = 1 / radius_soi

    # create body
    earth_circle = plt.Circle((0, 0), 6700 * 1e3 * scale, color='green')
    fig, ax = plt.subplots()  # note we must use plt.subplots, not plt.subplot
    ax.add_artist(earth_circle)
    ax.set_aspect(aspect="equal")

    # plot body velocity vector
    vp = planetary_departure.central_body_state.velocity_vector
    vp = vp / np.linalg.norm(vp)
    ax.arrow(0, 0, vp[0], vp[1], head_width=0.05, head_length=0.1, fc='k', ec='k')

    # plot parent body unit vector
    r_sun = planetary_departure.central_body_state.position_vector
    r_sun = - r_sun / np.linalg.norm(r_sun)
    ax.arrow(0, 0, r_sun[0], r_sun[1], head_width=0.05, head_length=0.1, fc='k', ec='k')

    # plot excess velocity vector
    v_out = planetary_departure.hyperbolic_excess_velocity
    v_out = v_out / np.linalg.norm(v_out)
    ax.arrow(-v_out[0], -v_out[1], v_out[0], v_out[1], head_width=0.05, head_length=0.1, fc='k', ec='k')

    # sample points for plotting sphere of influence
    theta = np.linspace(0, np.pi * 2, 300)
    x, y = radius_soi * np.cos(theta) * scale, radius_soi * np.sin(theta) * scale

    # plot sphere of influence
    ax.plot(x, y, linestyle=":")

    ax.plot(trajectory_1_cartesian[:, 0] * scale, trajectory_1_cartesian[:, 1] * scale)
    ax.plot(trajectory_2_cartesian[:, 0] * scale, trajectory_2_cartesian[:, 1] * scale)

    plt.show()


if __name__ == "__main__":

    spice_kernels = spice_interface.get_standard_kernels() + [
        "/home/ggarrett/Downloads/de430.bsp",
        "/home/ggarrett/Downloads/de431_part-1.bsp",
        "/home/ggarrett/Downloads/de431_part-2.bsp"]

    for kernel in spice_kernels:
        spice_interface.load_kernel(kernel)

    # define reference frame
    frame_origin = "SSB"
    frame_orientation = "ECLIPJ2000"

    # mission system parent body
    parent_body = "Sun"
    parent_body_mu = spice_interface.get_body_gravitational_parameter(parent_body)

    # Mission planetary departure parameter definition
    departure_date = "July 30, 2020"
    departure_body = "Earth"
    departure_r_limit = 1.048

    # Mission planetary arrival parameter definition
    rendezvous_date = "18 February, 2021"
    rendezvous_body = "Mars"
    rendezvous_r_limit = 1.05

    # Calculate ephemeris time and time of flight parameters
    departure_reference_ephemeris_time = spice_interface.convert_date_string_to_ephemeris_time(departure_date)
    rendezvous_reference_ephemeris_time = spice_interface.convert_date_string_to_ephemeris_time(rendezvous_date)

    # Determine constant departure parameters
    departure_body_mu = spice_interface.get_body_gravitational_parameter(departure_body)
    departure_body_radius = spice_interface.get_average_radius(departure_body)

    # Determine constant arrival geometry
    rendezvous_body_mu = spice_interface.get_body_gravitational_parameter(rendezvous_body)
    rendezvous_body_radius = spice_interface.get_average_radius(rendezvous_body)

    # Define departure body state at departure conic reference
    departure_body_state_at_reference = CartesianState.from_spice(
        body_name=departure_body,
        ephemeris_time=departure_reference_ephemeris_time,
        frame_origin=frame_origin,
        frame_orientation=frame_orientation,
        aberration_corrections="none"
    )

    # Determine rendezvous body state at rendezvous conic reference
    rendezvous_body_state_at_reference = CartesianState.from_spice(
        body_name=rendezvous_body,
        ephemeris_time=rendezvous_reference_ephemeris_time,
        frame_origin=frame_origin,
        frame_orientation=frame_orientation,
        aberration_corrections="none")

    # Set initial tof for departure and arrival conic
    tof_departure_conic = 0
    tof_rendezvous_conic = 0

    # Set initial planetocentric conic exit positions
    coast_conic_entry_position = departure_body_state_at_reference.position_vector
    coast_conic_exit_position = rendezvous_body_state_at_reference.position_vector

    # Define tolerance for patched conics analytical solution at conic boundaries
    velocity_error = velocity_tolerance = 1e-3
    position_error = position_tolerance = 1e-3

    while ((velocity_error >= velocity_tolerance) |
           (position_error >= position_tolerance)):
        # Calculate interplanetary coast leg
        rendezvous_conic_entry_time = rendezvous_reference_ephemeris_time + tof_rendezvous_conic  # returns negative, hence +
        departure_conic_exit_time = departure_reference_ephemeris_time + tof_departure_conic
        interplanetary_coast_tof = rendezvous_conic_entry_time - departure_conic_exit_time

        coast_conic_entry_velocity, coast_conic_exit_velocity = LambertTargeterIzzo(
            departure_position=coast_conic_entry_position,
            arrival_position=coast_conic_exit_position,
            time_of_flight=interplanetary_coast_tof,
            gravitational_parameter=parent_body_mu
        ).get_velocity_vectors()

        # Define departure body state at conic exit state
        departure_body_state_at_exit = CartesianState.from_spice(
            body_name=departure_body,
            ephemeris_time=departure_conic_exit_time,
            frame_origin=frame_origin,
            frame_orientation=frame_orientation,
            aberration_corrections="none")

        # Calculate departure parameters
        planetary_departure = PlanetaryDeparture(
            outgoing_velocity=coast_conic_entry_velocity,
            departure_body_state_at_reference=departure_body_state_at_reference,
            departure_body_state_at_exit=departure_body_state_at_exit,
            periapsis_distance=departure_body_radius * departure_r_limit,
            gravitational_parameter=departure_body_mu,
            prograde_orbit=True)

        # Determine reference rendezvous parameters
        rendezvous_body_state_at_entry = CartesianState.from_spice(
            body_name=rendezvous_body,
            ephemeris_time=rendezvous_conic_entry_time,
            frame_origin=frame_origin,
            frame_orientation=frame_orientation,
            aberration_corrections="none")

        # Calculate arrival parameters
        planetary_rendezvous = PlanetaryRendezvous(
            incoming_velocity=coast_conic_exit_velocity,
            rendezvous_body_state_at_reference=rendezvous_body_state_at_reference,
            rendezvous_body_state_at_entry=rendezvous_body_state_at_entry,
            periapsis_distance=rendezvous_r_limit * rendezvous_body_radius,
            gravitational_parameter=rendezvous_body_mu,
            prograde_orbit=True
        )

        # Retrieve departure and rendezvous conic exit & entry states
        tof_departure_conic = planetary_departure.time_of_flight_to_exit
        tof_rendezvous_conic = planetary_rendezvous.time_of_flight_from_entrance
        departure_conic_exit_state = planetary_departure.outgoing_cartesian_state.to_keplerian().propagate(
            tof_departure_conic).to_cartesian()
        rendezvous_conic_entry_state = planetary_rendezvous.incoming_cartesian_state.to_keplerian().propagate(
            tof_rendezvous_conic).to_cartesian()
        updated_coast_conic_entry_velocity = (departure_conic_exit_state.velocity_vector
                                              + departure_body_state_at_exit.velocity_vector)
        updated_coast_conic_entry_position = (departure_conic_exit_state.position_vector
                                              + departure_body_state_at_exit.position_vector)
        updated_coast_conic_exit_velocity = (rendezvous_conic_entry_state.velocity_vector
                                             + rendezvous_body_state_at_entry.velocity_vector)
        updated_coast_conic_exit_position = (rendezvous_conic_entry_state.position_vector
                                             + rendezvous_body_state_at_entry.position_vector)

        velocity_error = (
                np.linalg.norm(updated_coast_conic_exit_velocity - coast_conic_exit_velocity)
                + np.linalg.norm(updated_coast_conic_entry_velocity - coast_conic_entry_velocity)
        )

        position_error = (
                np.linalg.norm(updated_coast_conic_exit_position - coast_conic_exit_position)
                + np.linalg.norm(updated_coast_conic_entry_position - coast_conic_entry_position)
        )

        print("tof_exit_conic: ", tof_departure_conic, " [s]")
        print("tof_entry_conic: ", tof_rendezvous_conic, " [s]")
        print("position_error: ", position_error, " [m]")
        print("velocity_error: ", velocity_error, " [m/s]")

        if ((velocity_error >= velocity_tolerance) |
                (position_error >= position_tolerance)):
            # Update rendezvous and departure conic exit positions
            coast_conic_exit_position = updated_coast_conic_exit_position
            coast_conic_entry_position = updated_coast_conic_entry_position
            coast_conic_exit_velocity = updated_coast_conic_exit_velocity
            coast_conic_entry_velocity = updated_coast_conic_entry_velocity

    plot_planetary_departure(planetary_departure)
    plot_planetary_rendezvous(planetary_rendezvous)
    print_departure(planetary_departure)
    print_rendezvous(planetary_rendezvous)

    spice_kernel_repr = "\n".join(spice_kernels)
    print(
        f"""
Spice Kernels Loaded
====================
{spice_kernel_repr}

""")
