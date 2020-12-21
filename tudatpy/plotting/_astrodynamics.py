import numpy as np
from ..kernel.astro.two_body_dynamics.iod import (
    PlanetaryDeparture,
    PlanetaryRendezvous
)
from ..kernel.astro.fundamentals import propagate_kepler_orbit
from ..kernel.astro.conversion import keplerian_to_cartesian
import matplotlib.pyplot as plt


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
    parking_orbit_period = 1/planetary_departure.incoming_cartesian_state.to_keplerian().mean_motion
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
        np.linspace(-parking_orbit_period, 0, 300)])

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
    if simulated_trajectory:
        ax.plot(simulated_trajectory[:, 0] * scale, simulated_trajectory[:, 1] * scale)

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
    parking_orbit_period = 1/planetary_rendezvous.outgoing_cartesian_state.to_keplerian().mean_motion
    radius_soi = planetary_rendezvous.sphere_of_influence_distance_at_entrance
    tof = planetary_rendezvous.time_of_flight_from_entrance

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
        np.linspace(0, parking_orbit_period, 100)])

    # scaling
    scale = 1 / radius_soi

    # create body
    earth_circle = plt.Circle((0, 0), 6700 * 1e3 * scale, color='green')
    fig, ax = plt.subplots()  # note we must use plt.subplots, not plt.subplot
    ax.add_artist(earth_circle)
    ax.set_aspect(aspect="equal")

    # plot body velocity vector
    vp = planetary_rendezvous.central_body_state.velocity_vector
    vp = vp / np.linalg.norm(vp)
    ax.arrow(0, 0, vp[0], vp[1], head_width=0.05, head_length=0.1, fc='k', ec='k')

    # plot parent body unit vector
    r_sun = planetary_rendezvous.central_body_state.position_vector
    r_sun = - r_sun / np.linalg.norm(r_sun)
    ax.arrow(0, 0, r_sun[0], r_sun[1], head_width=0.05, head_length=0.1, fc='k', ec='k')

    # plot excess velocity vector
    v_out = planetary_rendezvous.hyperbolic_excess_velocity
    v_out = v_out / np.linalg.norm(v_out)
    ax.arrow(-v_out[0], -v_out[1], v_out[0], v_out[1], head_width=0.05, head_length=0.1, fc='k', ec='k')

    # sample points for plotting sphere of influence
    theta = np.linspace(0, np.pi * 2, 300)
    x, y = radius_soi * np.cos(theta) * scale, radius_soi * np.sin(theta) * scale

    # plot sphere of influence
    ax.plot(x, y, linestyle=":")

    ax.plot(trajectory_1_cartesian[:, 0] * scale, trajectory_1_cartesian[:, 1] * scale)
    ax.plot(trajectory_2_cartesian[:, 0] * scale, trajectory_2_cartesian[:, 1] * scale)
    if simulated_trajectory:
        ax.plot(simulated_trajectory[:, 0] * scale, simulated_trajectory[:, 1] * scale)

    plt.show()


def sample_keplerian_trajectory(state, t_sample):
    return np.vstack([keplerian_to_cartesian(
        keplerian_elements=propagate_kepler_orbit(
            initial_keplerian_state=state.to_keplerian().state_vector,
            gravitational_parameter=state.gravitational_parameter,
            propagation_time=t),
        gravitational_parameter=state.gravitational_parameter) for t in
        t_sample])
