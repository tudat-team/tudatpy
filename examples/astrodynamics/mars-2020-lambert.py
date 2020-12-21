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
from tudatpy.kernel.astro.two_body_dynamics import (
    cartesian_state_from_spice,
    LambertTargeterIzzo,
    CartesianState
)

from tudatpy.plotting import (
    sample_keplerian_trajectory)

import numpy as np
import matplotlib.pyplot as plt

from tudatpy.kernel import constants

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
    departure_body_state_at_reference = cartesian_state_from_spice(
        body_name=departure_body,
        ephemeris_time=departure_reference_ephemeris_time,
        gravitational_parameter=parent_body_mu,
        frame_origin=frame_origin,
        frame_orientation=frame_orientation,
        aberration_corrections="none")

    # Determine rendezvous body state at rendezvous conic reference
    rendezvous_body_state_at_reference = cartesian_state_from_spice(
        body_name=rendezvous_body,
        ephemeris_time=rendezvous_reference_ephemeris_time,
        gravitational_parameter=parent_body_mu,
        frame_origin=frame_origin,
        frame_orientation=frame_orientation,
        aberration_corrections="none")

    # Calculate coast interplanetary leg
    departure_velocity, rendezvous_velocity = LambertTargeterIzzo(
        departure_position=departure_body_state_at_reference.position_vector,
        arrival_position=rendezvous_body_state_at_reference.position_vector,
        time_of_flight=rendezvous_reference_ephemeris_time - departure_reference_ephemeris_time,
        gravitational_parameter=parent_body_mu
    ).get_velocity_vectors()

    start_epoch = departure_reference_ephemeris_time
    end_epoch = rendezvous_reference_ephemeris_time

    ###########################################################################
    # PLOT ANALYTICAL ORBIT ###################################################
    ###########################################################################

    n_samples = 400
    time_samples = np.linspace(0, end_epoch - start_epoch, n_samples)

    spacecraft_trajectory = sample_keplerian_trajectory(
        CartesianState(
            position_vector=departure_body_state_at_reference.position_vector,
            velocity_vector=departure_velocity,
            gravitational_parameter=parent_body_mu
        ),
        time_samples)

    departure_body_trajectory = sample_keplerian_trajectory(
        departure_body_state_at_reference,
        time_samples)

    rendezvous_body_trajectory = sample_keplerian_trajectory(
        rendezvous_body_state_at_reference,
        -time_samples[::-1])

    scale = constants.AU

    sun_circle = plt.Circle((0, 0), 0.2)
    departure_body_circle = plt.Circle(
        (departure_body_state_at_reference.position_vector[0] / scale,  # x0
         departure_body_state_at_reference.position_vector[1] / scale),  # y0
        departure_body_radius / scale)  # radius
    rendezvous_body_circle = plt.Circle(
        (rendezvous_body_state_at_reference.position_vector[0] / scale,  # x0
         rendezvous_body_state_at_reference.position_vector[1] / scale),  # y0
        rendezvous_body_radius / scale)  # radius

    fig = plt.gcf()
    ax = fig.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.add_artist(sun_circle)
    ax.add_artist(departure_body_circle)
    ax.add_artist(rendezvous_body_circle)
    ax.plot(spacecraft_trajectory[:, 0] / scale, spacecraft_trajectory[:, 1] / scale)
    ax.plot(departure_body_trajectory[:, 0] / scale, departure_body_trajectory[:, 1] / scale)
    ax.plot(rendezvous_body_trajectory[:, 0] / scale, rendezvous_body_trajectory[:, 1] / scale)
    ax.plot([0], [0])
    plt.show()

    spice_kernel_repr = "\n".join(spice_kernels)
    print(
        f"""
Spice Kernels Loaded
====================
{spice_kernel_repr}

""")
