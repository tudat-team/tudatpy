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
from tudatpy.kernel.simulation import environment_setup
from tudatpy.kernel.astro.two_body_dynamics import (
    calculate_departure_parameters,
    calculate_arrival_parameters,
    calculate_gravity_assist_parameters,
    cartesian_state_from_spice,
    LambertTargeterIzzo,
    PlanetaryDeparture,
    PlanetaryRendezvous,
    KeplerianState,
    CartesianState
)

import numpy as np
import matplotlib.pyplot as plt
from tudatpy.kernel.astro.conversion import keplerian_to_cartesian


if __name__ == "__main__":
    spice_kernels = spice_interface.get_standard_kernels() + [
        "/home/ggarrett/Downloads/de430.bsp",
        "/home/ggarrett/Downloads/de431_part-1.bsp",
        "/home/ggarrett/Downloads/de431_part-2.bsp"]

    for kernel in spice_kernels:
        spice_interface.load_kernel(kernel)

    import matplotlib.pyplot as plt
    fig = plt.gcf()
    ax = fig.gca()

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
    #
    # bodies_to_create = [departure_body, rendezvous_body, parent_body]
    #
    # # Create default body settings for bodies_to_create, with "SSB"/"J2000" as
    # # global frame origin and orientation
    # body_settings = environment_setup.get_default_body_settings(
    #     bodies_to_create, frame_origin, frame_orientation
    # )
    #
    # # Create system of bodies (in this case only Earth)
    # bodies = environment_setup.create_system_of_bodies(body_settings)

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
        aberration_corrections="none"
    )

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
        time_of_flight=rendezvous_reference_ephemeris_time-departure_reference_ephemeris_time,
        gravitational_parameter=parent_body_mu
    ).get_velocity_vectors()

    # Calculate departure parameters
    planetary_departure = PlanetaryDeparture(
        outgoing_velocity=departure_velocity,
        central_body_state=departure_body_state_at_reference,
        periapsis_distance=departure_body_radius * departure_r_limit + 500 * 1e3,
        gravitational_parameter=departure_body_mu,
        prograde_orbit=True)

    # Calculate planetary rendezvous geometry
    planetary_rendezvous = PlanetaryRendezvous(
        incoming_velocity=rendezvous_velocity,
        central_body_state=rendezvous_body_state_at_reference,
        periapsis_distance=rendezvous_r_limit * rendezvous_body_radius + 500 * 1e3,
        gravitational_parameter=rendezvous_body_mu,
        prograde_orbit=True
    )

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