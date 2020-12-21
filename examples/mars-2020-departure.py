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
from tudatpy.kernel.simulation import propagation_setup
from tudatpy.kernel.astro.two_body_dynamics.iod import (
    LambertTargeterIzzo,
    PlanetaryDeparture,
    PlanetaryRendezvous,

)

from tudatpy.kernel.astro import frames
from tudatpy.kernel.astro.frames import ReferenceFrame

from tudatpy.kernel.astro.two_body_dynamics.state import (
    KeplerianState,
    CartesianState,
)

from tudatpy.kernel.bodies.spice.solar_system import (
    Earth,
    Mars,
    Sun
)

from tudatpy.plotting import (
    plot_planetary_departure,
    plot_planetary_rendezvous,
    print_departure,
    print_rendezvous,
    sample_keplerian_trajectory)

import numpy as np
import matplotlib.pyplot as plt
from tudatpy.kernel.astro.conversion import keplerian_to_cartesian

from tudatpy.kernel import constants

if __name__ == "__main__":
    spice_kernels = spice_interface.get_standard_kernels() + [
        "/home/ggarrett/Downloads/de430.bsp",
        "/home/ggarrett/Downloads/de431_part-1.bsp",
        "/home/ggarrett/Downloads/de431_part-2.bsp"]

    for kernel in spice_kernels:
        spice_interface.load_kernel(kernel)

    # LEGACY
    # define reference frame
    # frame_origin = "SSB"
    # frame_orientation = "ECLIPJ2000"

    # NEW DESIGN: DEFAULT
    reference_frame = frames.SSB_ECLIPJ2000
    print("frame_origin: ", reference_frame.origin)
    print("frame_orientation: ", reference_frame.orientation)

    # NEW DESIGN: [ALTERNATIVE] DEFINITION OF NEW FRAME - EXTRACTED FROM SPICE
    # reference_frame = ReferenceFrame(  # SSB_ECLIPJ2000
    #     origin=frame_origin,
    #     orientation=frame_orientation)
    print(reference_frame)

    # mission system parent body
    parent_body = Sun

    # Mission planetary departure parameter definition
    departure_date = "July 30, 2020"
    departure_body = Earth

    # Mission planetary arrival parameter definition
    rendezvous_date = "18 February, 2021"
    rendezvous_body = Mars

    # Calculate ephemeris time and time of flight parameters
    departure_reference_ephemeris_time = spice_interface.convert_date_string_to_ephemeris_time(departure_date)
    rendezvous_reference_ephemeris_time = spice_interface.convert_date_string_to_ephemeris_time(rendezvous_date)

    # Define departure body state at departure conic reference
    departure_body_state_at_reference = CartesianState.from_spice(
        body=departure_body,
        ephemeris_time=departure_reference_ephemeris_time,
        # gravitational_parameter=parent_body_mu,
        reference_frame=reference_frame,
        aberration_corrections="none"
    )

    # departure_body_state_at_reference = CartesianState.from_spice(
    #     body_name=departure_body,
    #     ephemeris_time=departure_reference_ephemeris_time,
    #     gravitational_parameter=parent_body_mu,
    #     frame_origin=frame_origin,
    #     frame_orientation=frame_orientation,
    #     aberration_corrections="none"
    # )

    # Determine rendezvous body state at rendezvous conic reference
    rendezvous_body_state_at_reference = CartesianState.from_spice(
        # body_name=rendezvous_body,
        central_body=rendezvous_body,
        ephemeris_time=rendezvous_reference_ephemeris_time,
        # gravitational_parameter=parent_body_mu,
        reference_frame=reference_frame,
        aberration_corrections="none")

    # Calculate coast interplanetary leg
    inteplanetary_coast = InterplanetaryCoast(
        central_body=Sun,

    )

    departure_velocity, rendezvous_velocity = LambertTargeterIzzo(
        departure_position=departure_body_state_at_reference.position_vector,
        arrival_position=rendezvous_body_state_at_reference.position_vector,
        time_of_flight=rendezvous_reference_ephemeris_time - departure_reference_ephemeris_time,
        gravitational_parameter=Sun.gravitational_parameter
    ).get_velocity_vectors()

    # Calculate departure parameters
    planetary_departure = PlanetaryDeparture(
        outgoing_velocity=departure_velocity,
        central_body_state=departure_body_state_at_reference,
        periapsis_distance=departure_body_radius * departure_r_limit + 100 * 1e3,
        gravitational_parameter=departure_body_mu,
        prograde_orbit=True)

    plot_planetary_departure(planetary_departure)
    print_departure(planetary_departure)

    start_epoch = departure_reference_ephemeris_time
    end_epoch = rendezvous_reference_ephemeris_time

    ###########################################################################
    # PLOT ANALYTICAL ORBIT ###################################################
    ###########################################################################

    # Internal
    # spacecraft_trajectory_departure = planetary_departure.outgoing_orbit.sample(reference_frame=reference_frame)
    # spacecraft_trajectory_interplanetary = interplanetary_coast.orbit.sample(reference_frame=reference_frame)
    # departure_body_during_departure = departure_body_orbit.sample(reference_frame=reference_frame)

    orbit_plotter = OrbitPlotter2D(reference_frame=departure_body.local_frame, backend="matplotlib")
    orbit_plotter.plot(planetary_departure.incoming_orbit)
    orbit_plotter.plot(planetary_departure.outgoing_orbit)
    orbit_plotter.save_fig("planetary_departure.pdf")

    orbit_plotter = OrbitPlotter2D(reference_frame=frames.SSB_ECLIPJ2000, backend="matplotlib")
    orbit_plotter.plot(inteplanetary_coast.orbit)
    orbit_plotter.save_fig("inteplanetary_coast.pdf")

    Orbit.from_spice(departure_body)

    departure_body_trajectory = sample_keplerian_trajectory(
        departure_body_state_at_reference,
        time_samples)

    departure_body_trajectory_during_departure = sample_keplerian_trajectory(
        departure_body_state_at_reference,
        time_samples_departure)

    rendezvous_body_trajectory = sample_keplerian_trajectory(
        rendezvous_body_state_at_reference,
        -time_samples[::-1])

    scale = constants.AU

    sun_circle = plt.Circle((0, 0), 0.2)
    departure_body_circle = plt.Circle(
        (departure_body_state_at_reference.position_vector[0] / scale,
         departure_body_state_at_reference.position_vector[1] / scale),
        departure_body_radius / scale)
    rendezvous_body_circle = plt.Circle(
        (rendezvous_body_state_at_reference.position_vector[0] / scale,
         rendezvous_body_state_at_reference.position_vector[1] / scale),
        rendezvous_body_radius / scale)

    fig = plt.gcf()
    ax = fig.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.add_artist(sun_circle)
    ax.add_artist(departure_body_circle)
    ax.add_artist(rendezvous_body_circle)
    ax.plot((departure_conic_spacecraft_trajectory + departure_body_trajectory_during_departure)[:, 0] / scale,
            (departure_body_trajectory_during_departure + departure_conic_spacecraft_trajectory)[:, 1] / scale)
    ax.plot(coast_conic_spacecraft_trajectory[:, 0] / scale, coast_conic_spacecraft_trajectory[:, 1] / scale)
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
