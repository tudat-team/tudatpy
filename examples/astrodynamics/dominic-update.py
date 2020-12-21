from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.astro.fundamentals import propagate_kepler_orbit
from tudatpy.kernel.simulation import environment_setup
from tudatpy.kernel.astro.two_body_dynamics import (
    # calculate_departure_parameters,
    # calculate_arrival_parameters,
    # calculate_gravity_assist_parameters,
    LambertTargeterIzzo,
    cartesian_state_from_spice,
    PlanetaryDeparture,
    PlanetaryRendezvous,
    KeplerianState,
    CartesianState
)

from decouple import config

from tudatpy.plotting import plot_planetary_rendezvous
from tudatpy.plotting import plot_planetary_departure

# retrieve spice kernel paths from .env file
SPICE_KERNEL_BSP_DE430 = config("SPICE_KERNEL_BSP_DE430")
SPICE_KERNEL_BSP_DE431_PART1 = config("SPICE_KERNEL_BSP_DE431_PART1")
SPICE_KERNEL_BSP_DE431_PART2 = config("SPICE_KERNEL_BSP_DE431_PART2")

# append the additional spice kernels to the standard ones
spice_kernels = spice_interface.get_standard_kernels() + [
    SPICE_KERNEL_BSP_DE430,
    SPICE_KERNEL_BSP_DE431_PART1,
    SPICE_KERNEL_BSP_DE431_PART2]

for kernel in spice_kernels:
    spice_interface.load_kernel(kernel)

# mission system parent body
parent_body = "Sun"
parent_body_mu = spice_interface.get_body_gravitational_parameter(parent_body)

# define reference frame
frame_origin = parent_body
# frame_origin = "SSB"
frame_orientation = "ECLIPJ2000"

# Mission planetary departure parameter definition
departure_date = "July 30, 2020"
departure_body = "Earth"
departure_r_limit = 5.0

# Mission planetary arrival parameter definition
rendezvous_date = "18 February, 2026"
rendezvous_body = "Jupiter"
rendezvous_r_limit = 5.0

# Calculate ephemeris time and time of flight parameters
departure_reference_ephemeris_time = spice_interface.convert_date_string_to_ephemeris_time(departure_date)
rendezvous_reference_ephemeris_time = spice_interface.convert_date_string_to_ephemeris_time(rendezvous_date)

# Determine constant departure parameters
departure_body_mu = spice_interface.get_body_gravitational_parameter(departure_body)
departure_body_radius = spice_interface.get_average_radius(departure_body)

# Determine constant arrival geometry
rendezvous_body_mu = spice_interface.get_body_gravitational_parameter(rendezvous_body)
rendezvous_body_radius = spice_interface.get_average_radius(rendezvous_body)


# Calculate departure body state at periapsis passage
def get_departure_state_new_api():
    departure_body_cartesian_state_vector = spice_interface.get_body_cartesian_state_at_epoch(
        target_body_name=departure_body,
        ephemeris_time=departure_reference_ephemeris_time,
        observer_body_name=frame_origin,
        reference_frame_name=frame_orientation,
        aberration_corrections="none"
    )

    departure_body_state_at_periapsis = CartesianState(
        position_vector=departure_body_cartesian_state_vector[:3],
        velocity_vector=departure_body_cartesian_state_vector[3:],
        gravitational_parameter=parent_body_mu,
    )
    return departure_body_state_at_periapsis


def get_departure_state_old_api():
    # Define departure body state at departure conic reference
    departure_body_state_at_periapsis = cartesian_state_from_spice(
        body_name=departure_body,
        ephemeris_time=departure_reference_ephemeris_time,
        gravitational_parameter=parent_body_mu,
        frame_origin=frame_origin,
        frame_orientation=frame_orientation,
        aberration_corrections="none"
    )
    return departure_body_state_at_periapsis

    # Ideally
    # departure_body_state_at_periapsis = CartesianState.from_spice(
    #     body_name=departure_body,
    #     ephemeris_time=departure_reference_ephemeris_time,
    #     gravitational_parameter=parent_body_mu,
    #     frame_origin=frame_origin,
    #     frame_orientation=frame_orientation,
    #     aberration_corrections="none"
    # )


departure_body_state_at_periapsis = get_departure_state_old_api()
departure_body_state_at_periapsis = get_departure_state_new_api()

# Define departure body state at departure conic reference
rendezvous_body_state_at_periapsis = cartesian_state_from_spice(
    body_name=rendezvous_body,
    ephemeris_time=rendezvous_reference_ephemeris_time,
    gravitational_parameter=parent_body_mu,
    frame_origin=frame_origin,
    frame_orientation=frame_orientation,
    aberration_corrections="none"
)

# Calculate coast interplanetary leg
departure_velocity, rendezvous_velocity = LambertTargeterIzzo(
    departure_position=departure_body_state_at_periapsis.position_vector,
    arrival_position=rendezvous_body_state_at_periapsis.position_vector,
    time_of_flight=rendezvous_reference_ephemeris_time - departure_reference_ephemeris_time,
    gravitational_parameter=parent_body_mu
).get_velocity_vectors()

# Calculate planetary departure leg
planetary_departure = PlanetaryDeparture(
    outgoing_velocity=departure_velocity,
    central_body_state=departure_body_state_at_periapsis,
    periapsis_distance=departure_body_radius * departure_r_limit,
    gravitational_parameter=departure_body_mu,
    prograde_orbit=True)

planetary_rendezvous = PlanetaryRendezvous(
    incoming_velocity=rendezvous_velocity,
    central_body_state=rendezvous_body_state_at_periapsis,
    periapsis_distance=rendezvous_r_limit * rendezvous_body_radius,
    gravitational_parameter=rendezvous_body_mu,
    prograde_orbit=True)

plot_planetary_departure(planetary_departure)
plot_planetary_rendezvous(planetary_rendezvous)
