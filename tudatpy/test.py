###############################################################################
# IMPORT STATEMENTS ###########################################################
###############################################################################
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.astro.two_body_dynamics import (
    calculate_departure_parameters,
    calculate_arrival_parameters,
    calculate_gravity_assist_parameters,
    LambertTargeterIzzo)

if __name__ == "__main__":
    spice_interface.load_standard_kernels()

    """
    Case Study: Mars 2020
    ---------------------
    Mars 2020 is a Mars rover mission by NASA's Mars Exploration Program that 
    includes the Perseverance rover and the Ingenuity helicopter drone. It was
    launched on 30 July 2020 at 11:50 UTC, and will touch down in Jezero
    crater on Mars on 18 February 2021.
    """

    # Mission departure parameter definition
    departure_date = "July 30, 2020"
    departure_body = "Earth"
    departure_altitude = 100 * 1e3

    # Mission arrival parameter definition
    arrival_date = "18 February, 2021"
    arrival_body = "Mars"
    arrival_altitude = 100 * 1e3

    # Calculate ephemeris time and time of flight parameters
    departure_ephemeris_time = spice_interface.convert_date_string_to_ephemeris_time(departure_date)
    arrival_ephemeris_time = spice_interface.convert_date_string_to_ephemeris_time(arrival_date)

    # Determine departure geometry
    departure_body_mu = spice_interface.get_body_gravitational_parameter(departure_body)
    departure_body_radius = spice_interface.get_average_radius(departure_body)
    departure_body_state = spice_interface.get_body_cartesian_state_at_epoch(
        target_body_name=departure_body,
        observer_body_name="SSB",
        reference_frame_name="J2000",
        aberration_corrections="none",
        ephemeris_time=departure_ephemeris_time)
    departure_body_velocity = departure_body_state[3:]
    departure_body_position = departure_body_state[:3]

    # Determine arrival geometry
    arrival_body_mu = spice_interface.get_body_gravitational_parameter(arrival_body)
    arrival_body_radius = spice_interface.get_average_radius(arrival_body)
    arrival_body_state = spice_interface.get_body_cartesian_state_at_epoch(
        target_body_name=arrival_body,
        observer_body_name="SSB",
        reference_frame_name="J2000",
        aberration_corrections="none",
        ephemeris_time=arrival_ephemeris_time)
    arrival_body_velocity = arrival_body_state[3:]
    arrival_body_position = arrival_body_state[:3]

    # Calculate interplanetary coast leg
    departure_velocity, arrival_velocity = LambertTargeterIzzo(
        departure_position=departure_body_position,
        arrival_position=arrival_body_position,
        time_of_flight=arrival_ephemeris_time - departure_ephemeris_time,
        gravitational_parameter=spice_interface.get_body_gravitational_parameter("Sun")
    ).get_velocity_vectors()

    # Calculate departure parameters
    departure_parameters = calculate_departure_parameters(
        central_body_gravitational_parameter=departure_body_mu,
        central_body_velocity_on_exit_soi=departure_body_velocity,
        outgoing_velocity=departure_velocity,
        departure_periapsis_distance=departure_body_radius + departure_altitude)

    # Calculate arrival parameters
    arrival_parameters = calculate_arrival_parameters(
        central_body_gravitational_parameter=arrival_body_mu,
        central_body_velocity_on_entry_soi=arrival_body_velocity,
        incoming_velocity=arrival_velocity,
        arrival_periapsis_distance=arrival_body_radius + arrival_altitude)

    print("departure_parameters.radius_of_periapsis: ", departure_parameters.radius_of_periapsis)
    print("departure_parameters.initial_velocity_at_periapsis: ", departure_parameters.initial_velocity_at_periapsis)
    print("departure_parameters.outgoing_semi_major_axis: ", departure_parameters.outgoing_semi_major_axis)
    print("departure_parameters.outgoing_eccentricity: ", departure_parameters.outgoing_eccentricity)
    print("departure_parameters.final_velocity_at_periapsis: ", departure_parameters.final_velocity_at_periapsis)
    print("departure_parameters.total_delta_v: ", departure_parameters.total_delta_v)
    print('\n')
    print("arrival_parameters.radius_of_periapsis: ", arrival_parameters.radius_of_periapsis)
    print("arrival_parameters.initial_velocity_at_periapsis: ", arrival_parameters.initial_velocity_at_periapsis)
    print("arrival_parameters.incoming_semi_major_axis: ", arrival_parameters.incoming_semi_major_axis)
    print("arrival_parameters.incoming_eccentricity: ", arrival_parameters.incoming_eccentricity)
    print("arrival_parameters.final_velocity_at_periapsis: ", arrival_parameters.final_velocity_at_periapsis)
    print("arrival_parameters.total_delta_v: ", arrival_parameters.total_delta_v)
