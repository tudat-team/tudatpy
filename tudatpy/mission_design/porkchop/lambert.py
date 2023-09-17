''' 
Copyright (c) 2010-2020, Delft University of Technology
All rigths reserved

This file is part of the Tudat. Redistribution and use in source and 
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
'''

# General imports
import numpy as np

# Tudat imports
from tudatpy.kernel.astro import two_body_dynamics
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.astro.time_conversion import DateTime
from tudatpy.kernel.numerical_simulation import environment


def calculate_lambert_arc_impulsive_DV(
        bodies: environment.SystemOfBodies,
        global_frame_orientation: str,
        departure_body: str,
        target_body: str,
        departure_epoch: int,
        arrival_epoch: int,
        central_body: str = 'Sun' ) -> float:

    """"
    This function solved Lambert's problem for a transfer from Earth (at departure epoch) to
    a target body (at arrival epoch), with the states of Earth and the target body defined
    by ephemerides stored inside the SystemOfBodies object (bodies). Note that this solver
    assumes that the transfer departs/arrives to/from the center of mass of Earth and the target body

    Parameters
    ----------
    bodies : Body objects defining the physical simulation environment

    target_body : The name (string) of the body to which the Lambert arc is to be computed

    departure_epoch : Epoch at which the departure from Earth's center of mass is to take place

    arrival_epoch : Epoch at which the arrival at he target body's center of mass is to take place

    Return
    ------
    ΔV required for insertion into the Lambert transfer arc and capture by the arrival body
    """

    # Gravitational parameter of the Sun
    central_body_gravitational_parameter = bodies.get_body(central_body).gravitational_parameter

    # Retrieve states of departure and arrival body
    initial_state = spice_interface.get_body_cartesian_state_at_epoch(
        target_body_name       = departure_body,
        observer_body_name     = central_body,
        reference_frame_name   = global_frame_orientation,
        aberration_corrections = "NONE",
        ephemeris_time         = departure_epoch)
    final_state = spice_interface.get_body_cartesian_state_at_epoch(
        target_body_name       = target_body,
        observer_body_name     = central_body,
        reference_frame_name   = global_frame_orientation,
        aberration_corrections = "NONE",
        ephemeris_time         = arrival_epoch)

    # Retrieve initial and final positions for Lambert targeter
    departure_position = initial_state[:3]
    arrival_position   = final_state[:3]

    # Create Lambert targeter
    lambertTargeter = two_body_dynamics.LambertTargeterIzzo(
        departure_position,
        arrival_position, 
        arrival_epoch - departure_epoch,
        central_body_gravitational_parameter)
    
    # Obtain initial and final states
    departure_velocity = lambertTargeter.get_departure_velocity()
    arrival_velocity   = lambertTargeter.get_arrival_velocity()

    # Retrieve initial and insertion velocities to calculate ΔV
    initial_velocity   = initial_state[3:]
    insertion_velocity = final_state[3:]

    # Calculate ΔV
    ΔV_launch  = np.linalg.norm(departure_velocity - initial_velocity)
    ΔV_arrival = np.linalg.norm(arrival_velocity - insertion_velocity)

    return ΔV_launch, ΔV_arrival