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
from ...astro import two_body_dynamics
from ...numerical_simulation import environment


def calculate_lambert_arc_impulsive_delta_v(
        bodies: environment.SystemOfBodies,
        departure_body: str,
        target_body: str,
        departure_epoch: int,
        arrival_epoch: int,
        central_body: str = 'Sun' ) -> tuple[float]:

    """"
    This function solved Lambert's problem for a transfer from the `departure_body` (at departure epoch)
    to a `target_body` (at arrival epoch), with the states of the `departure_body` and the `target_body`
    defined by ephemerides stored inside the `bodies` (`SystemOfBodies` instance). Note that this solver
    assumes that the transfer departs/arrives to/from the center of mass of the departure and the target body.

    Parameters
    ----------
    bodies: environment.SystemOfBodies
        Body objects defining the physical simulation environment
    departure_body: str
        The name of the body from which the transfer is to be computed
    target_body: str
        The name of the body to which the transfer is to be computed
    departure_epoch: int
        Epoch at which the departure from the `target_body`'s center of mass is to take place
    arrival_epoch: int
        Epoch at which the arrival at he target body's center of mass is to take place

    Output
    ------
    ΔV_launch: float
        ΔV required for insertion into the Lambert transfer arc
    ΔV_arrival: float
        ΔV required for capture by the arrival body
    """

    # Gravitational parameter of the Sun
    central_body_gravitational_parameter = bodies.get_body(central_body).gravitational_parameter

    # Retrieve states of departure and arrival body
    initial_state = bodies.get_body(departure_body).state_in_base_frame_from_ephemeris(departure_epoch)
    final_state = bodies.get_body(target_body).state_in_base_frame_from_ephemeris(arrival_epoch)

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
