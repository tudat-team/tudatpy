''' 
Copyright (c) 2010-2023, Delft University of Technology
All rigths reserved

This file is part of the Tudat. Redistribution and use in source and 
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
'''

# General imports
import numpy as np
from tqdm import tqdm
from numbers import Number

# Tudat imports
from tudatpy.kernel import constants
from tudatpy.kernel.astro import time_conversion
from tudatpy.kernel.numerical_simulation import environment
from tudatpy.trajectory_design.porkchop import plot_porkchop
from tudatpy.trajectory_design.porkchop.lambert import calculate_lambert_arc_impulsive_delta_v


def determine_shape_of_delta_v(
        bodies: environment.SystemOfBodies,
        global_frame_orientation: str,
        departure_body: str,
        target_body: str,
        departure_epoch: float,
        arrival_epoch: float,
        function_to_calculate_delta_v: calculate_lambert_arc_impulsive_delta_v
    ):
    """
    Determine whether `function_to_calculate_delta_v` returns ΔV as
    
    * A single float representing the total ΔV of the transfer
    * A `list`/`tuple`/`np.ndarray` containing [departure ΔV, arrival ΔV]

    Arguments
    ---------
    * 

    Output
    ------
    * `shape`: 1 if the ΔV returned by the `function_to_calculate_delta_v`
               is a `float`, 2 if it is a list/tuple/NumPy array
               containing [departure ΔV, arrival ΔV]
    """

    ΔV_sample = function_to_calculate_delta_v(
        bodies, global_frame_orientation,
        departure_body,
        target_body,
        departure_epoch,
        arrival_epoch
    )

    # Assert ΔV is either a numeric type, or a list, tuple or NumPy array
    ΔV_is_numeric          = isinstance(ΔV_sample, Number)
    ΔV_is_dep_arr_ΔV_tuple = type(ΔV_sample) in [tuple, list, np.ndarray]
    if not (ΔV_is_numeric or ΔV_is_dep_arr_ΔV_tuple):
        raise TypeError(
            f'The return type of `function_to_calculate_delta_v` is [{type(ΔV_sample).__name__}], which is invalid. '
            f'Ensure the return type of `function_to_calculate_delta_v` is either:\n\n'
            f'    - Numeric, representing the **TOTAL ΔV** required by the transfer\n'
            f'    - A list/tuple/NumPy array of length 2 containing **[DEPARTURE ΔV, ARRIVAL ΔV]**\n'
        )

    # If ΔV is being returned as a tuple, list or NumPy array, assert its length is either 1 or 2
    if ΔV_is_dep_arr_ΔV_tuple and len(ΔV_sample) != 2:
        raise ValueError(
            f'The length of the output of `function_to_calculate_delta_v` is invalid ({len(ΔV_sample)}). '
            f'Ensure the return type of `function_to_calculate_delta_v` is either:\n\n'
            f'    - Numeric, representing the **TOTAL ΔV** required by the transfer\n'
            f'    - A list/tuple/NumPy array of length 2 containing **[DEPARTURE ΔV, ARRIVAL ΔV]**\n'
        )

    return 1 if ΔV_is_numeric else 2


def calculate_delta_v_time_map(
        bodies: environment.SystemOfBodies,
        global_frame_orientation: str,
        departure_body: str,
        target_body: str,
        earliest_departure_time: time_conversion.DateTime,
        latest_departure_time: time_conversion.DateTime,
        earliest_arrival_time: time_conversion.DateTime,
        latest_arrival_time: time_conversion.DateTime,
        time_resolution: float,
        function_to_calculate_delta_v: calculate_lambert_arc_impulsive_delta_v
    ):
    """
    Creates an array containing the ΔV of all coordinates of the grid of departure/arrival epochs.

    Arguments
    ---------

    Output
    ------

    """

    # Input validation
    supported_global_frame_orientations = ['J2000', 'ECLIPJ2000']
    assert global_frame_orientation in supported_global_frame_orientations, \
        f'\n    The `global_frame_orientation` provided ("{global_frame_orientation}") is not supported by Tudat.\n' \
        f'\n    Please provide a `global_frame_orientation` in {supported_global_frame_orientations}'

    # Departure and arrival epoch discretizations
    n_dep = int(
        (latest_departure_time.epoch() - earliest_departure_time.epoch()) / (time_resolution * constants.JULIAN_DAY)
    )
    n_arr = int(
        (latest_arrival_time.epoch() - earliest_arrival_time.epoch()) / (time_resolution * constants.JULIAN_DAY)
    )
    departure_epochs = np.linspace(
        earliest_departure_time.epoch(), 
        latest_departure_time.epoch(),
        n_dep
    )
    arrival_epochs = np.linspace(
        earliest_arrival_time.epoch(), 
        latest_arrival_time.epoch(), 
        n_arr
    )

    # Determine the shape of the ΔV returned by the user-provided `function_to_calculate_delta_v`
    ΔV_shape = determine_shape_of_delta_v(
        bodies, global_frame_orientation,
        departure_body,
        target_body,
        departure_epochs[0],
        arrival_epochs[0],
        function_to_calculate_delta_v
    )

    # Pre-allocate ΔV array
    ΔV = np.full((n_dep, n_arr, ΔV_shape), np.nan)

    for i_dep in tqdm(range(n_dep)):
        for i_arr in range(n_arr):
            # Calculate ΔV, only if the arrival epoch is greater than the departure epoch
            if arrival_epochs[i_arr] > departure_epochs[i_dep]:
                ΔV[i_dep, i_arr, :] = function_to_calculate_delta_v(
                    bodies, global_frame_orientation,
                    departure_body,
                    target_body,
                    departure_epochs[i_dep],
                    arrival_epochs[i_arr]
                )

    return [departure_epochs, arrival_epochs, ΔV]


def porkchop(
        # ΔV calculation arguments
        bodies: environment.SystemOfBodies,
        global_frame_orientation: str,
        departure_body: str,
        target_body: str,
        earliest_departure_time: time_conversion.DateTime,
        latest_departure_time: time_conversion.DateTime,
        earliest_arrival_time: time_conversion.DateTime,
        latest_arrival_time: time_conversion.DateTime,
        time_resolution: float,
        function_to_calculate_delta_v: calculate_lambert_arc_impulsive_delta_v,
        # Plot arguments
        C3: bool = False,
        total: bool = False,
        threshold: float = 10,
        upscale: bool = False,
        number_of_levels: int = 10,
        # Figure arguments
        percent_margin: float = 5,
        figsize: tuple[int, int] = (8, 8),
        show: bool = True,
        save: bool = False,
        filename: str = 'porkchop.png',
    ) -> None:
    """Tudat ΔV/C3 porkchop mission design plot.
    
    Arguments
    ---------

    Output
    ------
    
    """

    # Calculate ΔV map
    [departure_epochs, arrival_epochs, ΔV] = calculate_delta_v_time_map(
        bodies                        = bodies,
        global_frame_orientation      = global_frame_orientation,
        departure_body                = departure_body,
        target_body                   = target_body,
        earliest_departure_time       = earliest_departure_time,
        latest_departure_time         = latest_departure_time,
        earliest_arrival_time         = earliest_arrival_time,
        latest_arrival_time           = latest_arrival_time,
        time_resolution               = time_resolution,
        function_to_calculate_delta_v = function_to_calculate_delta_v
    )

    # Plot porkchop
    plot_porkchop(
        departure_body   = departure_body,
        target_body      = target_body,
        departure_epochs = departure_epochs,
        arrival_epochs   = arrival_epochs, 
        delta_v          = ΔV,
        C3               = C3,
        total            = total,
        threshold        = threshold,
        # Plot arguments
        upscale          = upscale,
        number_of_levels = number_of_levels,
        # Figure arguments
        percent_margin   = percent_margin,
        figsize          = figsize,
        show             = show, 
        save             = save, 
        filename         = filename,
    )

    return [departure_epochs, arrival_epochs, ΔV]