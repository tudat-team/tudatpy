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
from tudatpy.trajectory_design.porkchop._plot_porkchop import plot_porkchop
from tudatpy.trajectory_design.porkchop._lambert import calculate_lambert_arc_impulsive_delta_v


def determine_shape_of_delta_v(
        bodies: environment.SystemOfBodies,
        departure_body: str,
        target_body: str,
        departure_epoch: float,
        arrival_epoch: float,
        function_to_calculate_delta_v: callable = calculate_lambert_arc_impulsive_delta_v
):
    """
    Determine whether `function_to_calculate_delta_v` returns ΔV as
    
    * A single float representing the total ΔV of the transfer
    * A `list`/`tuple`/`np.ndarray` containing [departure ΔV, arrival ΔV]

    Parameters
    ----------
    bodies: environment.SystemOfBodies
        Body objects defining the physical simulation environment
    departure_body: str
        The name of the body from which the transfer is to be computed
    target_body: str
        The name of the body to which the transfer is to be computed
    departure_epoch: float
        Epoch at which the departure from the `target_body`'s center of mass is to take place
    arrival_epoch: float
        Epoch at which the arrival at he target body's center of mass is to take place
    function_to_calculate_delta_v: callable = calculate_lambert_arc_impulsive_delta_v
        Function with which the manoeuvre's required ΔV will be calculated

    Output
    ------
    shape: int
        1 if the ΔV returned by the `function_to_calculate_delta_v` is a `float`, 
        2 if it is a list/tuple/NumPy array containing [departure ΔV, arrival ΔV]
    """

    ΔV_sample = function_to_calculate_delta_v(
        bodies,
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
        departure_body: str,
        target_body: str,
        earliest_departure_time: time_conversion.DateTime,
        latest_departure_time: time_conversion.DateTime,
        earliest_arrival_time: time_conversion.DateTime,
        latest_arrival_time: time_conversion.DateTime,
        time_resolution: float,
        function_to_calculate_delta_v: callable = calculate_lambert_arc_impulsive_delta_v
):
    """
    Creates an array containing the ΔV of all coordinates of the grid of departure/arrival epochs.

    Parameters
    ----------
    bodies: environment.SystemOfBodies
        Body objects defining the physical simulation environment
    departure_body: str
        The name of the body from which the transfer is to be computed
    target_body: str
        The name of the body to which the transfer is to be computed
    earliest_departure_time: time_conversion.DateTime
        Earliest epoch of the departure window
    latest_departure_time: time_conversion.DateTime
        Latest epoch of the departure window
    earliest_arrival_time: time_conversion.DateTime
        Earliest epoch of the arrival window
    latest_arrival_time: time_conversion.DateTime
        Latest epoch of the arrival window
    time_resolution: float
        Resolution used to discretize the departure/arrival time windows
    function_to_calculate_delta_v: callable = calculate_lambert_arc_impulsive_delta_v
        Function with which the manoeuvre's required ΔV will be calculated

    Output
    ------
    departure_epochs: 
        Discretized departure time window
    arrival_epochs: 
        Discretized arrival time window
    ΔV: 
        Array containing the ΔV of all coordinates of the grid of departure/arrival epochs
    """

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
        bodies,
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
                    bodies,
                    departure_body,
                    target_body,
                    departure_epochs[i_dep],
                    arrival_epochs[i_arr]
                )

    return [departure_epochs, arrival_epochs, ΔV]


def porkchop(
        # ΔV calculation arguments
        bodies: environment.SystemOfBodies,
        departure_body: str,
        target_body: str,
        earliest_departure_time: time_conversion.DateTime,
        latest_departure_time: time_conversion.DateTime,
        earliest_arrival_time: time_conversion.DateTime,
        latest_arrival_time: time_conversion.DateTime,
        time_resolution: float,
        function_to_calculate_delta_v: callable = calculate_lambert_arc_impulsive_delta_v,
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
    """
    Calculates and displays ΔV/C3 porkchop mission design plots.
    
    Parameters
    ----------
    bodies: environment.SystemOfBodies
        Body objects defining the physical simulation environment
    departure_body: str
        The name of the body from which the transfer is to be computed
    target_body: str
        The name of the body to which the transfer is to be computed
    earliest_departure_time: time_conversion.DateTime
        Earliest epoch of the departure window
    latest_departure_time: time_conversion.DateTime
        Latest epoch of the departure window
    earliest_arrival_time: time_conversion.DateTime
        Earliest epoch of the arrival window
    latest_arrival_time: time_conversion.DateTime
        Latest epoch of the arrival window
    time_resolution: float
        Resolution used to discretize the departure/arrival time windows
    function_to_calculate_delta_v: callable = calculate_lambert_arc_impulsive_delta_v
        Function with which the manoeuvre's required ΔV will be calculated
    C3: bool = False
        Whether to plot C3 (specific energy) instead of ΔV
    total: bool = False
        Whether to plot departure and arrival ΔV/C3, or only the total ΔV/C3. This option is only respected if the ΔV map obtained from
    threshold: float = 10
        Upper threshold beyond which ΔV/C3 is not plotted. This is useful to mask regions of the plot where the ΔV/C3 is too high to be of interest.
    upscale: bool = False
        Whether to use interpolation to increase the resolution of the plot. This is not always reliable, and the detail generated cannot be relied upon for analysis. Its only purpose is aesthetic improvement.
    number_of_levels: int = 10
        The number of levels in the ΔV/C3 contour plot
    percent_margin: float = 5
        Empty margin between the axes of the plot and the plotted data
    figsize: tuple[int, int] = (8, 8)
        Size of the figure
    show: bool bool = True
        Whether to show the plot        
    save: bool = False
        Whether to save the plot
    filename: str = 'porkchop.png'
        The filename used for the saved plot
    
    Output
    ------
    departure_epochs: np.ndarray
        Discretized departure time window
    arrival_epochs: np.ndarray
        Discretized arrival time window
    ΔV: np.ndarray
        Array containing the ΔV of all coordinates of the grid of departure/arrival epochs
    """

    # Calculate ΔV map
    [departure_epochs, arrival_epochs, ΔV] = calculate_delta_v_time_map(
        bodies                        = bodies,
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