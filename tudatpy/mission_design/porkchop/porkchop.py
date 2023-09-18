# General imports
import numpy as np
from tqdm import tqdm

# Tudat imports
from tudatpy.kernel import constants
from tudatpy.kernel.astro import time_conversion
from tudatpy.kernel.numerical_simulation import environment

# Custom imports
from tudatpy.mission_design.porkchop import plot_porkchop
from tudatpy.budgeting.delta_v.lambert import calculate_lambert_arc_impulsive_delta_v


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
        function_to_calculate_delta_v: callable = calculate_lambert_arc_impulsive_delta_v
    ):

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

    # Pre-allocate ΔV and C3 array
    ΔV = np.full((n_dep, n_arr, 2), np.nan)

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
        function_to_calculate_delta_v: callable = calculate_lambert_arc_impulsive_delta_v,
        # Plot arguments
        C3: bool = False,
        total: bool = False,
        threshold: float = 10,
        upscale: bool = True,
        number_of_levels: int = 10,
        # Figure arguments
        percent_margin: float = 5,
        figsize: tuple[int, int] = (8, 8),
        show: bool = True,
        save: bool = False,
        filename: str = 'porkchop.png',
    ) -> None:
    """
    Create a porkchop plot from numerical propagations.
    """

    # Calculate ΔV map
    [departure_epochs, arrival_epochs, ΔV] = calculate_delta_v_time_map(
        bodies                   = bodies,
        global_frame_orientation = global_frame_orientation,
        departure_body           = departure_body,
        target_body              = target_body,
        earliest_departure_time  = earliest_departure_time,
        latest_departure_time    = latest_departure_time,
        earliest_arrival_time    = earliest_arrival_time,
        latest_arrival_time      = latest_arrival_time,
        time_resolution          = time_resolution,
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