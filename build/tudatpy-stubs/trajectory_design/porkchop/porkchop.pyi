import typing
from tudatpy import constants as constants
from tudatpy.astro import time_representation as time_representation
from tudatpy.dynamics import environment as environment
from tudatpy.trajectory_design.porkchop._lambert import calculate_lambert_arc_impulsive_delta_v as calculate_lambert_arc_impulsive_delta_v
from tudatpy.trajectory_design.porkchop._plot_porkchop import plot_porkchop as plot_porkchop

def determine_shape_of_delta_v(bodies: environment.SystemOfBodies, departure_body: str, target_body: str, departure_epoch: float, arrival_epoch: float, function_to_calculate_delta_v: callable=...):
    """Determine whether `function_to_calculate_delta_v` returns ΔV as
    
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
        2 if it is a list/tuple/NumPy array containing [departure ΔV, arrival ΔV]"""

def calculate_delta_v_time_map(bodies: environment.SystemOfBodies, departure_body: str, target_body: str, earliest_departure_time: time_representation.DateTime, latest_departure_time: time_representation.DateTime, earliest_arrival_time: time_representation.DateTime, latest_arrival_time: time_representation.DateTime, time_resolution: float, function_to_calculate_delta_v: callable=...):
    """Creates an array containing the ΔV of all coordinates of the grid of departure/arrival epochs.
    
    Parameters
    ----------
    bodies: environment.SystemOfBodies
        Body objects defining the physical simulation environment
    departure_body: str
        The name of the body from which the transfer is to be computed
    target_body: str
        The name of the body to which the transfer is to be computed
    earliest_departure_time: time_representation.DateTime
        Earliest epoch of the departure window
    latest_departure_time: time_representation.DateTime
        Latest epoch of the departure window
    earliest_arrival_time: time_representation.DateTime
        Earliest epoch of the arrival window
    latest_arrival_time: time_representation.DateTime
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
        Array containing the ΔV of all coordinates of the grid of departure/arrival epochs"""

def porkchop(bodies: environment.SystemOfBodies, departure_body: str, target_body: str, earliest_departure_time: time_representation.DateTime, latest_departure_time: time_representation.DateTime, earliest_arrival_time: time_representation.DateTime, latest_arrival_time: time_representation.DateTime, time_resolution: float, function_to_calculate_delta_v: callable=..., C3: bool=False, total: bool=False, threshold: float=10, upscale: bool=False, number_of_levels: int=10, percent_margin: float=5, figsize: tuple[int, int]=(8, 8), show: bool=True, save: bool=False, filename: str='porkchop.png') -> None:
    """Calculates and displays ΔV/C3 porkchop mission design plots.
    
    Parameters
    ----------
    bodies: environment.SystemOfBodies
        Body objects defining the physical simulation environment
    departure_body: str
        The name of the body from which the transfer is to be computed
    target_body: str
        The name of the body to which the transfer is to be computed
    earliest_departure_time: time_representation.DateTime
        Earliest epoch of the departure window
    latest_departure_time: time_representation.DateTime
        Latest epoch of the departure window
    earliest_arrival_time: time_representation.DateTime
        Earliest epoch of the arrival window
    latest_arrival_time: time_representation.DateTime
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
        Array containing the ΔV of all coordinates of the grid of departure/arrival epochs"""