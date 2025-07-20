import numpy
import pybind11_stubgen.typing_ext
from ....dynamics import environment
import typing
__all__ = ['compute_target_angles_and_range', 'compute_target_angles_and_range_vectors']

def compute_target_angles_and_range(bodies: environment.SystemOfBodies, station_id: tuple[str, str], target_body: str, observation_times: list[float], is_station_transmitting: bool) -> dict[float, typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]:
    """Function to compute the azimuth angle, elevation angle and range at a ground station.
    
    Function to compute the azimuth angle, elevation angle and range at a ground station. This functions is provided as a function of
    convenience, to prevent users having to manually define the relevant settings for this often-needed functionality. This function
    takes an observing station and a target body as input, and provides the observed angles and current range (without correction for aberrations, with correction for light time)
    as observed at that station
    
    
    Parameters
    ----------
    bodies : SystemOfBodies
        System of bodies that defines the full physical environment
    
    station_id : tuple[ str, str]
        Identifier for the observing station, as a pair of strings: the body name and the station name.
    
    target_body : str
        Name of body which is observed by ground station
    
    observation_times : list[float]
        List of times at which the ground station observations are to be analyzed
    
    is_station_transmitting : bool
        Boolean defining whether the observation times define times at which the station is transmitting to, or receiving from, the ground station.
        This has an impact on the whether the light-time is computed forward or backward in time from the ground station to the target
    
    Returns
    -------
    dict[float,numpy.ndarray[numpy.float64[3, 1]]]
        Dictionary with the required output. Key defines the observation time, the value is an array of size three containing entry 0 - elevation angle, entry 1 - azimuth angle, entry 2 - range"""

def compute_target_angles_and_range_vectors(bodies: environment.SystemOfBodies, station_id: tuple[str, str], target_body: str, observation_times: list[float], is_station_transmitting: bool) -> tuple[list[float], list[typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 1)]]]:
    """Function to compute the azimuth angle, elevation angle and range at a ground station.
    
    Function to compute the azimuth angle, elevation angle and range at a ground station. This functions is provided as a function of
    convenience, to prevent users having to manually define the relevant settings for this often-needed functionality. This function
    takes an observing station and a target body as input, and provides the observed angles and current range (without correction for aberrations, with correction for light time)
    as observed at that station
    
    
    Parameters
    ----------
    bodies : SystemOfBodies
        System of bodies that defines the full physical environment
    
    station_id : tuple[ str, str]
        Identifier for the observing station, as a pair of strings: the body name and the station name.
    
    target_body : str
        Name of body which is observed by ground station
    
    observation_times : list[float]
        List of times at which the ground station observations are to be analyzed
    
    is_station_transmitting : bool
        Boolean defining whether the observation times define times at which the station is transmitting to, or receiving from, the ground station.
        This has an impact on the whether the light-time is computed forward or backward in time from the ground station to the target
    
    Returns
    -------
    dict[float,numpy.ndarray[numpy.float64[3, 1]]]
        Dictionary with the required output. Key defines the observation time, the value is an array of size three containing entry 0 - elevation angle, entry 1 - azimuth angle, entry 2 - range"""