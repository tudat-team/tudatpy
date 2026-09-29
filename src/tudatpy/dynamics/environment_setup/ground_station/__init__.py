from tudatpy.kernel.dynamics.environment_setup.ground_station import *

from ._station_positions import JPL_RADAR_STATION_POSITIONS

_mpc_station_settings = optical_telescope_stations


def jpl_radar_stations():
    """Return settings for the stations in the JPL radar catalog."""
    return [
        basic_station(f"JPL:{station_code}", position)
        for station_code, position in JPL_RADAR_STATION_POSITIONS.items()
    ]


def optical_telescope_stations():
    """Return settings for the MPC and JPL small-body tracking stations."""
    return _mpc_station_settings() + jpl_radar_stations()
