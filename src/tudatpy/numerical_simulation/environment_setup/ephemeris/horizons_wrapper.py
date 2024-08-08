import datetime
from typing import Union

from ....data.horizons import HorizonsQuery


def jpl_horizons(
    horizons_query: str,
    horizons_location: str,
    frame_origin: str,
    frame_orientation: str = "ECLIPJ2000",
    query_type: str = "default",
    epoch_start: Union[datetime.datetime, float, None] = None,
    epoch_end: Union[datetime.datetime, float, None] = None,
    epoch_step: Union[str, None] = None,
    epoch_list: Union[list, None] = None,
    extended_query: bool = False,
    aberations: str = "geometric",
):
    """Factory function for creating ephemeris model settings from tabulated JPL Horizons vectors.

    JPL Horizons provides access to highly accurate ephemerides for many solar system objects,
    including asteriuds, comets, planets, moons and select spacecraft.

    This function is a wrapper for the tudatpy.data.horizons functionality.
    That api is not available on the api documentation yet.
    For now, visit the HorizonsQuery souce code for extensive documentation:
    https://github.com/tudat-team/tudatpy/blob/master/tudatpy/data/horizons.py

    For more information on the Horizons System, visit: https://ssd.jpl.nasa.gov/horizons/manual.html

    Examples
    ----------
    Add Ephemerides of JUICE to the body_settings

        >>> juice_eph_settings = jpl_horizons(
                horizons_query="-121", #-121 is the query code for JUICE
                horizons_location="500@399", # Geocentre@Earth
                frame_origin="Earth", #tudat frame origin and orientation
                frame_orientation="ECLIPJ2000",
                epoch_start=datetime.datetime(2018, 10, 21),
                epoch_end=datetime.datetime(2023, 9, 1),
                epoch_step="1d",
                extended_query=True,
            )

        >>> body_settings.add_empty_settings("JUICE")
        >>> body_settings.get("JUICE").ephemeris_settings = juice_eph_settings
    """
    query = HorizonsQuery(
        query_id=horizons_query,
        location=horizons_location,
        query_type=query_type,
        epoch_start=epoch_start,
        epoch_end=epoch_end,
        epoch_step=epoch_step,
        epoch_list=epoch_list,
        extended_query=extended_query,
    )

    eph = query.create_ephemeris_tabulated(
        frame_origin=frame_origin,
        frame_orientation=frame_orientation,
        aberations=aberations,
    )

    return eph
