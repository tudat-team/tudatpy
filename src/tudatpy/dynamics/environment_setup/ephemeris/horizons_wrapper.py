import datetime
from typing import Union

from tudatpy.data_input.environment_data.horizons import HorizonsBatch, HorizonsQuery


def horizons_query_to_ephemeris_settings(
    query: HorizonsQuery,
    frame_origin: str,
    frame_orientation: str = "ECLIPJ2000",
    aberations: str = "geometric",
):
    """Create tabulated ephemeris settings from a Horizons query."""
    if frame_orientation not in ["ECLIPJ2000", "J2000"]:
        raise ValueError("refplane parameter must be one of: " + '"ECLIPJ2000", "J2000"')

    vector = query.cartesian(frame_orientation=frame_orientation, aberations=aberations)
    table = {x[0]: x[1:7] for x in vector}

    from tudatpy.dynamics.environment_setup import ephemeris

    return ephemeris.tabulated(
        body_state_history=table,
        frame_origin=frame_origin,
        frame_orientation=frame_orientation,
    )


def add_horizons_batch_ephemerides(
    batch: HorizonsBatch,
    body_settings,
    frame_origin: str,
    frame_orientation: str = "ECLIPJ2000",
    aberations: str = "geometric",
) -> None:
    """Add tabulated Horizons ephemerides for all queries in a batch."""
    names = []

    for query in batch._query_objects.values():
        eph = horizons_query_to_ephemeris_settings(
            query,
            frame_origin=frame_origin,
            frame_orientation=frame_orientation,
            aberations=aberations,
        )
        name = query.name
        names.append(name)

        body_settings.add_empty_settings(name)
        body_settings.get(name).ephemeris_settings = eph

    batch._names = names


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
    """Factory function for ephemeris model settings from JPL Horizons vectors."""
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

    return horizons_query_to_ephemeris_settings(
        query,
        frame_origin=frame_origin,
        frame_orientation=frame_orientation,
        aberations=aberations,
    )
