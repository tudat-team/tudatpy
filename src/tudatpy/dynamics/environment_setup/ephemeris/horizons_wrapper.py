import datetime
import warnings
from typing import Union

from tudatpy.data_input.environment_data.horizons import (
    HorizonsBatch as _HorizonsBatch,
    HorizonsQuery as _HorizonsQuery,
)


class HorizonsQuery(_HorizonsQuery):
    def __init__(self, *args, **kwargs):
        warnings.warn(
            "`tudatpy.dynamics.environment_setup.ephemeris.HorizonsQuery` is deprecated. "
            "Use `tudatpy.data_input.environment_data.horizons.HorizonsQuery` instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(*args, **kwargs)

    def create_ephemeris_tabulated(
        self,
        frame_origin: str,
        frame_orientation: str = "ECLIPJ2000",
        aberations: str = "geometric",
    ):
        warnings.warn(
            "`HorizonsQuery.create_ephemeris_tabulated` is deprecated. "
            "Use `tudatpy.dynamics.environment_setup.ephemeris.jpl_horizons_from_query` instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return jpl_horizons_from_query(
            self,
            frame_origin=frame_origin,
            frame_orientation=frame_orientation,
            aberations=aberations,
        )


class HorizonsBatch(_HorizonsBatch):
    def __init__(self, *args, **kwargs):
        warnings.warn(
            "`tudatpy.dynamics.environment_setup.ephemeris.HorizonsBatch` is deprecated. "
            "Use `tudatpy.data_input.environment_data.horizons.HorizonsBatch` instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        super().__init__(*args, **kwargs)

    def add_batch_ephemerides(
        self,
        body_settings,
        frame_origin: str,
        frame_orientation: str = "ECLIPJ2000",
        aberations: str = "geometric",
    ) -> None:
        warnings.warn(
            "`HorizonsBatch.add_batch_ephemerides` is deprecated. "
            "Use `tudatpy.dynamics.environment_setup.ephemeris.add_horizons_batch_ephemerides` instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return add_horizons_batch_ephemerides(
            self,
            body_settings=body_settings,
            frame_origin=frame_origin,
            frame_orientation=frame_orientation,
            aberations=aberations,
        )


def jpl_horizons_from_query(
    query: HorizonsQuery,
    frame_origin: str,
    frame_orientation: str = "ECLIPJ2000",
    aberations: str = "geometric",
):
    """Create tabulated ephemeris settings from an existing Horizons query.

    This function is equivalent to :func:`jpl_horizons`, but takes an existing
    :class:`~tudatpy.data_input.environment_data.horizons.HorizonsQuery`
    object instead of creating the query internally. Use it when the query needs
    to be inspected or modified before converting it to ephemeris settings.
    """
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
        eph = jpl_horizons_from_query(
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
    query = _HorizonsQuery(
        query_id=horizons_query,
        location=horizons_location,
        query_type=query_type,
        epoch_start=epoch_start,
        epoch_end=epoch_end,
        epoch_step=epoch_step,
        epoch_list=epoch_list,
        extended_query=extended_query,
    )

    return jpl_horizons_from_query(
        query,
        frame_origin=frame_origin,
        frame_orientation=frame_orientation,
        aberations=aberations,
    )
