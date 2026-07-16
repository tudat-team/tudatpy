"""Radar ground-station helpers.

The station catalogue in this module is intentionally small and explicit. It
contains radar stations for which the current TudatPy radar readers can assign
known geodetic positions without asking the user for additional station
metadata.
"""

from collections.abc import Iterable
from typing import Any

import numpy as np

from tudatpy.astro import element_conversion
from tudatpy.dynamics import environment_setup

JPL_TO_MPC_RADAR_STATIONS = {
    "-1": "251",
    "-14": "253",
}

_POSITIONS_BY_MPC_CODE = {
    "251": np.array([453.34, np.deg2rad(18.3442199), np.deg2rad(293.2473068)]),
    "253": np.array([1001.39, np.deg2rad(35.4259009), np.deg2rad(243.1104618)]),
}

RADAR_STATION_GEODETIC = {
    **_POSITIONS_BY_MPC_CODE,
    # Keep source-prefixed aliases in the catalogue so data readers can retain
    # source-native station IDs without requiring users to remap stations.
    **{f"MPC:{code}": position for code, position in _POSITIONS_BY_MPC_CODE.items()},
    **{
        f"JPL:{jpl_code}": _POSITIONS_BY_MPC_CODE[mpc_code]
        for jpl_code, mpc_code in JPL_TO_MPC_RADAR_STATIONS.items()
    },
}


def normalize_radar_station_id(source: str, raw_station_id: str | int) -> str:
    """Return the canonical station ID used by the radar interfaces.

    Parameters
    ----------
    source : str
        Source catalogue name. MPC station IDs are normalized to three
        characters; other sources are prefixed as ``SOURCE:id``.
    raw_station_id : str | int
        Station identifier as provided by the source catalogue.

    Returns
    -------
    str
        Canonical station identifier used in TudatPy radar data tables.
    """
    source = str(source).strip().upper()
    station_id = str(raw_station_id).strip()
    if ":" in station_id:
        return station_id
    if source == "MPC":
        return station_id.zfill(3)
    return f"{source}:{station_id}"


def get_radar_station_geodetic_positions(include_aliases: bool = True) -> dict[str, np.ndarray]:
    """Return known radar station geodetic positions as copied arrays.

    Parameters
    ----------
    include_aliases : bool, default True
        If True, include MPC and source-prefixed aliases such as ``MPC:253``
        and ``JPL:-14``. If False, return only the MPC station-code keys.

    Returns
    -------
    dict[str, numpy.ndarray]
        Mapping from station ID to copied geodetic position vector.
    """
    positions = RADAR_STATION_GEODETIC if include_aliases else _POSITIONS_BY_MPC_CODE
    return {station: position.copy() for station, position in positions.items()}


def radar_ground_station_settings(
    station_ids: Iterable[str | int] | None = None,
    station_positions: dict[str, np.ndarray] | None = None,
    include_aliases: bool = True,
):
    """Create ground-station settings for known radar stations.

    Parameters
    ----------
    station_ids : iterable[str | int] | None, default None
        Radar station IDs to include. If None, settings are created for all
        known radar stations.
    station_positions : dict[str, numpy.ndarray] | None, default None
        Optional geodetic position overrides or extensions.
    include_aliases : bool, default True
        Whether the built-in catalogue includes aliased station IDs when
        ``station_ids`` is None.

    Returns
    -------
    list
        Ground-station settings objects for use with Tudat environment setup.

    Raises
    ------
    KeyError
        If a requested station has no known or user-provided geodetic position.
    """
    positions = get_radar_station_geodetic_positions(include_aliases=include_aliases)
    if station_positions is not None:
        positions.update(
            {str(station): np.asarray(position) for station, position in station_positions.items()}
        )

    station_ids = (
        sorted(positions)
        if station_ids is None
        else sorted({str(station) for station in station_ids})
    )
    settings = []
    for station_id in station_ids:
        if station_id not in positions:
            raise KeyError(
                f"No geodetic radar station position is available for '{station_id}'. "
                "Pass station_positions or add the station to the environment manually."
            )
        settings.append(
            environment_setup.ground_station.basic_station(
                station_id,
                positions[station_id],
                element_conversion.geodetic_position_type,
            )
        )
    return settings


def add_radar_ground_stations(
    bodies: Any,
    station_ids: Iterable[str | int],
    station_body: str = "Earth",
    station_positions: dict[str, np.ndarray] | None = None,
) -> None:
    """Add selected known radar ground stations to an existing Tudat body.

    Parameters
    ----------
    bodies : SystemOfBodies
        Tudat system of bodies to which stations are added.
    station_ids : iterable[str | int]
        Radar station IDs to add.
    station_body : str, default "Earth"
        Body on which the station reference points are defined.
    station_positions : dict[str, numpy.ndarray] | None, default None
        Optional geodetic position overrides or extensions.

    Returns
    -------
    None
    """
    station_ids = sorted({str(station) for station in station_ids})
    settings = radar_ground_station_settings(
        station_ids=station_ids,
        station_positions=station_positions,
    )
    for station_name, station_settings in zip(station_ids, settings):
        if station_name in bodies.get(station_body).ground_station_list:
            continue
        environment_setup.add_ground_station(bodies.get_body(station_body), station_settings)


def add_all_radar_ground_stations(
    bodies: Any,
    station_body: str = "Earth",
    station_positions: dict[str, np.ndarray] | None = None,
    include_aliases: bool = True,
) -> None:
    """Add all known radar ground stations to an existing Tudat body.

    Parameters
    ----------
    bodies : SystemOfBodies
        Tudat system of bodies to which stations are added.
    station_body : str, default "Earth"
        Body on which the station reference points are defined.
    station_positions : dict[str, numpy.ndarray] | None, default None
        Optional geodetic position overrides or extensions.
    include_aliases : bool, default True
        Whether source-prefixed aliases are added as separate station reference
        points.

    Returns
    -------
    None
    """
    station_ids = get_radar_station_geodetic_positions(include_aliases=include_aliases).keys()
    add_radar_ground_stations(
        bodies,
        station_ids=station_ids,
        station_body=station_body,
        station_positions=station_positions,
    )
