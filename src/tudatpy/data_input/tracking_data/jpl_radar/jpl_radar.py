"""JPL Small-Body Radar Astrometry API (https://ssd-api.jpl.nasa.gov/doc/sb_radar.html)."""

import functools

import numpy as np
import pandas as pd
import requests
from astropy import units as u

from tudatpy.astro import time_representation
from tudatpy.data_input.tracking_data.radar_utilities import (
    empty_radar_table,
    filter_radar_data,
    radar_data_from_raw,
    radar_data_to_tracking_data,
)

_API_URL = "https://ssd-api.jpl.nasa.gov/sb_radar.api"
_JPL_TO_MPC_STATION = {
    "-1": "251",
    "-2": "254",
    "-9": "256",
    "-13": "252",
    "-14": "253",
    "-25": "257",
    "-38": "255",
    "-73": "259",
}


def _query(params: dict, timeout: float) -> dict:
    response = requests.get(_API_URL, params=params, timeout=timeout)
    response.raise_for_status()
    return response.json()


def _as_frame(content: dict) -> pd.DataFrame:
    return pd.DataFrame(content.get("data", []), columns=content.get("fields"))


def _station_id(jpl_code) -> str:
    jpl_code = str(jpl_code).strip()
    return _JPL_TO_MPC_STATION.get(jpl_code, f"JPL:{jpl_code}")


def get_available_radar_targets(timeout: float = 30.0) -> list[str]:
    """Return the JPL designations of all small bodies with radar astrometry.

    Parameters
    ----------
    timeout : float, default 30.0
        Request timeout [s].

    Returns
    -------
    list[str]
        Sorted, unique designations.

    Raises
    ------
    requests.HTTPError
        If the JPL API request is unsuccessful.
    """
    data = _as_frame(_query({}, timeout))
    return sorted(data["des"].str.strip().unique()) if not data.empty else []


class JPLRadarQuery:
    """Radar astrometry for one small body from the JPL Small-Body Radar API.

    Parameters
    ----------
    target : str | int
        JPL small-body designation, for example ``433`` or ``"2004 VB"``.
    timeout : float, default 30.0
        Request timeout [s].
    """

    def __init__(self, target: str | int, timeout: float = 30.0) -> None:
        self.target = str(target)
        self.timeout = timeout

    @functools.cached_property
    def _content(self) -> dict:
        return _query({"des": self.target, "coords": 1}, self.timeout)

    @property
    def raw_data(self) -> pd.DataFrame:
        """**read-only**

        Records returned by the JPL API, with one row per measurement.

        :type: pandas.DataFrame
        """
        return _as_frame(self._content)

    def station_geodetic_positions(self) -> dict[str, np.ndarray]:
        """Return the geodetic positions of stations used by the observations.

        Returns
        -------
        dict[str, numpy.ndarray]
            MPC station codes mapped to arrays containing altitude [m],
            latitude [rad] and longitude [rad], in that order. Unmapped JPL
            stations use an identifier of the form ``"JPL:<code>"``.

        Raises
        ------
        requests.HTTPError
            If the JPL API request is unsuccessful.
        """
        return {
            _station_id(code): np.array(
                [
                    (float(station["altitude"]) * u.Unit(station["alt_units"])).to_value(u.m),
                    np.deg2rad(float(station["latitude"])),
                    np.deg2rad(float(station["longitude"])),
                ]
            )
            for code, station in self._content.get("coords", {}).items()
        }

    def to_radar_data(
        self,
        target_body=None,
        epoch_start=None,
        epoch_end=None,
        target_point="C",
    ) -> pd.DataFrame:
        """Return the measurements as a canonical radar table.

        Delays [us] become round-trip ranges [m]; Doppler shifts [Hz] become
        received frequencies. Known JPL station codes are mapped to MPC codes.

        Parameters
        ----------
        target_body : str | None, default None
            Tudat body name for the target; defaults to the query designation.
        epoch_start, epoch_end : float | datetime.datetime | None, default None
            Inclusive epoch bounds in UTC seconds since J2000 or datetime-like.
        target_point : str | None, default "C"
            Bounce point to keep; None keeps all points.

        Returns
        -------
        pandas.DataFrame
            Canonical radar table with one row per delay or Doppler measurement.

        Raises
        ------
        requests.HTTPError
            If the JPL API request is unsuccessful.
        RuntimeError
            If the API returns a radar measurement with unsupported units.
        """
        raw = self.raw_data
        if raw.empty:
            return empty_radar_table()
        unknown_units = set(raw["units"]) - {"us", "Hz"}
        if unknown_units:
            raise RuntimeError(
                f"Unsupported JPL radar units for {self.target}: {sorted(unknown_units)}"
            )
        value = raw["value"].astype(float)
        sigma = raw["sigma"].astype(float)
        is_delay = raw["units"] == "us"
        table = radar_data_from_raw(
            pd.DataFrame(
                {
                    "target_body": str(target_body or self.target),
                    "epoch_seconds_UTC": [
                        time_representation.iso_string_to_epoch(epoch) for epoch in raw["epoch"]
                    ],
                    "transmitter": raw["xmit"].map(_station_id),
                    "receiver": raw["rcvr"].map(_station_id),
                    "target_point": raw["bp"],
                    "transmitter_frequency_hz": raw["freq"].astype(float) * 1.0e6,
                    "delay_us": value.where(is_delay),
                    "delay_sigma_us": sigma.where(is_delay),
                    "doppler_hz": value.mask(is_delay),
                    "doppler_sigma_hz": sigma.mask(is_delay),
                }
            ),
            source="JPL",
        )
        return filter_radar_data(
            table,
            epoch_start=epoch_start,
            epoch_end=epoch_end,
            target_point=target_point,
        )


def read_jpl_radar_data(
    target: str | int,
    target_body=None,
    epoch_start=None,
    epoch_end=None,
    timeout: float = 30.0,
):
    """Retrieve JPL radar astrometry and return Tudat tracking data.

    Parameters
    ----------
    target : str | int
        JPL small-body designation, for example ``433`` or ``"2004 VB"``.
    target_body : str | None, default None
        Tudat body name for the target; defaults to the query designation.
    epoch_start : float | DateTime | Time | datetime.datetime | None, default None
        Inclusive lower epoch bound. Numeric values are UTC seconds since J2000.
    epoch_end : float | DateTime | Time | datetime.datetime | None, default None
        Inclusive upper epoch bound. Numeric values are UTC seconds since J2000.
    timeout : float, default 30.0
        HTTP request timeout [s].

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        Tracking data and transmitter-frequency supplementary data.

    Raises
    ------
    requests.HTTPError
        If the JPL API request is unsuccessful.
    RuntimeError
        If the API returns a radar measurement with unsupported units.
    """
    table = JPLRadarQuery(target, timeout).to_radar_data(
        target_body=target_body,
        epoch_start=epoch_start,
        epoch_end=epoch_end,
    )
    return radar_data_to_tracking_data(table)
