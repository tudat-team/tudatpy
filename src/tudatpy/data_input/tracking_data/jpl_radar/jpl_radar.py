"""JPL Small-Body Radar API retrieval."""

from __future__ import annotations

import datetime
import json
import urllib.parse
import urllib.request
from typing import Any

import numpy as np
import pandas as pd

from tudatpy.astro.time_representation import DateTime
from tudatpy.constants import SPEED_OF_LIGHT
from tudatpy.data_input.tracking_data.radar_utilities import (
    DOPPLER_OBSERVABLE,
    RANGE_OBSERVABLE,
    empty_radar_table,
    validate_radar_data,
)
from tudatpy.data_input.tracking_data.radar_utilities.stations import (
    JPL_TO_MPC_RADAR_STATIONS,
    normalize_radar_station_id,
)

_API_URL = "https://ssd-api.jpl.nasa.gov/sb_radar.api"


def _validate_response(response: Any, context: str) -> dict[str, Any]:
    """Validate the decoded JPL API response before row extraction."""
    if not isinstance(response, dict):
        raise RuntimeError(
            f"JPL radar API returned an unexpected response for {context}: "
            f"expected a JSON object, got {type(response).__name__}."
        )

    validated_response = dict(response)
    for key in ["fields", "data"]:
        value = validated_response.get(key, [])
        if value is None:
            value = []
        if not isinstance(value, list):
            raise RuntimeError(
                f"JPL radar API returned invalid '{key}' for {context}: "
                f"expected a list, got {type(value).__name__}."
            )
        validated_response[key] = value
    return validated_response


def _fetch_json(query_parameters: dict[str, str] | None, timeout: float, context: str) -> dict:
    """Fetch a JPL radar API JSON response and validate its top-level shape."""
    query = urllib.parse.urlencode(query_parameters or {})
    url = f"{_API_URL}?{query}" if query else _API_URL
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return _validate_response(json.load(response), context)


def get_available_radar_targets(timeout: float = 30.0) -> list[str]:
    """Return JPL designations for small bodies with radar data.

    Parameters
    ----------
    timeout : float, default 30.0
        Network timeout, in seconds, used for the JPL API request.

    Returns
    -------
    list[str]
        Sorted list of unique JPL small-body designations for which the API
        reports radar data.

    Raises
    ------
    RuntimeError
        If the API response does not contain the expected designation field.
    """
    response = _fetch_json(None, timeout, "all radar targets")
    fields = response.get("fields", [])
    if "des" not in fields:
        raise RuntimeError("JPL radar API response does not contain a 'des' field.")

    designation_index = fields.index("des")
    designations = {
        str(row[designation_index]).strip()
        for row in response.get("data", [])
        if len(row) > designation_index and str(row[designation_index]).strip()
    }
    return sorted(designations)


class JPLRadarQuery:
    """Query JPL small-body radar data and normalize it to TudatPy tables."""

    API_URL = _API_URL

    def __init__(self, target: str | int, timeout: float = 30.0) -> None:
        """Create a query for a JPL small-body target designation.

        Parameters
        ----------
        target : str | int
            JPL small-body designation.
        timeout : float, default 30.0
            Network timeout, in seconds, used for the JPL API request.
        """
        self.target = str(target)
        self.timeout = timeout
        self._response_cache = None

    def _validate_response(self, response: Any) -> dict[str, Any]:
        """Validate a decoded response in the context of this query target."""
        return _validate_response(response, f"target {self.target}")

    def _fetch_json(self) -> dict[str, Any]:
        """Fetch and cache the raw JPL radar API JSON response."""
        if self._response_cache is not None:
            return self._validate_response(self._response_cache)

        self._response_cache = _fetch_json(
            {"des": self.target},
            self.timeout,
            f"target {self.target}",
        )
        return self._response_cache

    def to_radar_data(
        self,
        target_body: str | None = None,
        epoch_start: datetime.datetime | None = None,
        epoch_end: datetime.datetime | None = None,
        target_point: str | None = "C",
        include_range: bool = True,
        include_doppler: bool = True,
        station_id_mode: str = "jpl",
    ) -> pd.DataFrame:
        """Return JPL radar data as a canonical radar-data table.

        JPL range rows in microseconds are converted to metres. Doppler rows
        are converted from frequency shift to measured received frequency by
        adding the transmitter frequency. Epochs are retained in UTC.

        Parameters
        ----------
        target_body : str | None, default None
            Target body name assigned to the returned radar rows. If None, the
            JPL query target is used.
        epoch_start : datetime.datetime | None, default None
            Optional inclusive lower bound on the UTC API epoch.
        epoch_end : datetime.datetime | None, default None
            Optional inclusive upper bound on the UTC API epoch.
        target_point : str | None, default "C"
            Body point code to retain. If None, all target points are retained.
        include_range : bool, default True
            Whether to include JPL delay rows.
        include_doppler : bool, default True
            Whether to include JPL Doppler rows.
        station_id_mode : str, default "jpl"
            Station ID convention. ``"jpl"`` preserves source-prefixed JPL
            station IDs; ``"mpc"`` maps known JPL station IDs to MPC radar
            station codes.

        Returns
        -------
        pandas.DataFrame
            Canonical radar data table. If no selected rows are available, an
            empty canonical radar table is returned.
        """
        response = self._fetch_json()
        fields = response.get("fields", [])
        rows = []
        target_body = str(target_body or self.target)

        for raw_row in response.get("data", []):
            row = dict(zip(fields, raw_row))
            if target_point is not None and row.get("bp") != target_point:
                continue

            # DateTime is used below only to preserve Tudat's UTC seconds since
            # J2000 convention. No UTC/TDB conversion is done in data input.
            epoch_datetime = datetime.datetime.strptime(row["epoch"], "%Y-%m-%d %H:%M:%S")
            if epoch_start is not None and epoch_datetime < epoch_start:
                continue
            if epoch_end is not None and epoch_datetime > epoch_end:
                continue

            units = str(row.get("units", "")).strip()
            if units == "us":
                if not include_range:
                    continue
                observable_type = RANGE_OBSERVABLE
            elif units == "Hz":
                if not include_doppler:
                    continue
                observable_type = DOPPLER_OBSERVABLE
            else:
                raise RuntimeError(
                    f"Unsupported JPL radar measurement unit '{units}' for target "
                    f"{self.target} at epoch {row.get('epoch')}. Expected 'us' or 'Hz'."
                )

            date_time = DateTime.from_iso_string(row["epoch"])
            epoch_utc = date_time.to_epoch()

            measurement = float(row["value"])
            sigma = float(row["sigma"])
            frequency_mhz = float(row["freq"])
            transmitter = self._station_id(row["xmit"], station_id_mode)
            receiver = self._station_id(row["rcvr"], station_id_mode)
            radar_row = {
                "target_body": target_body,
                "epoch_seconds_UTC": epoch_utc,
                "transmitter": transmitter,
                "receiver": receiver,
                "target_point": row.get("bp"),
                "transmitter_frequency_hz": np.nan,
                "source": "JPL",
                "raw_units": units,
                "jpl_transmitter": str(row["xmit"]).strip(),
                "jpl_receiver": str(row["rcvr"]).strip(),
            }

            if observable_type == RANGE_OBSERVABLE:
                # JPL delay values are given in microseconds; Tudat range
                # TrackingData expects metres.
                radar_row.update(
                    {
                        "observable_type": RANGE_OBSERVABLE,
                        "value": SPEED_OF_LIGHT * measurement * 1.0e-6,
                        "sigma": SPEED_OF_LIGHT * sigma * 1.0e-6,
                        "raw_delay_us": measurement,
                        "raw_delay_sigma_us": sigma,
                    }
                )
            else:
                transmitter_frequency_hz = frequency_mhz * 1.0e6
                # JPL Doppler values are frequency shifts. The Tudat observable
                # is the measured received frequency, so the shift is added to
                # the transmitter frequency.
                radar_row.update(
                    {
                        "observable_type": DOPPLER_OBSERVABLE,
                        "value": transmitter_frequency_hz + measurement,
                        "sigma": sigma,
                        "transmitter_frequency_hz": transmitter_frequency_hz,
                        "raw_doppler_shift_hz": measurement,
                    }
                )
            rows.append(radar_row)

        if not rows:
            return empty_radar_table()
        return validate_radar_data(pd.DataFrame(rows))

    @staticmethod
    def _station_id(raw_station_id: str | int, station_id_mode: str) -> str:
        """Normalize a JPL station ID or map it to MPC compatibility mode."""
        station_id_mode = station_id_mode.lower()
        raw_station_id = str(raw_station_id).strip()
        if station_id_mode == "jpl":
            return normalize_radar_station_id("JPL", raw_station_id)
        if station_id_mode == "mpc":
            if raw_station_id not in JPL_TO_MPC_RADAR_STATIONS:
                raise ValueError(
                    "No MPC radar station mapping is defined for JPL station " f"{raw_station_id}."
                )
            return JPL_TO_MPC_RADAR_STATIONS[raw_station_id]
        raise ValueError("station_id_mode must be 'jpl' or 'mpc'.")
