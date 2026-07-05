from __future__ import annotations

import datetime
import json
import urllib.parse
import urllib.request
from typing import Any

import numpy as np
import pandas as pd

from tudatpy.astro import time_representation
from tudatpy.astro.time_representation import DateTime
from tudatpy.constants import SPEED_OF_LIGHT
from tudatpy.data.radar import (
    DOPPLER_OBSERVABLE,
    RANGE_OBSERVABLE,
    RadarTrackingData,
    empty_radar_table,
)
from tudatpy.data.radar.stations import JPL_TO_MPC_RADAR_STATIONS, normalize_radar_station_id


class JPLRadarQuery:
    """Query and convert JPL small-body radar tracking data."""

    API_URL = "https://ssd-api.jpl.nasa.gov/sb_radar.api"

    def __init__(self, target: str | int, timeout: float = 30.0) -> None:
        """Create a query for a JPL small-body target designation."""
        self.target = str(target)
        self.timeout = timeout
        self._response_cache = None

    def _validate_response(self, response: Any) -> dict[str, Any]:
        """Validate the decoded API response before row extraction."""
        if not isinstance(response, dict):
            raise RuntimeError(
                "JPL radar API returned an unexpected response for target "
                f"{self.target}: expected a JSON object, got {type(response).__name__}."
            )

        validated_response = dict(response)
        for key in ["fields", "data"]:
            value = validated_response.get(key, [])
            if value is None:
                value = []
            if not isinstance(value, list):
                raise RuntimeError(
                    f"JPL radar API returned invalid '{key}' for target {self.target}: "
                    f"expected a list, got {type(value).__name__}."
                )
            validated_response[key] = value
        return validated_response

    def _fetch_json(self) -> dict[str, Any]:
        """Fetch and cache the raw JPL radar API JSON response."""
        if self._response_cache is not None:
            return self._validate_response(self._response_cache)

        query = urllib.parse.urlencode({"des": self.target})
        with urllib.request.urlopen(
            f"{self.API_URL}?{query}",
            timeout=self.timeout,
        ) as response:
            self._response_cache = self._validate_response(json.load(response))
        return self._response_cache

    def to_radar_tracking_data(
        self,
        target_body: str | None = None,
        epoch_start: datetime.datetime | None = None,
        epoch_end: datetime.datetime | None = None,
        target_point: str | None = "C",
        include_range: bool = True,
        include_doppler: bool = True,
        station_id_mode: str = "jpl",
    ) -> RadarTrackingData:
        """Return JPL radar data as a generic ``RadarTrackingData`` object.

        JPL range rows in microseconds are converted to meters. Doppler rows are
        converted from frequency shift to measured received frequency by adding
        the transmitter frequency. Epochs are retained in UTC and converted to
        TDB for Tudat observation processing.
        """
        response = self._validate_response(self._fetch_json())
        fields = response.get("fields", [])
        rows = []
        target_body = str(target_body or self.target)

        for raw_row in response.get("data", []):
            row = dict(zip(fields, raw_row))
            # Filter by target point, epoch, and observable type before unit conversion.
            if target_point is not None and row.get("bp") != target_point:
                continue

            epoch_datetime = datetime.datetime.strptime(row["epoch"], "%Y-%m-%d %H:%M:%S")
            if epoch_start is not None and epoch_datetime < epoch_start:
                continue
            if epoch_end is not None and epoch_datetime > epoch_end:
                continue

            units = row.get("units")
            if (
                units not in {"us", "Hz"}
                or (units == "us" and not include_range)
                or (units == "Hz" and not include_doppler)
            ):
                continue

            date_time = DateTime.from_iso_string(row["epoch"])
            epoch_utc = date_time.to_epoch()
            epoch_tdb = time_representation.default_time_scale_converter().convert_time(
                input_scale=time_representation.utc_scale,
                output_scale=time_representation.tdb_scale,
                input_value=epoch_utc,
            )

            measurement = float(row["value"])
            sigma = float(row["sigma"])
            frequency_mhz = float(row["freq"])
            transmitter = self._station_id(row["xmit"], station_id_mode)
            receiver = self._station_id(row["rcvr"], station_id_mode)
            # Keep raw JPL values for traceability while storing Tudat-ready values.
            radar_row = {
                "target_body": target_body,
                "epoch_seconds_UTC": epoch_utc,
                "epoch_seconds_TDB": epoch_tdb,
                "transmitter": transmitter,
                "receiver": receiver,
                "target_point": row.get("bp"),
                "transmitter_frequency_hz": np.nan,
                "source": "JPL",
                "raw_value": measurement,
                "raw_sigma": sigma,
                "raw_units": units,
                "jpl_transmitter": str(row["xmit"]).strip(),
                "jpl_receiver": str(row["rcvr"]).strip(),
            }

            if units == "us":
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
            return RadarTrackingData(empty_radar_table())
        return RadarTrackingData(pd.DataFrame(rows))

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
