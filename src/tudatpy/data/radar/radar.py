from __future__ import annotations

import warnings
from collections.abc import Callable, Iterable
from typing import Any

import numpy as np
import pandas as pd

from tudatpy.dynamics import environment
from tudatpy.estimation import observations
from tudatpy.estimation.observable_models_setup import links, model_settings
from tudatpy.estimation.observations_setup import ancillary_settings

from .stations import add_radar_ground_stations

RANGE_OBSERVABLE = "n_way_range"
DOPPLER_OBSERVABLE = "two_way_doppler_frequency"
RADAR_TABLE_META_KEY = "tudatpy_radar_tracking_table"

RADAR_COLUMNS = tuple(
    "target_body epoch_seconds_UTC epoch_seconds_TDB observable_type value sigma "
    "transmitter receiver target_point transmitter_frequency_hz source".split()
)
RADAR_OPTIONAL_COLUMNS = tuple(
    "epoch radar_frequency_mhz raw_units raw_delay_us raw_delay_sigma_us "
    "raw_doppler_shift_hz jpl_transmitter jpl_receiver".split()
)

MPC_RADAR_SPECS = (
    (RANGE_OBSERVABLE, "radar_range", "radar_range_sigma", None),
    (
        DOPPLER_OBSERVABLE,
        "radar_doppler_frequency",
        "radar_doppler_frequency_sigma",
        "radar_frequency_mhz",
    ),
)

REQUIRED_RADAR_COLUMNS = set(RADAR_COLUMNS)
ALLOWED_RADAR_COLUMNS = REQUIRED_RADAR_COLUMNS | set(RADAR_OPTIONAL_COLUMNS)


def empty_radar_table() -> pd.DataFrame:
    """Return an empty table with the canonical radar-tracking schema."""
    return pd.DataFrame(columns=RADAR_COLUMNS)


def _as_string_list(value: str | int | Iterable[str | int] | None) -> list[str] | None:
    """Normalize optional scalar/list filters to strings."""
    if value is None:
        return None
    if isinstance(value, (str, int)):
        return [str(value)]
    return [str(item) for item in value]


def _normalize_station_filter(value: str | int) -> str:
    """Normalize station IDs used for radar table filtering."""
    station_id = str(value).strip()
    return station_id if ":" in station_id else station_id.zfill(3)


def _to_pandas(table: Any) -> pd.DataFrame:
    """Return a pandas copy for pandas/Astropy-like table inputs."""
    if hasattr(table, "to_pandas"):
        return table.to_pandas()
    return table.copy()


def _add_observation_set_to_dataset(
    observation_dataset: observations.ObservationDataset,
    observable_type: model_settings.ObservableType,
    link_ends: dict,
    observation_values: np.ndarray,
    observation_times: np.ndarray,
    reference_link_end: links.LinkEndType,
    ancillary_observation_settings: Any = None,
) -> int:
    """Create a single Tudat observation set and append it to an existing dataset."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        single_observation_set = observations.create_single_observation_set(
            observable_type,
            link_ends,
            observation_values,
            observation_times,
            reference_link_end,
            ancillary_observation_settings,
        )
        single_set_dataset = observations.create_observation_dataset_from_single_observation_set(
            single_observation_set
        )

    return observation_dataset.add_observation_set_from_dataset(single_set_dataset, 0)


def _radar_frequency_band_from_hz(frequency_hz: float) -> Any:
    """Map a radar transmitter frequency to the Tudat frequency-band enum."""
    frequency_mhz = frequency_hz / 1.0e6
    if 2000.0 <= frequency_mhz <= 4000.0:
        return ancillary_settings.FrequencyBands.s_band
    # IEEE radar-band tables define Ku as starting at 12 GHz, so X is half-open.
    if 8000.0 <= frequency_mhz < 12000.0:
        return ancillary_settings.FrequencyBands.x_band
    if 26000.0 <= frequency_mhz <= 40000.0:
        return ancillary_settings.FrequencyBands.ka_band
    if 12000.0 <= frequency_mhz <= 18000.0:
        return ancillary_settings.FrequencyBands.ku_band
    raise ValueError(f"Unsupported radar Doppler transmitter frequency: {frequency_hz} Hz")


def _validate_radar_table(table: pd.DataFrame) -> pd.DataFrame:
    """Validate and normalize the generic radar table before it is stored."""
    if not isinstance(table, pd.DataFrame):
        raise TypeError("Radar tracking data must be provided as a pandas DataFrame.")

    unsupported_columns = set(table.columns).difference(ALLOWED_RADAR_COLUMNS)
    if unsupported_columns:
        raise ValueError(
            "Radar tracking data contains unsupported columns: " f"{sorted(unsupported_columns)}."
        )

    missing_columns = REQUIRED_RADAR_COLUMNS.difference(table.columns)
    if missing_columns:
        raise ValueError(
            "Radar tracking data is missing required columns: " f"{sorted(missing_columns)}."
        )

    validated = table.copy()
    string_columns = "target_body observable_type transmitter receiver target_point source".split()
    validated[string_columns] = validated[string_columns].astype(str)

    invalid_observable_types = set(validated["observable_type"].unique()).difference(
        {RANGE_OBSERVABLE, DOPPLER_OBSERVABLE}
    )
    if invalid_observable_types:
        raise ValueError(
            "Unsupported radar observable types: " f"{sorted(invalid_observable_types)}."
        )

    numeric_columns = (
        "epoch_seconds_UTC epoch_seconds_TDB value sigma transmitter_frequency_hz".split()
    )
    validated[numeric_columns] = validated[numeric_columns].apply(pd.to_numeric)

    if validated["sigma"].isna().any() or (validated["sigma"] <= 0.0).any():
        raise ValueError("Radar sigmas must be finite and strictly positive.")

    doppler_mask = validated["observable_type"] == DOPPLER_OBSERVABLE
    if validated.loc[doppler_mask, "transmitter_frequency_hz"].isna().any():
        raise ValueError("Doppler radar rows must define transmitter_frequency_hz.")

    return validated


class RadarTrackingData:
    """Container for radar range and Doppler tracking data.

    The stored table uses a radar-native schema and deliberately excludes angular
    astrometry columns such as RA and DEC. Values are in Tudat units: n-way range
    in meters and measured Doppler frequency in hertz.
    """

    def __init__(self, table: pd.DataFrame) -> None:
        """Create a radar data container from a table with ``RADAR_COLUMNS``."""
        self._table = _validate_radar_table(table)

    @classmethod
    def empty(cls) -> "RadarTrackingData":
        """Return an empty radar data container."""
        return cls(empty_radar_table())

    @classmethod
    def concatenate(cls, radar_data_sets: Iterable["RadarTrackingData"]) -> "RadarTrackingData":
        """Concatenate multiple radar data containers."""
        frames = [radar_data.table for radar_data in radar_data_sets if len(radar_data) > 0]
        if not frames:
            return cls.empty()
        return cls(pd.concat(frames, ignore_index=True, sort=False))

    @property
    def table(self) -> pd.DataFrame:
        """Return a copy of the stored radar observation table."""
        return self._table.copy()

    def __len__(self) -> int:
        return len(self._table)

    def filter(
        self,
        *,
        epoch_start: float | None = None,
        epoch_end: float | None = None,
        target_body: str | Iterable[str] | None = None,
        target_point: str | Iterable[str] | None = None,
        observable_type: str | Iterable[str] | None = None,
        station_ids: Iterable[str | int] | None = None,
        exclude_station_ids: Iterable[str | int] | None = None,
    ) -> "RadarTrackingData":
        """Return a filtered radar data container.

        Epoch filters are applied to TDB seconds since J2000. The remaining
        filters are exact string matches after normalizing the requested values.
        """
        filtered = self._table.copy()
        for column, value in {
            "target_body": target_body,
            "target_point": target_point,
            "observable_type": observable_type,
        }.items():
            allowed_values = _as_string_list(value)
            if allowed_values is not None:
                filtered = filtered.loc[filtered[column].isin(allowed_values)]

        for station_filter, keep_matches in [
            (station_ids, True),
            (exclude_station_ids, False),
        ]:
            if station_filter is None:
                continue
            normalized_ids = {_normalize_station_filter(station) for station in station_filter}
            station_mask = filtered["transmitter"].isin(normalized_ids) | filtered["receiver"].isin(
                normalized_ids
            )
            filtered = filtered.loc[station_mask if keep_matches else ~station_mask]

        if epoch_start is not None:
            filtered = filtered.loc[filtered["epoch_seconds_TDB"] >= epoch_start]
        if epoch_end is not None:
            filtered = filtered.loc[filtered["epoch_seconds_TDB"] <= epoch_end]
        return RadarTrackingData(filtered.reset_index(drop=True))

    def with_target_body(self, target_body: str) -> "RadarTrackingData":
        """Return a copy with all rows assigned to one target body name."""
        table = self.table
        if table.empty:
            return RadarTrackingData.empty()
        table["target_body"] = str(target_body)
        return RadarTrackingData(table)

    def station_id_list(self) -> list[str]:
        """Return transmitter/receiver station IDs in first-appearance order."""
        if self._table.empty:
            return []
        station_ids = self._table[["transmitter", "receiver"]].stack().dropna().astype(str)
        return list(dict.fromkeys(station_ids))

    def target_body_list(self) -> list[str]:
        """Return target-body names in first-appearance order."""
        return list(dict.fromkeys(self._table["target_body"].dropna().astype(str)))

    def station_observation_counts(self) -> pd.DataFrame:
        """Return per-station observation counts, counting each radar row once per station."""
        if self._table.empty:
            return pd.DataFrame(columns=["observatory", "count"])
        return (
            self._table.loc[:, ["transmitter", "receiver"]]
            .reset_index()
            .melt(id_vars="index", value_name="observatory")
            .drop_duplicates(["index", "observatory"])
            .groupby("observatory")
            .size()
            .rename("count")
            .reset_index(drop=False)
        )

    def required_station_ids(self) -> set[str]:
        """Return all transmitter and receiver station IDs used by the data."""
        return set(self._table[["transmitter", "receiver"]].stack().dropna().astype(str))

    def required_target_bodies(self) -> set[str]:
        """Return all target-body names used by the data."""
        return set(self._table["target_body"].dropna().astype(str))

    def to_tudat(
        self,
        bodies: Any,
        station_body: str = "Earth",
        add_ground_stations: bool = True,
        station_positions: dict[str, np.ndarray] | None = None,
    ) -> observations.ObservationDataset:
        """Convert the radar data to a Tudat ``ObservationDataset``."""
        observation_dataset = observations.ObservationDataset()
        self.add_to_observation_dataset(
            observation_dataset,
            bodies,
            station_body=station_body,
            add_ground_stations=add_ground_stations,
            station_positions=station_positions,
        )
        return observation_dataset

    def add_to_observation_dataset(
        self,
        observation_dataset: observations.ObservationDataset,
        bodies: Any,
        station_body: str = "Earth",
        add_ground_stations: bool = True,
        station_positions: dict[str, np.ndarray] | None = None,
    ) -> None:
        """Append the radar data to an existing Tudat observation dataset.

        Only center-of-mass radar returns (``target_point == "C"``) are added.
        The method optionally creates known radar ground stations, creates empty
        target bodies if needed, configures Doppler transmitter/transponder
        settings, and then adds one observation set per unique link geometry.
        """
        radar_table = self._table.loc[self._table["target_point"] == "C"].copy()
        if radar_table.empty:
            return

        # The station positions are needed by Tudat before link ends can be created.
        if add_ground_stations:
            add_radar_ground_stations(
                bodies,
                set(radar_table["transmitter"].astype(str)).union(
                    set(radar_table["receiver"].astype(str))
                ),
                station_body=station_body,
                station_positions=station_positions,
            )

        # Radar observations can be loaded before an ephemeris is attached to the target.
        for target_body in sorted(set(radar_table["target_body"].astype(str))):
            if not bodies.does_body_exist(target_body):
                bodies.create_empty_body(target_body)

        # Doppler observations require transmitter frequencies and target turnaround ratios.
        self._configure_doppler_environment(radar_table, bodies, station_body)
        for radar_type, tudat_type, ancillary_factory in [
            (RANGE_OBSERVABLE, model_settings.n_way_range_type, None),
            (
                DOPPLER_OBSERVABLE,
                model_settings.doppler_measured_frequency_type,
                self._doppler_ancillary_settings,
            ),
        ]:
            self._add_radar_sets(
                radar_table,
                observation_dataset,
                station_body,
                radar_type,
                tudat_type,
                ancillary_factory,
            )

    @staticmethod
    def _configure_doppler_environment(
        radar_table: pd.DataFrame, bodies: Any, station_body: str
    ) -> None:
        """Attach Doppler frequency models required by Tudat's Doppler observable."""
        doppler_table = radar_table.loc[
            (radar_table["observable_type"] == DOPPLER_OBSERVABLE)
            & radar_table["transmitter_frequency_hz"].notna()
        ].copy()
        if doppler_table.empty:
            return

        for transmitter, transmitter_rows in doppler_table.groupby("transmitter"):
            transmitter_frequencies_hz = (
                transmitter_rows["transmitter_frequency_hz"].dropna().unique()
            )
            if len(transmitter_frequencies_hz) != 1:
                raise ValueError(
                    "Radar Doppler conversion currently supports one transmitter "
                    f"frequency per station; station {transmitter} has "
                    f"{len(transmitter_frequencies_hz)} frequencies."
                )
            bodies.get_body(station_body).get_ground_station(
                transmitter
            ).set_transmitting_frequency_calculator(
                environment.ConstantTransmittingFrequencyCalculator(
                    float(transmitter_frequencies_hz[0])
                )
            )

        # The small body acts as a unit-ratio transponder for center-of-mass returns.
        for target_body, target_rows in doppler_table.groupby("target_body"):
            turnaround_ratios = {}
            for frequency_hz in target_rows["transmitter_frequency_hz"].dropna().unique():
                frequency_band = _radar_frequency_band_from_hz(float(frequency_hz))
                turnaround_ratios[(frequency_band, frequency_band)] = 1.0
            vehicle_systems = environment.VehicleSystems()
            vehicle_systems.set_transponder_turnaround_ratio(turnaround_ratios)
            bodies.get_body(str(target_body)).system_models = vehicle_systems

    @staticmethod
    def _doppler_ancillary_settings(group: dict[str, Any]):
        """Build the Doppler ancillary settings for a grouped radar link."""
        transmitter_frequency_band = _radar_frequency_band_from_hz(
            float(group["transmitter_frequency_hz"])
        )
        return ancillary_settings.doppler_measured_frequency_ancillary_settings(
            [transmitter_frequency_band, transmitter_frequency_band]
        )

    @staticmethod
    def _add_radar_sets(
        radar_table: pd.DataFrame,
        observation_dataset: observations.ObservationDataset,
        station_body: str,
        radar_observable_type: str,
        tudat_observable_type: model_settings.ObservableType,
        ancillary_factory: Callable[[dict[str, Any]], Any] | None = None,
    ) -> None:
        """Group radar rows into Tudat observation sets for one observable type."""
        observable_table = radar_table.loc[
            radar_table["observable_type"] == radar_observable_type
        ].copy()
        if observable_table.empty:
            return

        group_columns = ["target_body", "transmitter", "receiver"] + (
            ["transmitter_frequency_hz"] if radar_observable_type == DOPPLER_OBSERVABLE else []
        )

        for group_values, observations_for_link in observable_table.groupby(
            group_columns, sort=False
        ):
            # Each group has a single link geometry and, for Doppler, a single frequency.
            group = dict(zip(group_columns, group_values))
            observations_for_link = observations_for_link.sort_values(
                "epoch_seconds_TDB", kind="stable"
            )
            link_ends = {
                links.transmitter: links.body_reference_point_link_end_id(
                    station_body, group["transmitter"]
                ),
                links.retransmitter: links.body_origin_link_end_id(group["target_body"]),
                links.receiver: links.body_reference_point_link_end_id(
                    station_body, group["receiver"]
                ),
            }
            ancillary_observation_settings = (
                ancillary_factory(group) if ancillary_factory is not None else None
            )
            # Tudat expects one scalar column and TDB observation times for each set.
            set_id = _add_observation_set_to_dataset(
                observation_dataset,
                tudat_observable_type,
                link_ends,
                observations_for_link.loc[:, ["value"]].to_numpy(),
                observations_for_link.loc[:, "epoch_seconds_TDB"].to_numpy(),
                links.receiver,
                ancillary_observation_settings,
            )
            observation_dataset.set_weight_vector_for_set(
                set_id,
                1.0 / observations_for_link.loc[:, "sigma"].to_numpy() ** 2,
            )


def _normalize_mpc_station_series(station_series: pd.Series) -> pd.Series:
    """Normalize MPC radar station IDs to three-character station codes."""
    return station_series.astype(str).str.strip().str.zfill(3)


def _mpc_radar_frame(
    table: pd.DataFrame,
    observable_type: str,
    value_column: str,
    sigma_column: str,
    frequency_column: str | None,
) -> pd.DataFrame:
    """Convert one MPC radar observable family to the generic radar schema."""
    # Only rows with a complete radar link, measurement, and sigma can be used.
    required_observation_columns = [
        "radar_target_point",
        "radar_transmitter",
        "radar_receiver",
        value_column,
        sigma_column,
    ]
    if frequency_column is not None:
        required_observation_columns.append(frequency_column)
    required_columns = {
        "number",
        "epoch_seconds_UTC",
        "epoch_seconds_TDB",
        *required_observation_columns,
    }
    if not required_columns.issubset(table.columns):
        return empty_radar_table()
    rows = table.loc[table.loc[:, required_observation_columns].notna().all(axis=1)].copy()
    if rows.empty:
        return empty_radar_table()

    return pd.DataFrame(
        {
            "target_body": rows["number"].astype(str),
            "epoch_seconds_UTC": rows["epoch_seconds_UTC"],
            "epoch_seconds_TDB": rows["epoch_seconds_TDB"],
            "observable_type": observable_type,
            "value": rows[value_column],
            "sigma": rows[sigma_column],
            "transmitter": _normalize_mpc_station_series(rows["radar_transmitter"]),
            "receiver": _normalize_mpc_station_series(rows["radar_receiver"]),
            "target_point": rows["radar_target_point"],
            "transmitter_frequency_hz": (
                pd.to_numeric(rows[frequency_column]) * 1.0e6
                if frequency_column is not None
                else np.nan
            ),
            "source": "MPC",
        }
    )


def radar_tracking_data_from_mpc_table(table: pd.DataFrame) -> RadarTrackingData:
    """Convert a BatchMPC-style table with radar columns to ``RadarTrackingData``."""
    frames = [_mpc_radar_frame(table, *spec) for spec in MPC_RADAR_SPECS]
    frames = [frame for frame in frames if not frame.empty]
    table = pd.concat(frames, ignore_index=True, sort=False) if frames else empty_radar_table()
    return RadarTrackingData(table)


def radar_tracking_data_from_table(table: Any) -> RadarTrackingData:
    """Extract radar tracking data from radar-native, metadata, or legacy MPC tables."""
    if hasattr(table, "meta") and RADAR_TABLE_META_KEY in table.meta:
        radar_table = table.meta[RADAR_TABLE_META_KEY]
        if isinstance(radar_table, pd.DataFrame):
            return RadarTrackingData(radar_table)
        return RadarTrackingData(_to_pandas(radar_table))

    dataframe = _to_pandas(table)
    if any(column in dataframe.columns for column in ["radar_range", "radar_doppler_frequency"]):
        return radar_tracking_data_from_mpc_table(dataframe)
    if set(RADAR_COLUMNS).issubset(dataframe.columns) and not {"RA", "DEC"}.intersection(
        dataframe.columns
    ):
        return RadarTrackingData(dataframe)
    return RadarTrackingData.empty()


def remove_radar_tracking_data_from_table(table: Any) -> pd.DataFrame:
    """Return table rows that are not legacy flat radar observations."""
    dataframe = _to_pandas(table)
    if set(RADAR_COLUMNS).issubset(dataframe.columns) and not {"RA", "DEC"}.intersection(
        dataframe.columns
    ):
        return pd.DataFrame()

    radar_columns = [
        column
        for column in ["radar_target_point", "radar_range", "radar_doppler_frequency"]
        if column in dataframe
    ]
    if not radar_columns:
        return dataframe
    return dataframe.loc[~dataframe.loc[:, radar_columns].notna().any(axis=1)].copy()
