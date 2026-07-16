"""Canonical radar tracking-data utilities."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import numpy as np
import pandas as pd

from tudatpy.data_input.tracking_data import (
    PiecewiseConstantFrequencySupplementaryData,
    TrackingData,
    TrackingSupplementaryData,
)

RANGE_OBSERVABLE = "n_way_range"
DOPPLER_OBSERVABLE = "doppler_measured_frequency"
RADAR_TABLE_META_KEY = "tudatpy_radar_tracking_table"

RADAR_COLUMNS = tuple(
    "target_body epoch_seconds_UTC observable_type value sigma "
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

_TRACKING_OBSERVABLE_TYPES = {
    RANGE_OBSERVABLE: "NWayRange",
    DOPPLER_OBSERVABLE: "DopplerMeasuredFrequency",
}


def empty_radar_table() -> pd.DataFrame:
    """Return an empty canonical radar-data table.

    Returns
    -------
    pandas.DataFrame
        Empty table with the columns used by all TudatPy radar readers. Epochs
        are UTC seconds since J2000, ranges are metres, Doppler observations
        are received frequencies in hertz, and sigmas use the same units as the
        corresponding observation values.
    """
    return pd.DataFrame(columns=RADAR_COLUMNS)


def radar_frequency_band_string_from_hz(frequency_hz: float) -> str:
    """Map a radar transmitter frequency to a Tudat frequency-band string.

    Parameters
    ----------
    frequency_hz : float
        Transmitter frequency in hertz.

    Returns
    -------
    str
        Tudat frequency-band label corresponding to ``frequency_hz``.

    Raises
    ------
    ValueError
        If the frequency is outside the currently supported radar bands.
    """
    frequency_mhz = float(frequency_hz) / 1.0e6
    if 2000.0 <= frequency_mhz <= 4000.0:
        return "S-band"
    if 8000.0 <= frequency_mhz < 12000.0:
        return "X-band"
    if 12000.0 <= frequency_mhz <= 18000.0:
        return "Ku-band"
    if 26000.0 <= frequency_mhz <= 40000.0:
        return "Ka-band"
    raise ValueError(f"Unsupported radar Doppler transmitter frequency: {frequency_hz} Hz")


def _as_string_list(value: str | int | Iterable[str | int] | None) -> list[str] | None:
    """Normalize scalar or iterable filter values to a list of strings."""
    if value is None:
        return None
    if isinstance(value, (str, int)):
        return [str(value)]
    return [str(item) for item in value]


def _normalize_station_filter(value: str | int) -> str:
    """Normalize station filters to the canonical table representation."""
    station_id = str(value).strip()
    return station_id if ":" in station_id else station_id.zfill(3)


def _to_pandas(table: Any) -> pd.DataFrame:
    """Convert supported table-like inputs to a pandas DataFrame."""
    if hasattr(table, "to_pandas"):
        return table.to_pandas()
    return table.copy()


def validate_radar_data(table: pd.DataFrame) -> pd.DataFrame:
    """Validate and normalize canonical radar data.

    Parameters
    ----------
    table : pandas.DataFrame
        Radar data table using the columns listed in ``RADAR_COLUMNS``. Extra
        columns are allowed only if they are listed in
        ``RADAR_OPTIONAL_COLUMNS``.

    Returns
    -------
    pandas.DataFrame
        Validated copy of ``table`` with normalized string columns, numeric
        columns coerced to numeric dtypes, and a reset row index.

    Raises
    ------
    TypeError
        If ``table`` is not a pandas DataFrame.
    ValueError
        If required columns are missing, unsupported columns are present, an
        observable type is unsupported, a sigma is invalid, or a Doppler row
        lacks its transmitter frequency.
    """
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

    numeric_columns = "epoch_seconds_UTC value sigma transmitter_frequency_hz".split()
    validated[numeric_columns] = validated[numeric_columns].apply(pd.to_numeric)

    if validated["sigma"].isna().any() or (validated["sigma"] <= 0.0).any():
        raise ValueError("Radar sigmas must be finite and strictly positive.")

    doppler_mask = validated["observable_type"] == DOPPLER_OBSERVABLE
    if validated.loc[doppler_mask, "transmitter_frequency_hz"].isna().any():
        raise ValueError("Doppler radar rows must define transmitter_frequency_hz.")

    return validated.reset_index(drop=True)


def filter_radar_data(
    radar_table: pd.DataFrame,
    *,
    epoch_start: float | None = None,
    epoch_end: float | None = None,
    target_body: str | Iterable[str] | None = None,
    target_point: str | Iterable[str] | None = None,
    observable_type: str | Iterable[str] | None = None,
    station_ids: Iterable[str | int] | None = None,
    exclude_station_ids: Iterable[str | int] | None = None,
) -> pd.DataFrame:
    """Return filtered canonical radar data.

    Epoch filters are applied to UTC seconds since J2000. Target, observable,
    and target-point filters are exact string matches. Station filters match
    either the transmitter or receiver station after normalizing bare numeric
    station identifiers to three-character MPC-style codes.

    Parameters
    ----------
    radar_table : pandas.DataFrame
        Canonical radar data table.
    epoch_start : float | None, default None
        Optional lower bound on ``epoch_seconds_UTC``.
    epoch_end : float | None, default None
        Optional upper bound on ``epoch_seconds_UTC``.
    target_body : str | iterable[str] | None, default None
        Target body or bodies to retain.
    target_point : str | iterable[str] | None, default None
        Radar target point or points to retain.
    observable_type : str | iterable[str] | None, default None
        Radar observable type or types to retain.
    station_ids : iterable[str | int] | None, default None
        Station identifiers to include.
    exclude_station_ids : iterable[str | int] | None, default None
        Station identifiers to exclude.

    Returns
    -------
    pandas.DataFrame
        Filtered canonical radar data table.
    """
    filtered = validate_radar_data(radar_table)
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
        filtered = filtered.loc[filtered["epoch_seconds_UTC"] >= epoch_start]
    if epoch_end is not None:
        filtered = filtered.loc[filtered["epoch_seconds_UTC"] <= epoch_end]
    return validate_radar_data(filtered.reset_index(drop=True))


def set_radar_target_body(radar_table: pd.DataFrame, target_body: str) -> pd.DataFrame:
    """Return radar data with all rows assigned to one target body name.

    Parameters
    ----------
    radar_table : pandas.DataFrame
        Canonical radar data table.
    target_body : str
        Target body name to assign to all radar rows.

    Returns
    -------
    pandas.DataFrame
        Copy of ``radar_table`` with the ``target_body`` column replaced.
    """
    validated = validate_radar_data(radar_table)
    if validated.empty:
        return empty_radar_table()
    validated["target_body"] = str(target_body)
    return validate_radar_data(validated)


def radar_station_id_list(radar_table: pd.DataFrame) -> list[str]:
    """Return transmitter/receiver station IDs in first-appearance order.

    Parameters
    ----------
    radar_table : pandas.DataFrame
        Canonical radar data table.

    Returns
    -------
    list[str]
        Unique station identifiers from the transmitter and receiver columns.
    """
    validated = validate_radar_data(radar_table)
    if validated.empty:
        return []
    station_ids = validated[["transmitter", "receiver"]].stack().dropna().astype(str)
    return list(dict.fromkeys(station_ids))


def radar_target_body_list(radar_table: pd.DataFrame) -> list[str]:
    """Return target-body names in first-appearance order.

    Parameters
    ----------
    radar_table : pandas.DataFrame
        Canonical radar data table.

    Returns
    -------
    list[str]
        Unique target-body names in their first-appearance order.
    """
    validated = validate_radar_data(radar_table)
    return list(dict.fromkeys(validated["target_body"].dropna().astype(str)))


def radar_station_observation_counts(radar_table: pd.DataFrame) -> pd.DataFrame:
    """Return per-station observation counts.

    A row is counted once for each distinct transmitter/receiver station that
    appears in that row. Monostatic observations therefore contribute one count
    to the station; bistatic observations contribute one count to each station.

    Parameters
    ----------
    radar_table : pandas.DataFrame
        Canonical radar data table.

    Returns
    -------
    pandas.DataFrame
        Table with columns ``observatory`` and ``count``.
    """
    validated = validate_radar_data(radar_table)
    if validated.empty:
        return pd.DataFrame(columns=["observatory", "count"])
    return (
        validated.loc[:, ["transmitter", "receiver"]]
        .reset_index()
        .melt(id_vars="index", value_name="observatory")
        .drop_duplicates(["index", "observatory"])
        .groupby("observatory")
        .size()
        .rename("count")
        .reset_index(drop=False)
    )


def _radar_frequency_supplementary_data(
    radar_table: pd.DataFrame,
    station_body: str,
) -> list[TrackingSupplementaryData]:
    """Create supplementary transmitter-frequency data for radar Doppler rows."""
    doppler_table = radar_table.loc[
        (radar_table["observable_type"] == DOPPLER_OBSERVABLE)
        & radar_table["transmitter_frequency_hz"].notna()
    ].copy()
    if doppler_table.empty:
        return []

    supplementary_data = []
    for transmitter, group in doppler_table.groupby("transmitter", sort=False):
        history = {
            float(epoch): float(frequency)
            for epoch, frequency in zip(
                group["epoch_seconds_UTC"],
                group["transmitter_frequency_hz"],
            )
        }
        frequency_data = PiecewiseConstantFrequencySupplementaryData(history)
        station_data = TrackingSupplementaryData(station_body, str(transmitter))
        station_data.set_frequency_supplementary_data([frequency_data])
        supplementary_data.append(station_data)
    return supplementary_data


def radar_data_to_tracking_data(
    radar_table: pd.DataFrame,
    station_body: str = "Earth",
    target_point: str | None = "C",
) -> tuple[list[TrackingData], list[TrackingSupplementaryData]]:
    """Convert canonical radar data to Tudat tracking and supplementary data.

    Parameters
    ----------
    radar_table : pandas.DataFrame
        Canonical radar data table. Range values must be expressed in metres.
        Doppler values must be expressed as measured received frequencies in
        hertz. Epochs are interpreted as UTC seconds since J2000.
    station_body : str, default "Earth"
        Body on which transmitter and receiver station reference points are
        defined.
    target_point : str | None, default "C"
        Target point to convert. If None, all target points are converted.

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        Tudat tracking data objects and supplementary data objects. Doppler
        transmitter-frequency histories are returned as supplementary data.
    """
    radar_table = validate_radar_data(radar_table)
    if target_point is not None:
        radar_table = radar_table.loc[radar_table["target_point"] == target_point].copy()
    if radar_table.empty:
        return [], []

    tracking_data = []
    group_columns = ["target_body", "transmitter", "receiver", "observable_type"]
    doppler_column = "transmitter_frequency_hz"

    for group_values, group in radar_table.groupby(group_columns, sort=False):
        target_body, transmitter, receiver, observable_type = group_values
        group = group.sort_values("epoch_seconds_UTC", kind="stable")
        link_ends = [
            ((station_body, str(transmitter)), "transmitter"),
            ((str(target_body), ""), "reflector_1"),
            ((station_body, str(receiver)), "receiver"),
        ]

        data_object = TrackingData(
            observable_type=_TRACKING_OBSERVABLE_TYPES[observable_type],
            link_ends=link_ends,
            observations=[np.array([value]) for value in group["value"]],
            epochs=group["epoch_seconds_UTC"].tolist(),
            reference_link_end="receiver",
            time_scale="UTC",
        )

        if observable_type == DOPPLER_OBSERVABLE:
            # Tudat currently expects a single frequency-band ancillary setting
            # per Doppler TrackingData object. Split by transmitter/receiver and
            # observable first, then reject mixed-band groups explicitly.
            bands = [
                radar_frequency_band_string_from_hz(frequency_hz)
                for frequency_hz in group[doppler_column].to_numpy(dtype=float)
            ]
            unique_bands = list(dict.fromkeys(bands))
            if len(unique_bands) != 1:
                raise ValueError(
                    "Doppler radar conversion currently groups one frequency band "
                    "per TrackingData object."
                )
            data_object.add_string_vector_ancillary_setting(
                "frequency bands",
                [unique_bands[0], unique_bands[0]],
            )

        tracking_data.append(data_object)

    return tracking_data, _radar_frequency_supplementary_data(
        radar_table,
        station_body=station_body,
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


def radar_data_from_mpc_table(table: pd.DataFrame) -> pd.DataFrame:
    """Convert a legacy flat MPC table to canonical radar data.

    Parameters
    ----------
    table : pandas.DataFrame
        Table containing MPC-style radar columns such as ``radar_range`` or
        ``radar_doppler_frequency``.

    Returns
    -------
    pandas.DataFrame
        Canonical radar data table.
    """
    frames = [_mpc_radar_frame(table, *spec) for spec in MPC_RADAR_SPECS]
    frames = [frame for frame in frames if not frame.empty]
    radar_table = (
        pd.concat(frames, ignore_index=True, sort=False) if frames else empty_radar_table()
    )
    return validate_radar_data(radar_table)


def radar_data_from_table(table: Any) -> pd.DataFrame:
    """Extract canonical radar data from a supported table representation.

    Parameters
    ----------
    table : Any
        Astropy table, pandas DataFrame, or compatible table-like object. Radar
        data may be provided directly with the canonical columns, stored in the
        ``RADAR_TABLE_META_KEY`` metadata entry, or represented by legacy flat
        MPC radar columns.

    Returns
    -------
    pandas.DataFrame
        Canonical radar data table. If no radar data are present, an empty
        canonical table is returned.
    """
    if hasattr(table, "meta") and RADAR_TABLE_META_KEY in table.meta:
        radar_table = table.meta[RADAR_TABLE_META_KEY]
        if isinstance(radar_table, pd.DataFrame):
            return validate_radar_data(radar_table)
        return validate_radar_data(_to_pandas(radar_table))

    dataframe = _to_pandas(table)
    if any(column in dataframe.columns for column in ["radar_range", "radar_doppler_frequency"]):
        return radar_data_from_mpc_table(dataframe)
    if set(RADAR_COLUMNS).issubset(dataframe.columns) and not {"RA", "DEC"}.intersection(
        dataframe.columns
    ):
        return validate_radar_data(dataframe)
    return empty_radar_table()
