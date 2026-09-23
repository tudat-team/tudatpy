"""Canonical radar tracking-data table and its conversion to Tudat tracking data.

The table has one row per observation: epochs in UTC seconds since J2000, ranges as
round-trip light time times the speed of light [m], Doppler as received frequency [Hz].
"""

import warnings

import numpy as np
import pandas as pd

from tudatpy import constants
from tudatpy.data_input.tracking_data import (
    PiecewiseConstantFrequencySupplementaryData,
    TrackingData,
    TrackingSupplementaryData,
)
from tudatpy.data_input.tracking_data.optical_utilities import datetime_to_utc_seconds

RANGE_OBSERVABLE = "NWayRange"
DOPPLER_OBSERVABLE = "DopplerMeasuredFrequency"
RADAR_TABLE_META_KEY = "tudatpy_radar_tracking_table"
RADAR_COLUMNS = [
    "target_body",
    "epoch_seconds_UTC",
    "observable_type",
    "value",
    "sigma",
    "transmitter",
    "receiver",
    "target_point",
    "transmitter_frequency_hz",
    "source",
]

_TUDAT_BANDS = {
    "S-band": (2.0e9, 4.0e9),
    "X-band": (7.0e9, 12.0e9),
    "Ku-band": (12.0e9, 18.0e9),
    "Ka-band": (26.5e9, 40.0e9),
}


def empty_radar_table() -> pd.DataFrame:
    """Return an empty canonical radar table.

    Returns
    -------
    pandas.DataFrame
        Empty table with the columns listed in :data:`RADAR_COLUMNS`.
    """
    return pd.DataFrame(columns=RADAR_COLUMNS)


def validate_radar_data(table: pd.DataFrame) -> pd.DataFrame:
    """Check a canonical radar table and return it with numeric columns as floats.

    Parameters
    ----------
    table : pandas.DataFrame
        Radar table to validate. It must contain every column listed in
        :data:`RADAR_COLUMNS`.

    Returns
    -------
    pandas.DataFrame
        Validated table with a consecutive integer index and numeric columns
        converted to ``float``.

    Raises
    ------
    ValueError
        If columns are missing, an observable type is unknown, a sigma is not
        positive, or a Doppler row has no transmitter frequency.
    """
    missing = set(RADAR_COLUMNS) - set(table.columns)
    if missing:
        raise ValueError(f"Radar table is missing columns: {sorted(missing)}")
    table = table.astype(
        {
            "epoch_seconds_UTC": float,
            "value": float,
            "sigma": float,
            "transmitter_frequency_hz": float,
        }
    )
    unknown = set(table["observable_type"]) - {RANGE_OBSERVABLE, DOPPLER_OBSERVABLE}
    if unknown:
        raise ValueError(f"Unsupported radar observable types: {sorted(unknown)}")
    if not (table["sigma"] > 0.0).all():
        raise ValueError("Radar sigmas must be strictly positive.")
    is_doppler = table["observable_type"] == DOPPLER_OBSERVABLE
    if table.loc[is_doppler, "transmitter_frequency_hz"].isna().any():
        raise ValueError("Doppler radar rows must define transmitter_frequency_hz.")
    return table.reset_index(drop=True)


def radar_data_from_raw(raw: pd.DataFrame, source: str) -> pd.DataFrame:
    """Convert delay and Doppler-shift measurements to the canonical radar table.

    Parameters
    ----------
    raw : pandas.DataFrame
        Columns ``target_body``, ``epoch_seconds_UTC``, ``transmitter``, ``receiver``,
        ``target_point``, ``transmitter_frequency_hz`` [Hz], ``delay_us``,
        ``delay_sigma_us`` [us], ``doppler_hz``, and ``doppler_sigma_hz`` [Hz].
        A row may contain a delay, a Doppler shift, or both; absent values are NaN.
    source : str
        Data source label, for example ``"MPC"`` or ``"JPL"``.

    Returns
    -------
    pandas.DataFrame
        Canonical radar table with one row per delay or Doppler measurement.
    """
    common = raw[
        [
            "target_body",
            "epoch_seconds_UTC",
            "transmitter",
            "receiver",
            "target_point",
            "transmitter_frequency_hz",
        ]
    ].assign(source=source)
    is_range = raw["delay_us"].notna()
    is_doppler = raw["doppler_hz"].notna()
    range_rows = common[is_range].assign(
        observable_type=RANGE_OBSERVABLE,
        value=constants.SPEED_OF_LIGHT * 1.0e-6 * raw.loc[is_range, "delay_us"],
        sigma=constants.SPEED_OF_LIGHT * 1.0e-6 * raw.loc[is_range, "delay_sigma_us"],
    )
    doppler_rows = common[is_doppler].assign(
        observable_type=DOPPLER_OBSERVABLE,
        value=(raw.loc[is_doppler, "transmitter_frequency_hz"] + raw.loc[is_doppler, "doppler_hz"]),
        sigma=raw.loc[is_doppler, "doppler_sigma_hz"],
    )
    if range_rows.empty and doppler_rows.empty:
        return empty_radar_table()
    return validate_radar_data(
        pd.concat([range_rows, doppler_rows], ignore_index=True)[RADAR_COLUMNS]
    )


def radar_data_from_table(table) -> pd.DataFrame:
    """Return the radar data stored in a parsed MPC 80-column table.

    Parameters
    ----------
    table : astropy.table.Table
        Table returned by an MPC 80-column parser.

    Returns
    -------
    pandas.DataFrame
        Canonical radar table stored under :data:`RADAR_TABLE_META_KEY`, or an
        empty canonical radar table if the metadata contain no radar data.
    """
    return getattr(table, "meta", {}).get(RADAR_TABLE_META_KEY, empty_radar_table())


def filter_radar_data(
    table: pd.DataFrame,
    *,
    epoch_start=None,
    epoch_end=None,
    target_body=None,
    target_point=None,
    observable_type=None,
    station_ids=None,
    exclude_station_ids=None,
) -> pd.DataFrame:
    """Filter a canonical radar table.

    Parameters
    ----------
    table : pandas.DataFrame
        Canonical radar table to filter.
    epoch_start : float | DateTime | Time | datetime.datetime | None, default None
        Inclusive lower epoch bound. Numeric values are UTC seconds since J2000.
    epoch_end : float | DateTime | Time | datetime.datetime | None, default None
        Inclusive upper epoch bound. Numeric values are UTC seconds since J2000.
    target_body : str | iterable[str] | None, default None
        Target body name or names to retain.
    target_point : str | iterable[str] | None, default None
        Radar bounce-point code or codes to retain.
    observable_type : str | iterable[str] | None, default None
        Tudat observable type or types to retain.
    station_ids : iterable[str | int] | None, default None
        Retain rows that use one of these transmitters or receivers.
    exclude_station_ids : iterable[str | int] | None, default None
        Remove rows that use one of these transmitters or receivers.

    Returns
    -------
    pandas.DataFrame
        Filtered canonical radar table with a consecutive integer index.

    Notes
    -----
    Bare numeric station codes are padded to three characters following the
    MPC convention.
    """
    keep = pd.Series(True, index=table.index)
    for column, allowed in [
        ("target_body", target_body),
        ("target_point", target_point),
        ("observable_type", observable_type),
    ]:
        if allowed is not None:
            keep &= table[column].isin([allowed] if isinstance(allowed, str) else list(allowed))

    def uses_station(ids):
        ids = {str(station).strip().zfill(3) for station in ids}
        return table["transmitter"].isin(ids) | table["receiver"].isin(ids)

    if station_ids is not None:
        keep &= uses_station(station_ids)
    if exclude_station_ids is not None:
        keep &= ~uses_station(exclude_station_ids)
    if epoch_start is not None:
        keep &= table["epoch_seconds_UTC"] >= datetime_to_utc_seconds(epoch_start)
    if epoch_end is not None:
        keep &= table["epoch_seconds_UTC"] <= datetime_to_utc_seconds(epoch_end)
    return table[keep].reset_index(drop=True)


def radar_frequency_band_string_from_hz(frequency_hz: float) -> str:
    """Return the Tudat frequency-band label nearest to a transmitter frequency.

    Tudat defines S, X, Ku and Ka bands. Frequencies outside those ranges are
    assigned to the nearest band on a logarithmic scale.

    Parameters
    ----------
    frequency_hz : float
        Positive transmitter frequency [Hz].

    Returns
    -------
    str
        One of ``"S-band"``, ``"X-band"``, ``"Ku-band"`` or ``"Ka-band"``.
    """
    distance = {
        band: max(np.log(lower / frequency_hz), np.log(frequency_hz / upper), 0.0)
        for band, (lower, upper) in _TUDAT_BANDS.items()
    }
    return min(distance, key=distance.get)


def radar_data_to_tracking_data(
    table: pd.DataFrame,
    station_body: str = "Earth",
    target_point: str | None = "C",
):
    """Convert a canonical radar table to Tudat tracking and supplementary data.

    Parameters
    ----------
    table : pandas.DataFrame
        Canonical radar table.
    station_body : str, default "Earth"
        Body on which the transmitting and receiving stations are defined.
    target_point : str | None, default "C"
        Bounce point to convert. None converts all rows.

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        One TrackingData object per target, transmitter, receiver, observable and
        frequency band, plus the transmitter-frequency histories used by Doppler.

    Notes
    -----
    Radar delays are measured with station clocks. When computing residuals for
    these data, create the range observation models with
    ``model_settings.n_way_range(..., time_scale_for_observable=time_representation.utc_scale)``.
    Import ``model_settings`` from the public
    ``tudatpy.estimation.observable_models_setup`` module. With the default
    (TDB), ranges are biased by about 1.5e-8 times the round-trip light time
    (roughly 0.5 km for a target at 0.1 au).
    """
    if target_point is not None:
        skipped = table["target_point"] != target_point
        if skipped.any():
            warnings.warn(
                f"Skipping {int(skipped.sum())} radar observations with bounce point "
                f"{sorted(table.loc[skipped, 'target_point'].unique())} "
                f"(converting only '{target_point}')."
            )
        table = table[~skipped]
    table = table.sort_values("epoch_seconds_UTC", kind="stable")
    doppler_frequency = table["transmitter_frequency_hz"].where(
        table["observable_type"] == DOPPLER_OBSERVABLE
    )
    table = table.assign(
        band=doppler_frequency.map(radar_frequency_band_string_from_hz, na_action="ignore").fillna(
            ""
        )
    )

    tracking_data = []
    keys = ["target_body", "transmitter", "receiver", "observable_type", "band"]
    for (target, transmitter, receiver, observable, band), group in table.groupby(keys, sort=False):
        data = TrackingData(
            observable_type=observable,
            link_ends=[
                ((station_body, transmitter), "transmitter"),
                ((target, ""), "reflector_1"),
                ((station_body, receiver), "receiver"),
            ],
            observations=[np.array([value]) for value in group["value"]],
            epochs=group["epoch_seconds_UTC"].tolist(),
            reference_link_end="receiver",
            time_scale="UTC",
        )
        data.set_observation_weights([np.array([sigma**-2]) for sigma in group["sigma"]])
        if band:
            data.add_string_vector_ancillary_setting("frequency bands", [band, band])
        tracking_data.append(data)
    return tracking_data, _frequency_supplementary_data(table, station_body)


def _frequency_supplementary_data(table: pd.DataFrame, station_body: str):
    doppler = table[table["observable_type"] == DOPPLER_OBSERVABLE]
    supplementary_data = []
    for transmitter, group in doppler.groupby("transmitter", sort=False):
        epochs = group["epoch_seconds_UTC"].to_numpy()
        frequencies = group["transmitter_frequency_hz"].to_numpy()
        changes = np.flatnonzero(np.diff(frequencies)) + 1
        starts = [epochs[0] - constants.JULIAN_DAY] + [
            0.5 * (epochs[index - 1] + epochs[index]) for index in changes
        ]
        history = {
            float(epoch): float(frequency)
            for epoch, frequency in zip(starts, frequencies[np.r_[0, changes]])
        }
        station_data = TrackingSupplementaryData(station_body, str(transmitter))
        station_data.set_frequency_supplementary_data(
            [PiecewiseConstantFrequencySupplementaryData(history)]
        )
        supplementary_data.append(station_data)
    return supplementary_data
