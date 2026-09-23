"""Canonical radar tracking-data table and its conversion to Tudat tracking data.

The table has one row per observation: epochs in UTC seconds since J2000, ranges as
round-trip light time times the speed of light [m], Doppler as received frequency [Hz].
"""

import numpy as np
import pandas as pd

from tudatpy import constants
from tudatpy.data_input.tracking_data import (
    PiecewiseConstantFrequencySupplementaryData,
    TrackingData,
    TrackingSupplementaryData,
)
from tudatpy.data_input.tracking_data.optical_utilities import datetime_to_utc_seconds
from tudatpy.dynamics import environment
from tudatpy.estimation.observations_setup.ancillary_settings import FrequencyBands

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
    """Return an empty canonical radar table."""
    return pd.DataFrame(columns=RADAR_COLUMNS)


def validate_radar_data(table: pd.DataFrame) -> pd.DataFrame:
    """Check a canonical radar table and return it with numeric columns as floats.

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
        value=(
            raw.loc[is_doppler, "transmitter_frequency_hz"]
            + raw.loc[is_doppler, "doppler_hz"]
        ),
        sigma=raw.loc[is_doppler, "doppler_sigma_hz"],
    )
    if range_rows.empty and doppler_rows.empty:
        return empty_radar_table()
    return validate_radar_data(
        pd.concat([range_rows, doppler_rows], ignore_index=True)[RADAR_COLUMNS]
    )


def radar_data_from_table(table) -> pd.DataFrame:
    """Return the radar table stored in a parsed MPC 80-column table."""
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

    Epoch bounds may be UTC seconds since J2000 or datetime-like objects. Station
    filters match the transmitter or receiver; bare numeric codes are padded to
    three characters.
    """
    keep = pd.Series(True, index=table.index)
    for column, allowed in [
        ("target_body", target_body),
        ("target_point", target_point),
        ("observable_type", observable_type),
    ]:
        if allowed is not None:
            keep &= table[column].isin(
                [allowed] if isinstance(allowed, str) else list(allowed)
            )

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

    Tudat defines S, X, Ku and Ka bands. For a passive reflector the band only
    selects a turnaround ratio, which is one for every band pair after calling
    :func:`set_reflector_turnaround_ratio`.
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
    """
    if target_point is not None:
        table = table[table["target_point"] == target_point]
    table = table.sort_values("epoch_seconds_UTC", kind="stable")
    doppler_frequency = table["transmitter_frequency_hz"].where(
        table["observable_type"] == DOPPLER_OBSERVABLE
    )
    table = table.assign(
        band=doppler_frequency.map(
            radar_frequency_band_string_from_hz, na_action="ignore"
        ).fillna("")
    )

    tracking_data = []
    keys = ["target_body", "transmitter", "receiver", "observable_type", "band"]
    for (target, transmitter, receiver, observable, band), group in table.groupby(
        keys, sort=False
    ):
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
        data.set_observation_weights(
            [np.array([sigma**-2]) for sigma in group["sigma"]]
        )
        if band:
            data.add_string_vector_ancillary_setting("frequency bands", [band, band])
        tracking_data.append(data)
    return tracking_data, _frequency_supplementary_data(table, station_body)


def set_reflector_turnaround_ratio(bodies, target_body: str) -> None:
    """Set a passive radar target's turnaround ratio to one for all band pairs."""
    body = bodies.get_body(target_body)
    systems = (
        body.system_models
        if body.system_models is not None
        else environment.VehicleSystems()
    )
    bands = [
        FrequencyBands.s_band,
        FrequencyBands.x_band,
        FrequencyBands.ku_band,
        FrequencyBands.ka_band,
    ]
    systems.set_transponder_turnaround_ratio(
        {(uplink, downlink): 1.0 for uplink in bands for downlink in bands}
    )
    body.system_models = systems


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
