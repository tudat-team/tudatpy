"""
Converts parsed CCSDS TDM messages into tudat-native `TrackingData` objects.
"""

import numpy as np
import pandas as pd

from tudatpy.astro.time_representation import DateTime
from tudatpy.data_input.tracking_data import (
    RampedFrequencySupplementaryData,
    TrackingData,
    TrackingSupplementaryData,
)

from .model.message import TDMMessage
from tudatpy import constants

_SPEED_OF_LIGHT = constants.SPEED_OF_LIGHT

# CCSDS TDM ANGLE_TYPE -> tudat's canonical ObservableType name
# (observation_models::getObservableName/getObservableType,
# createObservationCollection.cpp::getObservableTypeFromTrackingDataString).

# "XSYE (X: South, Y: East) and XEYN (X:East, Y:North)" have no equivalent tudat observable and are not supported yet.
_ANGLE_TYPE_TO_OBSERVABLE_NAME = {
    "RADEC": "AngularPosition",
    "AZEL": "AzimuthElevation",
}

## TDM scalar data keyword (RANGE, DOPPLER_INSTANTANEOUS, ETC...) -> tudat's canonical ObservableType name. Doppler
## is handled separately (_convert_doppler_group below), since (unlike
## RANGE) it needs a reference frequency before it means anything.
# _SCALAR_KEYWORD_TO_OBSERVABLE_NAME = {
#    "RANGE": "OneWayRange",
# }


def _tdm_epoch_to_utc_seconds(time_tag: str) -> float:
    """
    Parses a CCSDS TDM time tag (ISO-8601 or day-of-year) into UTC seconds
    since J2000, tudat's native epoch representation.
    """
    try:
        t = pd.to_datetime(time_tag)
    except ValueError:
        clean = time_tag.rstrip("Z")
        fmt = "%Y-%jT%H:%M:%S.%f" if "." in clean else "%Y-%jT%H:%M:%S"
        t = pd.to_datetime(clean, format=fmt)
    return DateTime(
        t.year,
        t.month,
        t.day,
        t.hour,
        t.minute,
        t.second + t.microsecond / 1e6 + t.nanosecond / 1e9,
    ).epoch()


def _resolve_link_ends(metadata) -> list:
    """
    Builds a `TrackingData`-shaped link_ends list from a TDM segment's
    PARTICIPANT_1/PARTICIPANT_2/PATH, matching the same "PATH starts with
    '1' -> PARTICIPANT_1 is the tracked object" convention already used by
    `TDMMessage.to_dataframe`. Only the two-participant case is supported
    -- `TDMMetadata` itself has no field for a third participant.
    """
    path = (metadata.PATH or "1,2").replace(" ", "")
    part_1, part_2 = metadata.PARTICIPANT_1, metadata.PARTICIPANT_2
    if path.startswith("1"):
        rso_id, station_name = part_1, part_2
    else:
        rso_id, station_name = part_2, part_1

    return rso_id, station_name


def _build_link_ends_list(rso_id: str, station_name: str) -> list:
    return [
        ((str(rso_id), ""), "transmitter"),
        (("Earth", str(station_name)), "receiver"),
    ]


############################################ FREQUENCY RAMPS ####################################
# def _build_transmit_frequency_ramp(
#    segment, sorted_time_tags: list[str]
# ) -> list[tuple[float, float, float]] | None:
#    """
#    Parses TRANSMIT_FREQ_1/TRANSMIT_FREQ_RATE_1 into a piecewise-linear
#    frequency ramp: a list of (start_time, start_frequency, rate) tuples,
#    sorted by start_time. Only band 1 is read -- TRANSMIT_FREQ_2..5 support
#    multi-band tracking setups this function doesn't handle yet.
#
#    Each entry's rate defaults to 0.0 (constant) if TRANSMIT_FREQ_RATE_1
#    isn't given at that same epoch.
#
#    Returns
#    -------
#    list[tuple[float, float, float]] | None
#        None if the segment has no TRANSMIT_FREQ_1 entries at all.
#    """
#    ramp_points = []
#    for time_tag in sorted_time_tags:
#        obs = segment.data[time_tag]
#        if "TRANSMIT_FREQ_1" not in obs:
#            continue
#        ramp_points.append(
#            (
#                _tdm_epoch_to_utc_seconds(time_tag),
#                float(obs["TRANSMIT_FREQ_1"]),
#                float(obs.get("TRANSMIT_FREQ_RATE_1", 0.0)),
#            )
#        )
#
#    if not ramp_points:
#        return None
#
#    ramp_points.sort(key=lambda entry: entry[0])
#    return ramp_points


# def _evaluate_frequency_ramp(ramp_points: list[tuple[float, float, float]], epoch: float) -> float:
#    """
#    Evaluates a piecewise-linear frequency ramp (as built by
#    `_build_transmit_frequency_ramp`) at `epoch`.
#
#    `epoch` before the first ramp point extrapolates backward using the
#    first point's rate; `epoch` after the last point extrapolates forward
#    using the last point's rate -- in both cases assuming that point's
#    ramp remained in effect, since there's nothing later/earlier to bound
#    it with.
#    """
#    applicable = ramp_points[0]
#    for point in ramp_points:
#        if point[0] > epoch:
#            break
#        applicable = point
#
#    start_time, start_frequency, rate = applicable
#    return start_frequency + rate * (epoch - start_time)


############################################ UNIT CONVERSIONS ####################################
def _convert_angle_group(
    segment, sorted_time_tags: list[str], link_ends: list
) -> "TrackingData | None":
    angle_type = segment.metadata.ANGLE_TYPE
    observable_name = _ANGLE_TYPE_TO_OBSERVABLE_NAME.get(angle_type)
    if observable_name is None:
        raise NotImplementedError(f"No tudat observable mapping for TDM ANGLE_TYPE={angle_type!r}.")

    epochs, observations, weights = [], [], []
    for time_tag in sorted_time_tags:
        obs = segment.data[time_tag]
        if "ANGLE_1" not in obs or "ANGLE_2" not in obs:
            continue
        epochs.append(_tdm_epoch_to_utc_seconds(time_tag))
        # Wrap ANGLE_1 (right ascension/azimuth) to (-pi, pi], same
        # convention used by optical_table_to_tracking_data for RA.
        # ANGLE_2 (declination/elevation) needs no wrapping: its CCSDS
        # range, [-90, 90] deg, already maps into [-pi/2, pi/2] rad.
        angle_1_rad = (np.radians(obs["ANGLE_1"]) + np.pi) % (2 * np.pi) - np.pi
        angle_2_rad = np.radians(obs["ANGLE_2"])
        observations.append(np.array([angle_1_rad, angle_2_rad]))
        weights_present = "SIGMA_ANGLE_1" in obs and "SIGMA_ANGLE_2" in obs
        if weights_present:
            sigma_1_rad = np.radians(obs["SIGMA_ANGLE_1"])
            sigma_2_rad = np.radians(obs["SIGMA_ANGLE_2"])
            weights.append(np.array([1.0 / sigma_1_rad**2, 1.0 / sigma_2_rad**2]))

    tracking_data = TrackingData(
        observable_name, link_ends, observations, epochs, "receiver", "UTC"
    )

    if len(weights) == len(observations):
        # Every epoch had SIGMA_ANGLE_1/SIGMA_ANGLE_2 -- set per-observation weights.
        # (This also covers the "no epoch had sigma" case: weights == observations == [].)
        tracking_data.set_observation_weights(weights)
    elif weights:
        # Some, but not all, epochs had sigma -- a partial weight vector would
        # silently misalign with the observations, so raise instead.
        raise ValueError(
            "Weights are present in the TDM, but their length is inconsistent with the number of observations."
        )
    return tracking_data


# def _convert_scalar_groups(segment, sorted_time_tags: list[str], link_ends: list) -> list:
#    tracking_data_objects = []
#    for keyword, observable_name in _SCALAR_KEYWORD_TO_OBSERVABLE_NAME.items():
#        if not any(keyword in obs for obs in segment.data.values()):
#            continue
#
#        epochs, observations = [], []
#        for time_tag in sorted_time_tags:
#            obs = segment.data[time_tag]
#            if keyword not in obs:
#                continue
#            epochs.append(_tdm_epoch_to_utc_seconds(time_tag))
#            # CCSDS TDM RANGE is conventionally given in km; tudat is SI (m).
#            observations.append(np.array([obs[keyword] * 1000.0]))
#
#        tracking_data_objects.append(
#            TrackingData(observable_name, link_ends, observations, epochs, "receiver", "UTC")
#        )
#    return tracking_data_objects


# def _convert_doppler_group(
#    segment, sorted_time_tags: list[str], link_ends: list, rso_id: str, station_name: str
# ):
#    """
#    Converts DOPPLER_INSTANTANEOUS into a "OneWayDoppler" `TrackingData`, in
#    velocity units (m/s) -- i.e. matching an observation-model configured
#    with `normalizeWithSpeedOfLight=False`, NOT tudat's dimensionless
#    default. This is the same physical quantity as a range-rate.
#
#    Requires TRANSMIT_FREQ_1 (optionally TRANSMIT_FREQ_RATE_1) to be
#    present in this segment -- CCSDS's raw Hz shift is meaningless without
#    knowing what frequency it's a shift from, and this function does not
#    accept an externally-supplied fallback. Raises if DOPPLER_INSTANTANEOUS
#    is present but TRANSMIT_FREQ_1 is not, rather than silently skipping,
#    since that's very likely lost data the caller would want to know about.
#
#    Caveats
#    -------
#    - Only DOPPLER_INSTANTANEOUS is handled; DOPPLER_INTEGRATED needs a
#      count-interval interpretation this function doesn't have and is not
#      converted.
#    - Only band 1 (TRANSMIT_FREQ_1/_RATE_1) is read; multi-band TDMs
#      (TRANSMIT_FREQ_2..5) are not supported.
#    - Assumes the standard convention DOPPLER_INSTANTANEOUS = f_received -
#      f_transmitted (unverified against a real reference TDM with known
#      truth Doppler -- worth checking before trusting this numerically).
#    - Uses the classical (non-relativistic) Doppler relation; tudat's own
#      OneWayDoppler model includes proper-time-rate correction terms this
#      conversion does not reproduce, so treat this as a first-order
#      approximation of what that model would predict.
#
#    Returns
#    -------
#    tuple[TrackingData, TrackingSupplementaryData] | tuple[None, None]
#        `(None, None)` if this segment has no DOPPLER_INSTANTANEOUS data.
#    """
#    if not any("DOPPLER_INSTANTANEOUS" in obs for obs in segment.data.values()):
#        return None, None
#
#    ramp_points = _build_transmit_frequency_ramp(segment, sorted_time_tags)
#    if ramp_points is None:
#        raise ValueError(
#            "Segment has DOPPLER_INSTANTANEOUS data but no TRANSMIT_FREQ_1 "
#            "entries to interpret it against -- cannot convert a raw Hz "
#            "Doppler shift without knowing the reference transmit "
#            "frequency. Supply a TDM with TRANSMIT_FREQ_1 (and, if the "
#            "uplink is ramped, TRANSMIT_FREQ_RATE_1)."
#        )
#
#    epochs, observations = [], []
#    for time_tag in sorted_time_tags:
#        obs = segment.data[time_tag]
#        if "DOPPLER_INSTANTANEOUS" not in obs:
#            continue
#        epoch = _tdm_epoch_to_utc_seconds(time_tag)
#        transmit_frequency_hz = _evaluate_frequency_ramp(ramp_points, epoch)
#        range_rate = -_SPEED_OF_LIGHT * obs["DOPPLER_INSTANTANEOUS"] / transmit_frequency_hz
#        epochs.append(epoch)
#        observations.append(np.array([range_rate]))
#
#    tracking_data = TrackingData(
#        "OneWayDoppler", link_ends, observations, epochs, "receiver", "UTC"
#    )
#
#    frequency_ramp_data = RampedFrequencySupplementaryData()
#    ramp_end_time = epochs[-1] if epochs else ramp_points[-1][0]
#    for i, (start_time, start_frequency, rate) in enumerate(ramp_points):
#        end_time = ramp_points[i + 1][0] if i + 1 < len(ramp_points) else ramp_end_time
#        frequency_ramp_data.add_frequency_ramp(start_time, end_time, start_frequency, rate)
#
#    supplementary_data = TrackingSupplementaryData("Earth", str(station_name))
#    supplementary_data.set_frequency_supplementary_data([frequency_ramp_data])
#
#    return tracking_data, supplementary_data


################################## CONVERSION TO TRACKING DATA ####################################
def tdm_message_to_tracking_data(tdm_message: TDMMessage) -> tuple[list, list]:
    """
    Converts every segment of a parsed TDM message into one or more
    `TrackingData` objects: one per distinct observable found in that
    segment's data. A segment carrying ANGLE_1/ANGLE_2 becomes one
    `AngularPosition`/`AzimuthElevation` `TrackingData`; a segment that
    also carries RANGE and/or DOPPLER_INSTANTANEOUS becomes additional,
    separate `TrackingData` objects for the same link ends.

    Per-observation angle uncertainty (`SIGMA_ANGLE_1`/`SIGMA_ANGLE_2` --
    see `implementations.tdm.TudatTDMDataKeywordExtension`) is attached via
    `TrackingData.set_observation_weights` (as `1 / sigma**2`, in rad^-2)
    when present at every epoch of the group; partial coverage (present at
    some epochs but not others) is treated as absent, since a partial
    weight vector would silently misalign with the observations.

    See `_convert_doppler_group` for the Doppler conversion's caveats
    (velocity-unit convention, classical approximation, TRANSMIT_FREQ_1
    only, DOPPLER_INTEGRATED unsupported).

    Parameters
    ----------
    tdm_message : TDMMessage
        An already-parsed TDM message.

    Returns
    -------
    tuple[list[TrackingData], list[TrackingSupplementaryData]]
        One `TrackingData` per (segment, observable) group found (empty if
        none); one `TrackingSupplementaryData` per segment that needed a
        transmit-frequency ramp to convert Doppler (i.e. usually empty
        unless DOPPLER_INSTANTANEOUS was present and convertible).
    """
    tracking_data_objects = []
    supplementary_data_objects = []

    for segment in tdm_message.segments:
        rso_id, station_name = _resolve_link_ends(segment.metadata)
        link_ends = _build_link_ends_list(rso_id, station_name)
        sorted_time_tags = sorted(segment.data.keys())

        has_angles = any("ANGLE_1" in obs and "ANGLE_2" in obs for obs in segment.data.values())
        if has_angles:
            tracking_data_objects.append(_convert_angle_group(segment, sorted_time_tags, link_ends))

        # UNCOMMENT WHENEVER WE WANT TO SUPPORT SCALAR GROUPS.
        # tracking_data_objects.extend(_convert_scalar_groups(segment, sorted_time_tags, link_ends))

        # UNCOMMENT WHENEVER WE WANT TO SUPPORT DOPPLER GROUP
        # doppler_tracking_data, doppler_supplementary_data = _convert_doppler_group(
        #    segment, sorted_time_tags, link_ends, rso_id, station_name
        # )
        # if doppler_tracking_data is not None:
        #    tracking_data_objects.append(doppler_tracking_data)
        #    supplementary_data_objects.append(doppler_supplementary_data)

    return tracking_data_objects, supplementary_data_objects
