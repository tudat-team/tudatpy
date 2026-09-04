from enum import Enum

from pydantic import BaseModel

# Metadata keywords that become mandatory once a given data keyword appears
# anywhere in a segment's data section. Requirements are additive per data
# keyword rather than exclusive per segment, since a single TDM segment may
# legally mix observable types (e.g. angles and range in the same pass).
# Enforced by `validate_tdm_metadata_for_data` below, not by `TDMMetadata`
# itself, since requiredness here depends on the data section, not on the
# metadata alone.
#
# Only ANGLE_1/ANGLE_2 -> ANGLE_TYPE is populated: CCSDS 503.0-B-2 3.5.4.1
# ties ANGLE_TYPE directly to interpreting ANGLE_1/ANGLE_2 ("The ANGLE_TYPE
# metadata keyword indicates how these two keywords should be
# interpreted") -- without it the angle values are ambiguous.
#
# RANGE and DOPPLER_INTEGRATED are deliberately NOT added here, even though
# they're the obvious next candidates -- the Blue Book uses much softer
# language for them than it does for ANGLE_TYPE, so treating them as hard
# requirements (raising like ANGLE_TYPE does) would reject TDMs the
# standard itself considers valid:
#   - RANGE (3.5.2.7): "The 'RANGE_UNITS' metadata keyword should always
#     be specified, but if it is not, the default (preferred) value shall
#     be 'km'." -- a recommendation with a defined fallback ('km'), not a
#     mandatory field.
#   - DOPPLER_INTEGRATED (3.5.2.3): interpretation depends on
#     TIMETAG_REF/INTEGRATION_REF/INTEGRATION_INTERVAL, but Table 3-3
#     marks all three optional with no stated conditional-on-
#     DOPPLER_INTEGRATED requirement.
#   - DOPPLER_INSTANTANEOUS (3.5.2.2): transmit/receive frequency "should
#     be supplied" for ionosphere/solar-plasma corrections -- again a
#     recommendation, and those are DATA keywords, not metadata, so they
#     don't fit this table's shape anyway.
#
# If RANGE/Doppler support is re-enabled in tracking_data_conversion.py
# and stricter-than-Blue-Book validation is wanted (e.g. always requiring
# RANGE_UNITS rather than silently assuming 'km'), that would need to be a
# deliberate project policy decision, not an implicit side effect of
# extending this table -- see the RANGE_UNITS default-handling note above
# before adding entries here.
_REQUIRED_METADATA_FOR_DATA_KEYWORD: dict[str, tuple[str, ...]] = {
    "ANGLE_1": ("ANGLE_TYPE",),
    "ANGLE_2": ("ANGLE_TYPE",),
}


class TDMDataKeyword(str, Enum):
    """
    Canonical CCSDS TDM data-section keywords (503.0-B-2 Table 3-5).

    Unlike OMM's data section (one fixed set of named fields per epoch --
    see `OMMData`), a TDM's data section has no single fixed schema: which
    of these keywords appear, and in what combination, varies per file
    depending on what's actually being tracked (e.g. an angles pass uses
    only ANGLE_1/ANGLE_2; a Doppler pass uses different keywords again).
    That's why TDM has no equivalent of `OMMData`/`OEMStateVector` wrapping
    `Segment.data` -- it stays a plain `dict[str, dict[str, float]]`.

    This enum exists purely as a reference/validation aid (e.g. for future
    keyword -> observable_type mapping); it is not enforced anywhere --
    `NDMParser`/`NDMWriter` accept any keyword.
    """

    ANGLE_1 = "ANGLE_1"
    ANGLE_2 = "ANGLE_2"
    CARRIER_POWER = "CARRIER_POWER"
    CLOCK_BIAS = "CLOCK_BIAS"
    CLOCK_DRIFT = "CLOCK_DRIFT"
    DOPPLER_INSTANTANEOUS = "DOPPLER_INSTANTANEOUS"
    DOPPLER_INTEGRATED = "DOPPLER_INTEGRATED"
    DOR = "DOR"
    PC_N0 = "PC_N0"
    PR_N0 = "PR_N0"
    PRESSURE = "PRESSURE"
    RANGE = "RANGE"
    RECEIVE_FREQ = "RECEIVE_FREQ"
    RECEIVE_FREQ_1 = "RECEIVE_FREQ_1"
    RECEIVE_FREQ_2 = "RECEIVE_FREQ_2"
    RECEIVE_FREQ_3 = "RECEIVE_FREQ_3"
    RECEIVE_FREQ_4 = "RECEIVE_FREQ_4"
    RECEIVE_FREQ_5 = "RECEIVE_FREQ_5"
    RHUMIDITY = "RHUMIDITY"
    STEC = "STEC"
    TEMPERATURE = "TEMPERATURE"
    TRANSMIT_FREQ_1 = "TRANSMIT_FREQ_1"
    TRANSMIT_FREQ_2 = "TRANSMIT_FREQ_2"
    TRANSMIT_FREQ_3 = "TRANSMIT_FREQ_3"
    TRANSMIT_FREQ_4 = "TRANSMIT_FREQ_4"
    TRANSMIT_FREQ_5 = "TRANSMIT_FREQ_5"
    TRANSMIT_FREQ_RATE_1 = "TRANSMIT_FREQ_RATE_1"
    TRANSMIT_FREQ_RATE_2 = "TRANSMIT_FREQ_RATE_2"
    TRANSMIT_FREQ_RATE_3 = "TRANSMIT_FREQ_RATE_3"
    TRANSMIT_FREQ_RATE_4 = "TRANSMIT_FREQ_RATE_4"
    TRANSMIT_FREQ_RATE_5 = "TRANSMIT_FREQ_RATE_5"
    TROPO_DRY = "TROPO_DRY"
    TROPO_WET = "TROPO_WET"
    VLBI_DELAY = "VLBI_DELAY"


class TudatTDMDataKeywordExtension(str, Enum):
    """
    Non-standard TDM data-section keywords used by Tudat to represent
    quantities the CCSDS TDM Blue Book (503.0-B-2 Table 3-5) has no keyword
    for, e.g. notably, `SIGMA_ANGLE_1`/`SIGMA_ANGLE_2`: CCSDS TDM has
    no uncertainty/covariance representation at all (unlike OEM/OMM, which
    have proper covariance blocks), so per-observation angle uncertainty is
    written as its own time-tagged keyword, at the same epoch as the
    corresponding `ANGLE_1`/`ANGLE_2` value, in the same units.

    NOTE the SIGMA fields are Read/written exactly like any
    `TDMDataKeyword` (`NDMParser`/`NDMWriter`make no distinction),
    but external CCSDS TDM consumers won't recognize these keywords.
    """

    SIGMA_ANGLE_1 = "SIGMA_ANGLE_1"
    SIGMA_ANGLE_2 = "SIGMA_ANGLE_2"


class TDMMetadata(BaseModel):
    """
    Metadata fields specific to Tracking Data Messages (TDM).

    Attributes:
        TIME_SYSTEM: The time system used for the message (e.g., UTC, TAI).
        START_TIME: The start time of the data in the message.
        STOP_TIME: The stop time of the data in the message.
        PARTICIPANT_1: The first participant in the tracking session.
        PARTICIPANT_2: The second participant in the tracking session.
        MODE: The tracking mode (e.g., SEQUENTIAL, SINGLE_DIFF).
        PATH: The signal path (e.g., 1,2).
        ANGLE_TYPE: The type of angles provided (e.g., AZEL, RADEC).
        REFERENCE_FRAME: The reference frame for the data.
        COMMENT: Optional comments regarding the metadata.
    """

    TIME_SYSTEM: str
    START_TIME: str | None = None
    STOP_TIME: str | None = None
    PARTICIPANT_1: str
    PARTICIPANT_2: str | None = None
    MODE: str | None = None
    PATH: str | None = None
    ANGLE_TYPE: str | None = None
    REFERENCE_FRAME: str | None = None
    COMMENT: str | None = None


def validate_tdm_metadata_for_data(metadata: TDMMetadata, data: dict) -> None:
    """
    Cross-checks a TDM segment's metadata against the data keywords
    actually present, raising if a keyword required (per
    `_REQUIRED_METADATA_FOR_DATA_KEYWORD`) by an observable present in
    `data` is missing from `metadata`.

    Parameters
    ----------
    metadata : TDMMetadata
        The segment's metadata.
    data : dict
        The segment's data section: `{time_tag: {keyword: value}}`.

    Raises
    ------
    ValueError
        If a required metadata field is missing for a data keyword that
        is present.
    """
    present_keywords = {keyword for epoch_data in data.values() for keyword in epoch_data}

    missing: dict[str, list[str]] = {}
    for data_keyword in present_keywords:
        for required_field in _REQUIRED_METADATA_FOR_DATA_KEYWORD.get(data_keyword, ()):
            if getattr(metadata, required_field, None) is None:
                missing.setdefault(required_field, []).append(data_keyword)

    if missing:
        details = "; ".join(
            f"{field} (required by {', '.join(sorted(keywords))})"
            for field, keywords in sorted(missing.items())
        )
        raise ValueError(f"TDM segment metadata is missing required field(s): {details}")
