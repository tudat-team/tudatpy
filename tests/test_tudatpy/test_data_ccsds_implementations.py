import pytest
from pydantic import ValidationError

from tudatpy.data_input.tracking_data.ccsds.implementations import common, oem, tdm, omm


############### ccsds.implementations ##############
def test_oem_state_vector_access_and_required_fields():
    """Dict-style access matches attribute access; a missing field raises."""
    state_vector = oem.OEMStateVector(x=1.0, y=2.0, z=3.0, vx=0.1, vy=0.2, vz=0.3)
    assert state_vector["x"] == state_vector.x == 1.0
    assert state_vector["vz"] == state_vector.vz == 0.3

    with pytest.raises(ValidationError, match="vz"):
        oem.OEMStateVector(x=1.0, y=2.0, z=3.0, vx=0.1, vy=0.2)


def test_oem_metadata_optional_fields_default_to_none():
    oem_metadata = oem.OEMMetadata(
        OBJECT_NAME="International Space Station",
        OBJECT_ID="1998-067A",
        CENTER_NAME="EARTH",
        REF_FRAME="EME2000",
        TIME_SYSTEM="UTC",
        START_TIME="2025-001T00:00:00Z",
        STOP_TIME="2025-002T00:00:00Z",
    )
    assert oem_metadata.INTERPOLATION is None
    assert oem_metadata.INTERPOLATION_DEGREE is None
    assert oem_metadata.REF_FRAME_EPOCH is None
    assert oem_metadata.COMMENT is None


def test_oem_metadata_interpolation_degree_required_with_interpolation():
    base_kwargs = dict(
        OBJECT_NAME="ISS",
        OBJECT_ID="1998-067A",
        CENTER_NAME="EARTH",
        REF_FRAME="EME2000",
        TIME_SYSTEM="UTC",
        START_TIME="2025-001T00:00:00Z",
        STOP_TIME="2025-002T00:00:00Z",
    )
    with pytest.raises(ValidationError, match="INTERPOLATION_DEGREE"):
        oem.OEMMetadata(**base_kwargs, INTERPOLATION="LAGRANGE")

    # Both together, or neither, are fine.
    oem.OEMMetadata(**base_kwargs, INTERPOLATION="LAGRANGE", INTERPOLATION_DEGREE=7)
    oem.OEMMetadata(**base_kwargs)


def test_tdm_metadata_optional_fields_default_to_none():
    tdm_metadata = tdm.TDMMetadata(
        TIME_SYSTEM="UTC",
        PARTICIPANT_1="TU_Delft_Rooftop",
        PARTICIPANT_2="1998-067A",
    )
    assert tdm_metadata.MODE is None
    assert tdm_metadata.PATH is None
    assert tdm_metadata.ANGLE_TYPE is None
    assert tdm_metadata.REFERENCE_FRAME is None
    assert tdm_metadata.COMMENT is None


def test_tdm_metadata_time_system_has_no_silent_default():
    """
    TIME_SYSTEM and MODE were deliberately given no default: silently
    assuming e.g. UTC when the source data is actually TAI would corrupt
    downstream analysis, so omitting TIME_SYSTEM must raise rather than
    fall back to a guessed value.
    """
    with pytest.raises(ValidationError, match="TIME_SYSTEM"):
        tdm.TDMMetadata(PARTICIPANT_1="TU_Delft_Rooftop", PARTICIPANT_2="1998-067A")


def test_tdm_data_keyword_values_match_ccsds_keywords():
    assert tdm.TDMDataKeyword.ANGLE_1.value == "ANGLE_1"
    assert tdm.TDMDataKeyword.RANGE.value == "RANGE"
    assert tdm.TDMDataKeyword.DOPPLER_INSTANTANEOUS.value == "DOPPLER_INSTANTANEOUS"
    assert not hasattr(tdm.TDMDataKeyword, "SIGMA_ANGLE_1")


def test_tudat_tdm_data_keyword_extension_is_separate_from_canonical_enum():
    assert tdm.TudatTDMDataKeywordExtension.SIGMA_ANGLE_1.value == "SIGMA_ANGLE_1"
    assert tdm.TudatTDMDataKeywordExtension.SIGMA_ANGLE_2.value == "SIGMA_ANGLE_2"
    assert not hasattr(tdm.TDMDataKeyword, "SIGMA_ANGLE_2")


def test_omm_metadata_optional_fields_default_to_none():
    omm_metadata = omm.OMMMetadata(
        OBJECT_NAME="International Space Station",
        OBJECT_ID="1998-067A",
        CENTER_NAME="EARTH",
        REF_FRAME="TEME",
        TIME_SYSTEM="UTC",
        MEAN_ELEMENT_THEORY="SGP4",
    )
    assert omm_metadata.REF_FRAME_EPOCH is None
    assert omm_metadata.COMMENT is None


def test_omm_data_requires_orbit_size_field():
    """OMMData._check_orbit_size_field: needs SEMI_MAJOR_AXIS or MEAN_MOTION, not neither."""
    base_kwargs = dict(
        EPOCH="2025-001T00:00:00Z",
        ECCENTRICITY=0.001,
        INCLINATION=51.6,
        RA_OF_ASC_NODE=100.0,
        ARG_OF_PERICENTER=90.0,
        MEAN_ANOMALY=0.0,
    )
    with pytest.raises(ValidationError, match="SEMI_MAJOR_AXIS or MEAN_MOTION"):
        omm.OMMData(**base_kwargs)

    omm.OMMData(**base_kwargs, SEMI_MAJOR_AXIS=7000.0)
    omm.OMMData(**base_kwargs, MEAN_MOTION=15.5)


def _make_omm_metadata(mean_element_theory: str) -> omm.OMMMetadata:
    return omm.OMMMetadata(
        OBJECT_NAME="ISS",
        OBJECT_ID="1998-067A",
        CENTER_NAME="EARTH",
        REF_FRAME="TEME",
        TIME_SYSTEM="UTC",
        MEAN_ELEMENT_THEORY=mean_element_theory,
    )


def _make_omm_data(**extra_fields) -> omm.OMMData:
    return omm.OMMData(
        EPOCH="2025-001T00:00:00Z",
        MEAN_MOTION=15.5,
        ECCENTRICITY=0.001,
        INCLINATION=51.6,
        RA_OF_ASC_NODE=100.0,
        ARG_OF_PERICENTER=90.0,
        MEAN_ANOMALY=0.0,
        **extra_fields,
    )


@pytest.mark.parametrize(
    "mean_element_theory,required_fields",
    [
        ("SGP4", ("BSTAR",)),
        ("SGP4-XP", ("BTERM", "AGOM")),
        ("SGP", ("MEAN_MOTION_DOT", "MEAN_MOTION_DDOT")),
        ("PPT3", ("MEAN_MOTION_DOT", "MEAN_MOTION_DDOT")),
    ],
)
def test_omm_validate_data_for_theory_raises_when_missing(mean_element_theory, required_fields):
    metadata = _make_omm_metadata(mean_element_theory)
    data = _make_omm_data()
    with pytest.raises(ValueError) as exc_info:
        omm.validate_omm_data_for_theory(metadata, data)
    for field in required_fields:
        assert field in str(exc_info.value)


def test_omm_validate_data_for_theory_passes_when_satisfied():
    metadata = _make_omm_metadata("SGP4")
    data = _make_omm_data(BSTAR=0.0001)
    omm.validate_omm_data_for_theory(metadata, data)  # should not raise


def test_omm_validate_data_for_theory_ignores_unlisted_theory():
    """DSST/USM/etc. aren't in the requirements registry, so nothing is enforced for them."""
    metadata = _make_omm_metadata("DSST")
    data = _make_omm_data()
    omm.validate_omm_data_for_theory(metadata, data)  # should not raise


def test_omm_keyword_sets_derived_from_model_fields():
    assert "CCSDS_OMM_VERS" in omm.OMM_HEADER_KEYWORDS
    assert "CCSDS_VERS" not in omm.OMM_HEADER_KEYWORDS
    assert "COMMENT" not in omm.OMM_HEADER_KEYWORDS
    assert "ORIGINATOR" in omm.OMM_HEADER_KEYWORDS

    assert "MEAN_ELEMENT_THEORY" in omm.OMM_METADATA_KEYWORDS
    assert "COMMENT" not in omm.OMM_METADATA_KEYWORDS


@pytest.mark.parametrize(
    "creation_date",
    ["2025-01-15T12:00:00Z", "2025-015T12:00:00.5", "2025-01-15T12:00:00"],
)
def test_ccsds_header_accepts_calendar_and_doy_dates(creation_date):
    header = common.CCSDSHeader(CREATION_DATE=creation_date)
    assert header.CREATION_DATE == creation_date
    assert header.ORIGINATOR == "Tudat"


def test_ccsds_header_rejects_malformed_date():
    with pytest.raises(ValidationError, match="Invalid CCSDS date format"):
        common.CCSDSHeader(CREATION_DATE="not-a-date")
