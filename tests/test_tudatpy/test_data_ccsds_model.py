import numpy as np
import pandas as pd
import pytest

from tudatpy.data_input.tracking_data.ccsds import ccsds_utils
from tudatpy.data_input.tracking_data.ccsds.implementations.common import CCSDSHeader
from tudatpy.data_input.tracking_data.ccsds.implementations.omm import OMMData, OMMMetadata
from tudatpy.data_input.tracking_data.ccsds.implementations.oem import OEMMetadata, OEMStateVector
from tudatpy.data_input.tracking_data.ccsds.implementations.tdm import TDMMetadata
from tudatpy.data_input.tracking_data.ccsds.model.message import (
    OEMMessage,
    OMMMessage,
    Segment,
    TDMMessage,
)
from tudatpy.data_input.tracking_data.ccsds.model.parser import NDMParser
from tudatpy.data_input.tracking_data.ccsds.model.writer import NDMWriter
from tudatpy.data_input.tracking_data.ccsds import tracking_data_conversion as tdc


def _header(**overrides) -> CCSDSHeader:
    kwargs = dict(CCSDS_VERS="2.0", ORIGINATOR="TEST", CREATION_DATE="2025-01-01T00:00:00Z")
    kwargs.update(overrides)
    return CCSDSHeader(**kwargs)


def _tdm_metadata(**overrides) -> TDMMetadata:
    kwargs = dict(
        TIME_SYSTEM="UTC",
        PARTICIPANT_1="SAT1",
        PARTICIPANT_2="STATION1",
        PATH="1,2",
        ANGLE_TYPE="RADEC",
    )
    kwargs.update(overrides)
    return TDMMetadata(**kwargs)


def _oem_metadata(**overrides) -> OEMMetadata:
    kwargs = dict(
        OBJECT_NAME="SAT1",
        OBJECT_ID="2025-001A",
        CENTER_NAME="EARTH",
        REF_FRAME="EME2000",
        TIME_SYSTEM="UTC",
        START_TIME="2025-01-01T00:00:00Z",
        STOP_TIME="2025-01-01T00:01:00Z",
    )
    kwargs.update(overrides)
    return OEMMetadata(**kwargs)


def _omm_metadata(**overrides) -> OMMMetadata:
    kwargs = dict(
        OBJECT_NAME="SAT1",
        OBJECT_ID="2025-001A",
        CENTER_NAME="EARTH",
        REF_FRAME="TEME",
        TIME_SYSTEM="UTC",
        MEAN_ELEMENT_THEORY="SGP4",
    )
    kwargs.update(overrides)
    return OMMMetadata(**kwargs)


def _omm_data(**overrides) -> OMMData:
    kwargs = dict(
        EPOCH="2025-001T00:00:00Z",
        MEAN_MOTION=15.5,
        ECCENTRICITY=0.001,
        INCLINATION=51.6,
        RA_OF_ASC_NODE=100.0,
        ARG_OF_PERICENTER=90.0,
        MEAN_ANOMALY=0.0,
    )
    kwargs.update(overrides)
    return OMMData(**kwargs)


############################### ccsds_utils ###############################
def test_covariance_matrix_from_named_keywords_returns_none_when_absent():
    assert ccsds_utils.covariance_matrix_from_named_keywords({"OBJECT_NAME": "SAT1"}) is None


def test_covariance_matrix_from_named_keywords_raises_when_partial():
    with pytest.raises(ValueError, match="all-or-nothing"):
        ccsds_utils.covariance_matrix_from_named_keywords({"CX_X": "1.0", "CY_X": "0.1"})


def test_covariance_matrix_from_named_keywords_builds_symmetric_matrix():
    data = {
        key: str(i) for i, (key, _, _) in enumerate(ccsds_utils._NAMED_COVARIANCE_KEYWORD_CELLS)
    }
    matrix = ccsds_utils.covariance_matrix_from_named_keywords(data)
    assert matrix.shape == (6, 6)
    assert np.array_equal(matrix, matrix.T)
    # CY_X (index 1) sits at (row=1, col=0)
    assert matrix[1, 0] == 1.0
    assert matrix[0, 1] == 1.0


def test_covariance_matrix_to_named_keywords_round_trips_with_from_named_keywords():
    original = np.diag([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).astype(float)
    original[1, 0] = original[0, 1] = 0.5
    named = ccsds_utils.covariance_matrix_to_named_keywords(original)
    rebuilt = ccsds_utils.covariance_matrix_from_named_keywords(named)
    assert np.allclose(rebuilt, original)


def test_reconstruct_covariance_builds_symmetric_lower_triangular_matrix():
    lines = [
        "1.0",
        "0.1 2.0",
        "0.2 0.3 3.0",
        "0.4 0.5 0.6 4.0",
        "0.7 0.8 0.9 1.0 5.0",
        "1.1 1.2 1.3 1.4 1.5 6.0",
    ]
    matrix = ccsds_utils.reconstruct_covariance(lines)
    assert matrix.shape == (6, 6)
    assert np.array_equal(matrix, matrix.T)
    assert matrix[0, 0] == 1.0
    assert matrix[5, 5] == 6.0
    assert matrix[2, 1] == matrix[1, 2] == 0.3


def test_make_lower_triangular_round_trips_with_reconstruct_covariance():
    original = np.arange(36, dtype=float).reshape(6, 6)
    original = (original + original.T) / 2.0  # force symmetry
    lines = ccsds_utils.make_lower_triangular(original)
    assert len(lines) == 6
    assert len(lines[0].split()) == 1
    assert len(lines[5].split()) == 6
    rebuilt = ccsds_utils.reconstruct_covariance(lines)
    assert np.allclose(rebuilt, original)


############################### model.parser: NDMParser ###############################
@pytest.mark.parametrize(
    "filename,expected_type",
    [
        ("foo.tdm", TDMMessage),
        ("foo.TDM", TDMMessage),
        ("foo.oem", OEMMessage),
        ("foo.omm", OMMMessage),
    ],
)
def test_get_message_type_from_extension(filename, expected_type):
    assert NDMParser()._get_message_type(filename) is expected_type


def test_get_message_type_raises_for_unsupported_extension():
    with pytest.raises(ValueError, match="Unsupported file type"):
        NDMParser()._get_message_type("foo.txt")


def test_normalize_version_keyword_renames_type_specific_key():
    normalized = NDMParser._normalize_version_keyword(
        {"CCSDS_TDM_VERS": "1.0", "ORIGINATOR": "TEST"}
    )
    assert normalized == {"CCSDS_VERS": "1.0", "ORIGINATOR": "TEST"}


def test_parse_tdm_basic(tmp_path):
    content = (
        "CCSDS_TDM_VERS = 2.0\n"
        "CREATION_DATE = 2025-01-01T00:00:00Z\n"
        "ORIGINATOR = TEST\n\n"
        "META_START\n"
        "TIME_SYSTEM = UTC\n"
        "PARTICIPANT_1 = SAT1\n"
        "PARTICIPANT_2 = STATION1\n"
        "PATH = 1,2\n"
        "ANGLE_TYPE = RADEC\n"
        "META_STOP\n\n"
        "DATA_START\n"
        "ANGLE_1 = 2025-001T00:00:00Z 10.0\n"
        "ANGLE_2 = 2025-001T00:00:00Z 20.0\n"
        "DATA_STOP\n"
    )
    path = tmp_path / "test.tdm"
    path.write_text(content)

    msg = NDMParser().parse(path)

    assert isinstance(msg, TDMMessage)
    assert msg.header.CCSDS_VERS == "2.0"
    assert len(msg.segments) == 1
    assert msg.segments[0].metadata.PARTICIPANT_1 == "SAT1"
    assert msg.segments[0].data == {"2025-001T00:00:00Z": {"ANGLE_1": 10.0, "ANGLE_2": 20.0}}


def test_parse_tdm_multi_segment(tmp_path):
    content = (
        "CCSDS_TDM_VERS = 2.0\n"
        "CREATION_DATE = 2025-01-01T00:00:00Z\n\n"
        "META_START\n"
        "TIME_SYSTEM = UTC\n"
        "PARTICIPANT_1 = SAT1\n"
        "PARTICIPANT_2 = STATION1\n"
        "PATH = 1,2\n"
        "ANGLE_TYPE = RADEC\n"
        "META_STOP\n\n"
        "DATA_START\n"
        "ANGLE_1 = 2025-001T00:00:00Z 10.0\n"
        "ANGLE_2 = 2025-001T00:00:00Z 20.0\n"
        "DATA_STOP\n\n"
        "META_START\n"
        "TIME_SYSTEM = UTC\n"
        "PARTICIPANT_1 = SAT1\n"
        "PARTICIPANT_2 = STATION2\n"
        "PATH = 1,2\n"
        "ANGLE_TYPE = AZEL\n"
        "META_STOP\n\n"
        "DATA_START\n"
        "ANGLE_1 = 2025-001T00:02:00Z 30.0\n"
        "ANGLE_2 = 2025-001T00:02:00Z 40.0\n"
        "DATA_STOP\n"
    )
    path = tmp_path / "multi.tdm"
    path.write_text(content)

    msg = NDMParser().parse(path)

    assert len(msg.segments) == 2
    assert msg.segments[0].metadata.PARTICIPANT_2 == "STATION1"
    assert msg.segments[1].metadata.PARTICIPANT_2 == "STATION2"
    assert msg.segments[1].data["2025-001T00:02:00Z"]["ANGLE_1"] == 30.0


def test_parse_oem_with_covariance(tmp_path):
    content = (
        "CCSDS_OEM_VERS = 2.0\n"
        "CREATION_DATE = 2025-01-01T00:00:00Z\n\n"
        "META_START\n"
        "OBJECT_NAME = SAT1\n"
        "OBJECT_ID = 2025-001A\n"
        "CENTER_NAME = EARTH\n"
        "REF_FRAME = EME2000\n"
        "TIME_SYSTEM = UTC\n"
        "START_TIME = 2025-01-01T00:00:00Z\n"
        "STOP_TIME = 2025-01-01T00:01:00Z\n"
        "META_STOP\n\n"
        "2025-01-01T00:00:00Z 1000.0 2000.0 3000.0 1.0 2.0 3.0\n"
        "2025-01-01T00:01:00Z 1001.0 2001.0 3001.0 1.1 2.1 3.1\n\n"
        "COVARIANCE_START\n"
        "EPOCH = 2025-01-01T00:00:00Z\n"
        "1.0\n"
        "0.1 2.0\n"
        "0.2 0.3 3.0\n"
        "0.4 0.5 0.6 4.0\n"
        "0.7 0.8 0.9 1.0 5.0\n"
        "1.1 1.2 1.3 1.4 1.5 6.0\n"
        "COVARIANCE_STOP\n"
    )
    path = tmp_path / "test.oem"
    path.write_text(content)

    msg = NDMParser().parse(path)

    assert isinstance(msg, OEMMessage)
    assert len(msg.segments) == 1
    state = msg.segments[0].data["2025-01-01T00:00:00Z"]
    assert (state.x, state.y, state.z) == (1000.0, 2000.0, 3000.0)
    assert msg.segments[0].covariance is not None
    cov = msg.segments[0].covariance["2025-01-01T00:00:00Z"]
    assert cov[0, 0] == 1.0
    assert cov[5, 5] == 6.0


def test_parse_oem_covariance_without_epoch_raises(tmp_path):
    content = (
        "CCSDS_OEM_VERS = 2.0\n"
        "CREATION_DATE = 2025-01-01T00:00:00Z\n\n"
        "META_START\n"
        "OBJECT_NAME = SAT1\n"
        "OBJECT_ID = 2025-001A\n"
        "CENTER_NAME = EARTH\n"
        "REF_FRAME = EME2000\n"
        "TIME_SYSTEM = UTC\n"
        "START_TIME = 2025-01-01T00:00:00Z\n"
        "STOP_TIME = 2025-01-01T00:01:00Z\n"
        "META_STOP\n\n"
        "2025-01-01T00:00:00Z 1000.0 2000.0 3000.0 1.0 2.0 3.0\n\n"
        "COVARIANCE_START\n"
        "1.0\n"
        "COVARIANCE_STOP\n"
    )
    path = tmp_path / "bad.oem"
    path.write_text(content)

    with pytest.raises(ValueError, match="EPOCH"):
        NDMParser().parse(path)


def test_parse_omm_with_covariance_and_user_defined(tmp_path):
    content = (
        "CCSDS_OMM_VERS = 2.0\n"
        "CREATION_DATE = 2025-01-01T00:00:00Z\n"
        "OBJECT_NAME = SAT1\n"
        "OBJECT_ID = 2025-001A\n"
        "CENTER_NAME = EARTH\n"
        "REF_FRAME = TEME\n"
        "TIME_SYSTEM = UTC\n"
        "MEAN_ELEMENT_THEORY = SGP4\n"
        "EPOCH = 2025-001T00:00:00Z\n"
        "MEAN_MOTION = 15.5\n"
        "ECCENTRICITY = 0.001\n"
        "INCLINATION = 51.6\n"
        "RA_OF_ASC_NODE = 100.0\n"
        "ARG_OF_PERICENTER = 90.0\n"
        "MEAN_ANOMALY = 0.0\n"
        "BSTAR = 0.0001\n"
        "CX_X = 1.0\n"
        "CY_X = 0.1\n"
        "CY_Y = 2.0\n"
        "CZ_X = 0.2\n"
        "CZ_Y = 0.3\n"
        "CZ_Z = 3.0\n"
        "CX_DOT_X = 0.4\n"
        "CX_DOT_Y = 0.5\n"
        "CX_DOT_Z = 0.6\n"
        "CX_DOT_X_DOT = 4.0\n"
        "CY_DOT_X = 0.7\n"
        "CY_DOT_Y = 0.8\n"
        "CY_DOT_Z = 0.9\n"
        "CY_DOT_X_DOT = 1.0\n"
        "CY_DOT_Y_DOT = 5.0\n"
        "CZ_DOT_X = 1.1\n"
        "CZ_DOT_Y = 1.2\n"
        "CZ_DOT_Z = 1.3\n"
        "CZ_DOT_X_DOT = 1.4\n"
        "CZ_DOT_Y_DOT = 1.5\n"
        "CZ_DOT_Z_DOT = 6.0\n"
        "USER_DEFINED_TEST = hello\n"
    )
    path = tmp_path / "test.omm"
    path.write_text(content)

    msg = NDMParser().parse(path)

    assert isinstance(msg, OMMMessage)
    assert msg.metadata.MEAN_ELEMENT_THEORY == "SGP4"
    assert msg.data.BSTAR == 0.0001
    assert msg.covariance is not None
    cov = msg.covariance["2025-001T00:00:00Z"]
    assert cov[0, 0] == 1.0
    assert cov[5, 5] == 6.0
    assert msg.user_defined == {"USER_DEFINED_TEST": "hello"}


############################### model.message ###############################
def test_segment_tdm_validation_raises_when_angle_type_missing():
    metadata = _tdm_metadata(ANGLE_TYPE=None)
    data = {"2025-001T00:00:00Z": {"ANGLE_1": 10.0, "ANGLE_2": 20.0}}
    with pytest.raises(ValueError, match="ANGLE_TYPE"):
        Segment(metadata, data)


def test_segment_omm_validation_raises_when_theory_requirement_missing():
    metadata = _omm_metadata(MEAN_ELEMENT_THEORY="SGP4")
    data = _omm_data()  # no BSTAR
    with pytest.raises(ValueError, match="BSTAR"):
        Segment(metadata, data)


def test_oem_message_has_covariance():
    seg = Segment(_oem_metadata(), {})
    msg = OEMMessage(header=_header(), segments=[seg])
    assert msg.has_covariance() is False

    seg.covariance = {"2025-01-01T00:00:00Z": np.eye(6)}
    assert msg.has_covariance() is True


def test_oem_message_time_bounds_across_segments():
    seg1 = Segment(
        _oem_metadata(),
        {"2025-01-01T00:00:00Z": OEMStateVector(x=1, y=1, z=1, vx=0, vy=0, vz=0)},
    )
    seg2 = Segment(
        _oem_metadata(),
        {"2025-01-01T00:05:00Z": OEMStateVector(x=2, y=2, z=2, vx=0, vy=0, vz=0)},
    )
    msg = OEMMessage(header=_header(), segments=[seg1, seg2])

    assert msg.time_bounds() == ["2025-01-01T00:00:00Z", "2025-01-01T00:05:00Z"]
    assert msg.time_bounds(segment_idx=1) == ["2025-01-01T00:05:00Z", "2025-01-01T00:05:00Z"]


def test_oem_message_to_ccsds_state_history_concatenates_and_sorts():
    seg1 = Segment(
        _oem_metadata(),
        {"2025-01-01T00:01:00Z": OEMStateVector(x=1, y=2, z=3, vx=4, vy=5, vz=6)},
    )
    seg2 = Segment(
        _oem_metadata(),
        {"2025-01-01T00:00:00Z": OEMStateVector(x=7, y=8, z=9, vx=10, vy=11, vz=12)},
    )
    msg = OEMMessage(header=_header(), segments=[seg1, seg2])

    history = msg.to_ccsds_state_history()
    epochs = list(history.keys())

    assert epochs == sorted(epochs)
    assert np.array_equal(history[epochs[0]], np.array([7, 8, 9, 10, 11, 12]))
    assert np.array_equal(history[epochs[1]], np.array([1, 2, 3, 4, 5, 6]))


def test_oem_message_to_ccsds_covariance_history_empty_when_none():
    seg = Segment(_oem_metadata(), {})
    msg = OEMMessage(header=_header(), segments=[seg])
    assert msg.to_ccsds_covariance_history() == {}

    seg.covariance = {"2025-01-01T00:00:00Z": np.eye(6) * 2.0}
    history = msg.to_ccsds_covariance_history()
    assert list(history.values())[0][0, 0] == 2.0


def test_tdm_message_number_of_segments_and_empty_time_bounds():
    msg = TDMMessage(header=_header(), segments=[])
    assert msg.number_of_segments() == 0
    assert msg.time_bounds() == ["N/A", "N/A"]


def test_tdm_message_time_bounds_actual_data_vs_metadata():
    metadata = _tdm_metadata(START_TIME="2025-01-01T00:00:00Z", STOP_TIME="2025-01-01T00:10:00Z")
    data = {"2025-01-01T00:05:00Z": {"ANGLE_1": 1.0, "ANGLE_2": 2.0}}
    seg = Segment(metadata, data)
    msg = TDMMessage(header=_header(), segments=[seg])

    assert msg.time_bounds(actual_data=True) == ["2025-01-01T00:05:00Z", "2025-01-01T00:05:00Z"]
    assert msg.time_bounds(actual_data=False) == ["2025-01-01T00:00:00Z", "2025-01-01T00:10:00Z"]


# def test_tdm_message_from_dataframe_wraps_ra_and_builds_metadata():
#    df = pd.DataFrame(
#        {
#            "Time": pd.to_datetime(["2025-01-01T00:00:00Z", "2025-01-01T00:01:00Z"]),
#            "Right_Ascension": [-10.0, 30.0],
#            "Declination": [5.0, 6.0],
#            "Sigma_Right_Ascension": [0.1, 0.1],
#            "Sigma_Declination": [0.2, 0.2],
#        }
#    )
#    msg = TDMMessage.from_dataframe(df, norad_id="25544", obs_code="DELFT")
#
#    seg = msg.segments[0]
#    assert seg.metadata.PARTICIPANT_1 == "25544"
#    assert seg.metadata.PARTICIPANT_2 == "DELFT"
#    assert seg.metadata.PATH == "1,2"
#    assert seg.metadata.ANGLE_TYPE == "RADEC"
#    assert seg.metadata.REFERENCE_FRAME == "ICRF"
#
#    first_epoch = min(seg.data.keys())
#    assert seg.data[first_epoch]["ANGLE_1"] == pytest.approx(350.0)  # -10 % 360
#    assert seg.data[first_epoch]["SIGMA_ANGLE_1"] == pytest.approx(0.1)
#    assert seg.data[first_epoch]["SIGMA_ANGLE_2"] == pytest.approx(0.2)


# def test_tdm_message_from_dataframe_raises_for_unsupported_angle_type():
#    df = pd.DataFrame(
#        {"Time": pd.to_datetime(["2025-01-01T00:00:00Z"]), "Right_Ascension": [1.0], "Declination": [2.0]}
#    )
#    with pytest.raises(ValueError, match="Unsupported angle type"):
#        TDMMessage.from_dataframe(df, norad_id="25544", obs_code="DELFT", angle_type="BOGUS")


# def test_tdm_message_to_dataframe_participant_mapping_follows_path():
#    metadata = _tdm_metadata(PARTICIPANT_1="STATION1", PARTICIPANT_2="SAT1", PATH="2,1")
#    data = {"2025-01-01T00:00:00Z": {"ANGLE_1": 10.0, "ANGLE_2": 20.0}}
#    seg = Segment(metadata, data)
#    msg = TDMMessage(header=_header(), segments=[seg])
#
#    df = msg.to_dataframe()
#
#    assert df.iloc[0]["PARTICIPANT_1"] == "SAT1"  # PATH "2,1" -> participant 2 is the RSO
#    assert df.iloc[0]["PARTICIPANT_2"] == "STATION1"
#    assert df.iloc[0]["RA"] == 10.0
#    assert df.iloc[0]["DEC"] == 20.0
#    assert "epoch_seconds_UTC" in df.columns
#    assert "epoch_seconds_TDB" in df.columns
#    # TDB trails UTC by the ~69s (2025-era) UTC-TDB offset, not the reverse.
#    assert df.iloc[0]["epoch_seconds_TDB"] > df.iloc[0]["epoch_seconds_UTC"]


# def test_tdm_message_to_dataframe_parses_doy_timestamps():
#    metadata = _tdm_metadata()
#    data = {"2025-001T00:00:00.000000": {"ANGLE_1": 10.0, "ANGLE_2": 20.0}}
#    seg = Segment(metadata, data)
#    msg = TDMMessage(header=_header(), segments=[seg])
#
#    df = msg.to_dataframe()
#
#    assert len(df) == 1
#    assert df.index[0] == pd.Timestamp("2025-01-01")


# def test_tdm_message_to_dataframe_empty_when_no_segment_data():
#    seg = Segment(_tdm_metadata(), {})
#    msg = TDMMessage(header=_header(), segments=[seg])
#    assert msg.to_dataframe().empty


def test_omm_message_segment_properties():
    metadata = _omm_metadata()
    data = _omm_data(BSTAR=0.0001)
    seg = Segment(metadata, data)
    seg.covariance = {data.EPOCH: np.eye(6)}
    seg.user_defined = {"USER_DEFINED_TEST": "hello"}
    msg = OMMMessage(header=_header(), segments=[seg])

    assert msg.metadata is metadata
    assert msg.data is data
    assert msg.covariance == {data.EPOCH: seg.covariance[data.EPOCH]}
    assert msg.user_defined == {"USER_DEFINED_TEST": "hello"}


############################### model.writer: NDMWriter ###############################
def test_validate_before_write_raises_when_stop_before_start():
    metadata = _tdm_metadata(START_TIME="2025-01-02T00:00:00Z", STOP_TIME="2025-01-01T00:00:00Z")
    seg = Segment(metadata, {})
    msg = TDMMessage(header=_header(), segments=[seg])
    with pytest.raises(ValueError, match="STOP_TIME is earlier than START_TIME"):
        NDMWriter().validate_before_write(msg)


def test_write_and_parse_tdm_round_trip(tmp_path):
    metadata = _tdm_metadata(START_TIME="2025-01-01T00:00:00Z", STOP_TIME="2025-01-01T00:01:00Z")
    data = {"2025-01-01T00:00:00Z": {"ANGLE_1": 350.0, "ANGLE_2": 5.0}}
    seg = Segment(metadata, data)
    msg = TDMMessage(header=_header(), segments=[seg])

    path = tmp_path / "out.tdm"
    NDMWriter().write(msg, str(path))
    parsed = NDMParser().parse(path)

    assert parsed.segments[0].metadata.PARTICIPANT_1 == "SAT1"
    assert parsed.segments[0].metadata.ANGLE_TYPE == "RADEC"
    assert parsed.segments[0].data == data


def test_write_and_parse_oem_round_trip_with_meters_conversion(tmp_path):
    metadata = _oem_metadata()
    data = {
        "2025-01-01T00:00:00Z": OEMStateVector(
            x=1.0e6, y=2.0e6, z=3.0e6, vx=1.0e3, vy=2.0e3, vz=3.0e3
        )
    }
    seg = Segment(metadata, data)
    covariance = np.diag([1e6, 1e6, 1e6, 1.0, 1.0, 1.0]).astype(float)
    seg.covariance = {"2025-01-01T00:00:00Z": covariance}
    msg = OEMMessage(header=_header(), segments=[seg])

    path = tmp_path / "out.oem"
    NDMWriter().write(msg, str(path), input_message_in_meters=True)
    parsed = NDMParser().parse(path)

    state = parsed.segments[0].data["2025-01-01T00:00:00Z"]
    assert (state.x, state.y, state.z) == pytest.approx((1000.0, 2000.0, 3000.0))
    assert (state.vx, state.vy, state.vz) == pytest.approx((1.0, 2.0, 3.0))

    parsed_cov = parsed.segments[0].covariance["2025-01-01T00:00:00Z"]
    assert np.diag(parsed_cov) == pytest.approx([1.0, 1.0, 1.0, 1e-6, 1e-6, 1e-6])


def test_write_and_parse_omm_round_trip_with_covariance_and_user_defined(tmp_path):
    metadata = _omm_metadata()
    data = _omm_data(BSTAR=0.0001)
    seg = Segment(metadata, data)
    covariance = np.diag([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).astype(float)
    seg.covariance = {data.EPOCH: covariance}
    seg.user_defined = {"USER_DEFINED_TEST": "hello"}
    msg = OMMMessage(header=_header(), segments=[seg])

    path = tmp_path / "out.omm"
    NDMWriter().write(msg, str(path))
    parsed = NDMParser().parse(path)

    assert parsed.data.BSTAR == pytest.approx(0.0001)
    assert np.diag(parsed.covariance[data.EPOCH]) == pytest.approx([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    assert parsed.user_defined == {"USER_DEFINED_TEST": "hello"}


############################### tracking_data_conversion ###############################
def test_tdm_epoch_to_utc_seconds_iso_and_doy_agree():
    iso_seconds = tdc._tdm_epoch_to_utc_seconds("2025-01-01T00:00:00Z")
    doy_seconds = tdc._tdm_epoch_to_utc_seconds("2025-001T00:00:00")
    assert iso_seconds == pytest.approx(doy_seconds)


def test_resolve_link_ends_path_variants():
    metadata = _tdm_metadata(PARTICIPANT_1="SAT1", PARTICIPANT_2="STATION1", PATH="1,2")
    assert tdc._resolve_link_ends(metadata) == ("SAT1", "STATION1")

    metadata_swapped = _tdm_metadata(PARTICIPANT_1="STATION1", PARTICIPANT_2="SAT1", PATH="2,1")
    assert tdc._resolve_link_ends(metadata_swapped) == ("SAT1", "STATION1")

    metadata_default = _tdm_metadata(PARTICIPANT_1="SAT1", PARTICIPANT_2="STATION1", PATH=None)
    assert tdc._resolve_link_ends(metadata_default) == ("SAT1", "STATION1")


def test_build_link_ends_list_shape():
    link_ends = tdc._build_link_ends_list("SAT1", "STATION1")
    assert link_ends == [
        (("SAT1", ""), "transmitter"),
        (("Earth", "STATION1"), "receiver"),
    ]


def test_tdm_message_to_tracking_data_radec():
    metadata = _tdm_metadata(ANGLE_TYPE="RADEC")
    data = {
        "2025-001T00:00:00Z": {"ANGLE_1": 10.0, "ANGLE_2": 20.0},
        "2025-001T00:01:00Z": {"ANGLE_1": 11.0, "ANGLE_2": 21.0},
    }
    seg = Segment(metadata, data)
    msg = TDMMessage(header=_header(), segments=[seg])

    tracking_data_objects, supplementary_data_objects = tdc.tdm_message_to_tracking_data(msg)

    assert supplementary_data_objects == []
    assert len(tracking_data_objects) == 1
    td = tracking_data_objects[0]
    assert td.observable_type == "AngularPosition"
    assert td.reference_link_end == "receiver"
    assert td.time_scale == "UTC"
    assert td.link_ends == [(("SAT1", ""), "transmitter"), (("Earth", "STATION1"), "receiver")]
    assert len(td.observations) == 2
    assert len(td.epochs) == 2


def test_tdm_message_to_tracking_data_azel():
    metadata = _tdm_metadata(ANGLE_TYPE="AZEL")
    data = {"2025-001T00:00:00Z": {"ANGLE_1": 10.0, "ANGLE_2": 20.0}}
    seg = Segment(metadata, data)
    msg = TDMMessage(header=_header(), segments=[seg])

    tracking_data_objects, _ = tdc.tdm_message_to_tracking_data(msg)

    assert tracking_data_objects[0].observable_type == "AzimuthElevation"


def test_tdm_message_to_tracking_data_unsupported_angle_type_raises():
    metadata = _tdm_metadata(ANGLE_TYPE="XSYE")
    data = {"2025-001T00:00:00Z": {"ANGLE_1": 10.0, "ANGLE_2": 20.0}}
    seg = Segment(metadata, data)
    msg = TDMMessage(header=_header(), segments=[seg])

    with pytest.raises(NotImplementedError, match="XSYE"):
        tdc.tdm_message_to_tracking_data(msg)


def test_tdm_message_to_tracking_data_wraps_angle_1_to_minus_pi_pi():
    metadata = _tdm_metadata(ANGLE_TYPE="RADEC")
    data = {"2025-001T00:00:00Z": {"ANGLE_1": 200.0, "ANGLE_2": 20.0}}
    seg = Segment(metadata, data)
    msg = TDMMessage(header=_header(), segments=[seg])

    tracking_data_objects, _ = tdc.tdm_message_to_tracking_data(msg)

    observation = tracking_data_objects[0].observations[0]
    assert observation[0] == pytest.approx(np.radians(200.0) - 2 * np.pi)
    assert observation[1] == pytest.approx(np.radians(20.0))


def test_tdm_message_to_tracking_data_sets_weights_only_when_sigma_present_everywhere():
    metadata = _tdm_metadata(ANGLE_TYPE="RADEC")

    # Sigma present at only one of two epochs: partial coverage, weights not set.
    partial_data = {
        "2025-001T00:00:00Z": {
            "ANGLE_1": 10.0,
            "ANGLE_2": 20.0,
            "SIGMA_ANGLE_1": 1.0,
            "SIGMA_ANGLE_2": 1.0,
        },
        "2025-001T00:01:00Z": {"ANGLE_1": 11.0, "ANGLE_2": 21.0},
    }
    partial_msg = TDMMessage(header=_header(), segments=[Segment(metadata, partial_data)])
    partial_td, _ = tdc.tdm_message_to_tracking_data(partial_msg)
    assert list(partial_td[0].get_observation_weights()) == []

    # Sigma present at every epoch: weights set as 1/sigma^2 (rad^-2).
    full_data = {
        "2025-001T00:00:00Z": {
            "ANGLE_1": 10.0,
            "ANGLE_2": 20.0,
            "SIGMA_ANGLE_1": 1.0,
            "SIGMA_ANGLE_2": 2.0,
        },
    }
    full_msg = TDMMessage(header=_header(), segments=[Segment(metadata, full_data)])
    full_td, _ = tdc.tdm_message_to_tracking_data(full_msg)
    weights = full_td[0].get_observation_weights()
    assert weights[0][0] == pytest.approx(1.0 / np.radians(1.0) ** 2)
    assert weights[0][1] == pytest.approx(1.0 / np.radians(2.0) ** 2)


def test_tdm_message_to_tracking_data_empty_when_no_angle_data():
    metadata = _tdm_metadata(ANGLE_TYPE=None)
    msg = TDMMessage(header=_header(), segments=[Segment(metadata, {})])

    tracking_data_objects, supplementary_data_objects = tdc.tdm_message_to_tracking_data(msg)

    assert tracking_data_objects == []
    assert supplementary_data_objects == []
