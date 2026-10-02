import tudatpy.data_input.tracking_data.obs_ADES.parsers as parser
from tudatpy.data_input.tracking_data.obs_ADES.parsers import MODE_DICT, CATALOG_DICT

import os
import pandas as pd
import pytest
from astropy.table import Table
import numpy as np

current_dir = os.getcwd()


# define the path of the file to be saved
@pytest.fixture
def file_path():
    return f"{current_dir}/obs_2.xml"


@pytest.fixture
def file_path_psv():
    return f"{current_dir}/obs_1.psv"


# define a dataframe with observations equal to file obs_1.psv
@pytest.fixture
def observations():
    return pd.DataFrame(
        {
            "permID": [1234456, 1234457, 1234458, 1234459],
            "provID": ["provID"] * 4,
            "trkSub": ["aa", "aa", "aa", "aa"],
            "artSat": [None, None, None, None],
            "mode": ["CCD", "CCD", "CCD", "CCD"],
            "stn": ["F51", "F51", "F51", "F51"],
            "obsTime": ["2016-08-29T12:32:34Z"] * 4,
            "ra": [0, 0.1, 0.1, 0.1],
            "dec": [90, 30.0, -30.0, 30.0],
            "rmsRA": [0.15, 0.15, 0.15, 0.15],
            "rmsDec": [0.13, 0.13, 0.13, 0.13],
            "rmsCorr": [0.21, 0.21, 0.21, 0.21],
            "astCat": ["2MASS", "2MASS", "2MASS", "2MASS"],
            "mag": [21.9, 21.9, 21.9, 21.9],
            "rmsMag": [0.25, 0.25, 0.25, 0.25],
            "band": ["w", "w", "w", "w"],
            "photCat": ["2MASS", "2MASS", "2MASS", "2MASS"],
            "logSNR": [0.775, 0.775, 0.775, 0.775],
            "notes": ["kmnt", "kmnt", "kmnt", "kmnt"],
            "remarks": ["High winds affected tracking"] * 4,
        }
    )


"""
Test all internal functions used by the parser 
"""


# test _first_present_column
def test_first_present_column(observations):
    result = parser._first_present_column(observations, "permID", "provID", "notes")
    # permID is the first present columns so should be returned
    pd.testing.assert_series_equal(result, observations["permID"])


def test_returns_none_when_no_names_are_present(observations):
    result = parser._first_present_column(observations, "missing1", "missing2")
    assert result is None


def test_error_if_no_names(observations):
    with pytest.raises(ValueError, match="At least one column name must be provided"):
        result = parser._first_present_column(observations)


# test _check_possible_failures
def test_rejects_missing_required_column(observations):
    invalid = observations.drop(columns=["mode"])
    with pytest.raises(ValueError, match="Missing required column: mode"):
        parser._check_possible_failures(invalid, "optical")


def test_rejects_missing_value_in_required_column(observations):
    invalid = observations.copy()
    invalid.loc[0, "ra"] = None
    with pytest.raises(ValueError, match="Column 'ra' contains missing values"):
        parser._check_possible_failures(invalid, "optical")


def test_rejects_row_without_identifier(observations):
    invalid = observations.copy()
    invalid.loc[0, ["permID", "provID", "artSat", "trkSub"]] = None
    with pytest.raises(ValueError, match="Each row must have at least one non-missing identifier"):
        parser._check_possible_failures(invalid, "optical")


def test_rejects_incomplete_location_group(observations):
    invalid = observations.copy()
    invalid["sys"] = ["ICRF"] * 4
    invalid["ctr"] = ["500"] * 4
    invalid["pos1"] = [1.0] * 4
    invalid["pos2"] = [2.0] * 4
    # pos3 is omitted, so the location group is incomplete
    with pytest.raises(
        ValueError, match="If any of sys, ctr, pos1, pos2, pos3 is provided, all must be provided"
    ):
        parser._check_possible_failures(invalid, "optical")


def test_rejects_position_that_is_not_numeric(observations):
    invalid = observations.copy()
    invalid["sys"] = ["ICRF"] * 4
    invalid["ctr"] = ["500"] * 4
    invalid["pos1"] = ["not-a-number"] * 4
    invalid["pos2"] = [2.0] * 4
    invalid["pos3"] = [3.0] * 4
    with pytest.raises(ValueError, match="pos1, pos2, and pos3 must be convertible to floats"):
        parser._check_possible_failures(invalid, "optical")


def test_rejects_ra_that_is_not_numeric(observations):
    invalid = observations.copy()
    invalid["ra"] = invalid["ra"].astype("object")
    invalid.loc[0, "ra"] = "not-a-number"
    with pytest.raises(ValueError, match="Column 'ra' must contain values convertible to floats"):
        parser._check_possible_failures(invalid, "optical")


def test_rejects_rms_value_that_is_not_numeric(observations):
    invalid = observations.copy()
    invalid["rmsRA"] = invalid["rmsRA"].astype("object")
    invalid.loc[0, "rmsRA"] = "not-a-number"
    with pytest.raises(ValueError, match="Column 'rmsRA' must be convertible to floats"):
        parser._check_possible_failures(invalid, "optical")


def test_rejects_invalid_observation_time(observations):
    invalid = observations.copy()
    invalid["obsTime"] = invalid["obsTime"].astype("object")
    invalid.loc[0, "obsTime"] = "not-a-time"
    with pytest.raises(ValueError, match="Column 'obsTime' must match yyyy-mm-ddThh:mm:ss"):
        parser._check_possible_failures(invalid, "optical")


# missing unit test for valid observatory code - could store the dataframe of the MPC locally to double check ??


def test_rejects_invalid_observation_kind(observations):
    with pytest.raises(ValueError, match="The input observation kind is invalid"):
        parser._check_possible_failures(observations.copy(), "unknown")


# check _epochs_UTC_to_seconds_UTC
def test_converts_observation_times_to_seconds_since_j2000(observations):
    result = parser._epochs_UTC_to_seconds_UTC(observations.copy())
    expected = (
        pd.to_datetime(observations["obsTime"], format="ISO8601", utc=True)
        - pd.Timestamp("2000-01-01T12:00:00Z")
    ).dt.total_seconds()
    pd.testing.assert_series_equal(
        result["epoch_seconds_UTC"],
        expected.rename("epoch_seconds_UTC"),
    )


def test_returns_dataframe_with_epoch_seconds_column(observations):
    result = parser._epochs_UTC_to_seconds_UTC(observations.copy())
    assert isinstance(result, pd.DataFrame)
    assert "epoch_seconds_UTC" in result.columns
    assert len(result) == len(observations)


# test _mode_to_note2
def test_maps_mode_to_note2(observations):
    result = parser._mode_to_note2(observations.copy())
    assert result["note2"].tolist() == ["C", "C", "C", "C"]


def test_returns_dataframe_with_note2_column(observations):
    result = parser._mode_to_note2(observations.copy())
    assert isinstance(result, pd.DataFrame)
    assert "note2" in result.columns
    assert len(result) == len(observations)


# test _astCat_to_catalog
def test_astCat_to_catalog(observations):
    result = parser._astCat_to_catalog(observations.copy())
    assert result["catalog"].tolist() == ["L", "L", "L", "L"]


def test_returns_dataframe_with_note2_column(observations):
    result = parser._astCat_to_catalog(observations.copy())
    assert isinstance(result, pd.DataFrame)
    assert "catalog" in result.columns
    assert len(result) == len(observations)


# test _roving_obs_columns
def test_adds_earth_roving_columns(observations):
    df = observations.copy()
    df["sys"] = ["WGS84"] * len(df)
    df["ctr"] = ["Earth"] * len(df)
    df["pos1"] = [1.0, 2.0, 3.0, 4.0]
    df["pos2"] = [5.0, 6.0, 7.0, 8.0]
    df["pos3"] = [9.0, 10.0, 11.0, 12.0]

    result_data = pd.DataFrame(index=range(len(df)))

    result = parser._roving_obs_columns(df, result_data)

    assert result["observatoryTipe"].tolist() == ["Earth Roving"] * len(df)
    assert result["roving_position_1"].tolist() == [1.0, 2.0, 3.0, 4.0]
    assert result["roving_position_2"].tolist() == [5.0, 6.0, 7.0, 8.0]
    assert result["roving_position_3"].tolist() == [9.0, 10.0, 11.0, 12.0]
    assert result["ctr"].tolist() == ["Earth"] * len(df)
    assert result["sys"].tolist() == ["WGS84"] * len(df)


def test_adds_spacecraft_roving_columns(observations):
    df = observations.copy()
    df["sys"] = ["ICRF_KM"] * len(df)
    df["pos1"] = [1.0] * len(df)
    df["pos2"] = [2.0] * len(df)
    df["pos3"] = [3.0] * len(df)
    df["vel1"] = [4.0] * len(df)
    df["vel2"] = [5.0] * len(df)
    df["vel3"] = [6.0] * len(df)

    result_data = pd.DataFrame(index=range(len(df)))

    result = parser._roving_obs_columns(df, result_data)

    assert result["observatoryTipe"].tolist() == ["Spacecraft"] * len(df)
    assert result["spacecraft_position_x"].tolist() == [1000.0] * len(df)
    assert result["spacecraft_position_y"].tolist() == [2000.0] * len(df)
    assert result["spacecraft_position_z"].tolist() == [3000.0] * len(df)
    assert result["spacecraft_velocity_vx"].tolist() == [4000.0] * len(df)
    assert result["spacecraft_velocity_vy"].tolist() == [5000.0] * len(df)
    assert result["spacecraft_velocity_vz"].tolist() == [6000.0] * len(df)


def test_rejects_none_result_data(observations):
    with pytest.raises(ValueError, match="The input result table is empty"):
        parser._roving_obs_columns(observations.copy(), None)


def test_rejects_different_number_of_rows(observations):
    result_data = pd.DataFrame(index=range(len(observations) - 1))
    with pytest.raises(ValueError, match="df and result_data must have the same number of rows"):
        parser._roving_obs_columns(observations.copy(), result_data)


# Final test on parse_ades_file
# as all used functions have been previously tested, use them in the file
def test_parses_optical_xml(file_path):
    result = parser.parse_ades_file(file_path)

    # check that it is an astropy table
    assert isinstance(result, Table)
    # check that there are 4 observations
    assert len(result) == 4

    assert result["number"][0] == "1234456"
    assert result["provisional_designation"].mask[0]
    assert result["epoch"][0] == "2016-08-29T12:32:34Z"
    assert result["observatory"][0] == "F51"
    assert result["mode"][0] == "CCD"
    assert result["note2"][0] == "C"

    # RA and DEC are converted from degrees to radians.
    assert result["RA"][1] == pytest.approx(np.deg2rad(0.1))
    assert result["DEC"][1] == pytest.approx(np.deg2rad(30.0))

    assert result["magnitude"][0] == pytest.approx(21.9)
    assert result["band"][0] == "w"
    assert result["rmsRA"][0] == pytest.approx(0.15)
    assert result["rmsDec"][0] == pytest.approx(0.13)
    assert result["rmsCorr"][0] == pytest.approx(0.21)

    # The epoch helper computes seconds from the J2000 UTC epoch.
    assert result["epoch_seconds_UTC"][0] == pytest.approx(
        (
            parser.pd.Timestamp("2016-08-29T12:32:34Z")
            - parser.pd.Timestamp("2000-01-01T12:00:00Z")
        ).total_seconds()
    )


# repeat the test giving a .psv file as input
def test_parses_optical_xml(file_path_psv):
    result = parser.parse_ades_file(file_path_psv)

    # check that it is an astropy table
    assert isinstance(result, Table)
    # check that there are 4 observations
    assert len(result) == 4

    assert result["number"][0] == "1234456"
    assert result["provisional_designation"].mask[0]
    assert result["epoch"][0] == "2016-08-29T12:32:34Z"
    assert result["observatory"][0] == "F51"
    assert result["mode"][0] == "CCD"
    assert result["note2"][0] == "C"

    # RA and DEC are converted from degrees to radians.
    assert result["RA"][1] == pytest.approx(np.deg2rad(0.1))
    assert result["DEC"][1] == pytest.approx(np.deg2rad(30.0))

    assert result["magnitude"][0] == pytest.approx(21.9)
    assert result["band"][0] == "w"
    assert result["rmsRA"][0] == pytest.approx(0.15)
    assert result["rmsDec"][0] == pytest.approx(0.13)
    assert result["rmsCorr"][0] == pytest.approx(0.21)

    # The epoch helper computes seconds from the J2000 UTC epoch.
    assert result["epoch_seconds_UTC"][0] == pytest.approx(
        (
            parser.pd.Timestamp("2016-08-29T12:32:34Z")
            - parser.pd.Timestamp("2000-01-01T12:00:00Z")
        ).total_seconds()
    )


def test_rejects_empty_file_path():
    with pytest.raises(ValueError, match="There is no input file"):
        parser.parse_ades_file("")


def test_rejects_unsupported_observation_type(tmp_path):
    path = tmp_path / "unsupported.xml"
    path.write_text(
        """\
<ades>
  <obsData>
    <unknownObservation><field>value</field></unknownObservation>
  </obsData>
</ades>
""",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="Unsupported observation type: unknownObservation",
    ):
        parser.parse_ades_file(str(path))
