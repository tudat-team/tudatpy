import tudatpy.data_input.tracking_data.obs_ADES.parsers as parser

import os
import pandas as pd
import pytest

current_dir = os.getcwd()

file_path = f"{current_dir}/obs_1.psv"


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
