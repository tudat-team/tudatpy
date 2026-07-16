import importlib
from unittest.mock import patch

import pandas as pd
import pytest
from tudatpy.data_input.tracking_data.mpc import (
    BatchMPC,
    read_mpc_data,
)
from tudatpy.data_input.tracking_data.optical_utilities import (
    create_augmented_optical_table,
    read_astropy_optical_data,
    read_pandas_optical_data,
    standardize_optical_dataframe,
    validate_optical_table,
)
import numpy as np
from astropy.table import Table
from tudatpy.data_input.tracking_data import obs_80_cols

read_80_column_data = obs_80_cols.read_80_column_data


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(autouse=True)
def fail_on_implicit_station_info():
    """BatchMPC construction and simple optical loading must not retrieve station data."""

    def unexpected_get_observatory_codes():
        raise AssertionError("Unexpected MPC observatory-code retrieval")

    with patch("astroquery.mpc.MPC.get_observatory_codes", unexpected_get_observatory_codes):
        yield


# A minimal optical table. RA/DEC are in degrees, epoch is Julian Day.
@pytest.fixture
def valid_obs_table():
    """Minimal valid table satisfying _req_cols, RA/DEC in degrees."""
    return pd.DataFrame(
        {
            "number": ["433", "433"],
            "epoch": [2451545.0, 2451546.0],  # JD
            "RA": [10.0, 200.0],
            "DEC": [5.0, -5.0],
            "band": ["V", "V"],
            "observatory": ["500", "673"],
        }
    )


def _batch_from_optical_table(table, in_degrees=True, custom_name=None):
    batch = BatchMPC()
    batch._table = create_augmented_optical_table(
        table, in_degrees=in_degrees, custom_name=custom_name
    )
    batch._refresh_metadata()
    return batch


# ---------------------------------------------------------------------------
# standardize_optical_dataframe: normalizes observatory codes and MPC numbers
# before they're stored, so downstream lookups/joins are consistent.
# ---------------------------------------------------------------------------
class TestStandardizeDataframe:
    def test_zero_pads_short_observatory_codes(self):
        # MPC observatory codes are always 3 chars (e.g. "89" -> "089");
        # alphanumeric codes like "C51" are already 3 chars and pass through.
        df = pd.DataFrame({"observatory": ["89", "673", "C51"], "number": [1, 2, 3]})
        out = standardize_optical_dataframe(df)
        assert out["observatory"].tolist() == ["089", "673", "C51"]

    def test_strips_whitespace_and_stringifies_number(self):
        df = pd.DataFrame({"number": [" 433 ", 1]})
        out = standardize_optical_dataframe(df)
        assert out["number"].tolist() == ["433", "1"]

    def test_does_not_mutate_input(self):
        # standardize_optical_dataframe must work on a copy; callers (e.g.
        # get_observations) rely on the original frame staying untouched.
        df = pd.DataFrame({"observatory": ["89"]})
        standardize_optical_dataframe(df)
        assert df["observatory"].iloc[0] == "89"  # original untouched


# ---------------------------------------------------------------------------
# validate_optical_table: the guard rail shared by optical data ingestion that
# checks reference frame, required columns, and non-emptiness before any
# data is ingested.
# ---------------------------------------------------------------------------
class TestValidateTable:
    def test_rejects_non_j2000_frame(self, valid_obs_table):
        # Only J2000 is currently supported (see class docstring notes).
        with pytest.raises(NotImplementedError):
            validate_optical_table(valid_obs_table, frame="FK5")

    def test_rejects_missing_required_columns(self):
        df = pd.DataFrame({"number": [1], "epoch": [2451545.0]})  # missing RA/DEC/band/observatory
        with pytest.raises(ValueError, match="mandatory columns"):
            validate_optical_table(df, frame="J2000")

    def test_rejects_empty_table(self, valid_obs_table):
        with pytest.raises(ValueError, match="zero rows"):
            validate_optical_table(valid_obs_table.iloc[0:0], frame="J2000")

    def test_rejects_unsupported_type(self):
        # validate_optical_table only understands (Q)Table and DataFrame inputs.
        with pytest.raises(TypeError):
            validate_optical_table([1, 2, 3], frame="J2000")


# ---------------------------------------------------------------------------
# RA/DEC unit handling: degrees must be
# converted to radians and wrapped into [-pi, pi]; radians must be passed
# through unchanged when in_degrees=False.
# ---------------------------------------------------------------------------
class TestAngleConversion:
    def test_degrees_converted_and_wrapped_to_minus_pi_pi(self, valid_obs_table):
        table = create_augmented_optical_table(valid_obs_table, in_degrees=True)
        ra = table["RA"].to_numpy()
        assert (ra >= -np.pi).all() and (ra <= np.pi).all()
        # 200 deg -> should wrap to a negative radian value, not stay > pi
        assert ra[1] < 0

    def test_radians_left_untouched_when_in_degrees_false(self, valid_obs_table):
        table_rad = valid_obs_table.copy()
        table_rad["RA"] = np.radians(table_rad["RA"])
        table_rad["DEC"] = np.radians(table_rad["DEC"])
        table = create_augmented_optical_table(table_rad, in_degrees=False)
        assert np.isclose(table["RA"].iloc[0], table_rad["RA"].iloc[0])


# ---------------------------------------------------------------------------
# Augmented table: the shared intermediate representation for optical data.
# Epochs remain UTC here; conversion to TDB is done later when creating an
# ObservationCollection from TrackingData.
# ---------------------------------------------------------------------------
def test_augmented_table_adds_utc_epochs_only(valid_obs_table):
    out = create_augmented_optical_table(valid_obs_table)
    assert "epoch_seconds_UTC" in out.columns
    assert "epoch_seconds_TDB" not in out.columns


def _check_tracking_data_pair(result, expected_sets):
    tracking_data, supplementary_data = result
    assert len(tracking_data) == expected_sets
    assert supplementary_data == []
    assert all(data.time_scale == "UTC" for data in tracking_data)
    assert all(data.observable_type == "AngularPosition" for data in tracking_data)


def test_read_pandas_optical_data_returns_tracking_data(valid_obs_table):
    _check_tracking_data_pair(read_pandas_optical_data(valid_obs_table), expected_sets=2)


def test_read_pandas_optical_data_marks_vfcc17_weighing_scheme(valid_obs_table):
    tracking_data, supplementary_data = read_pandas_optical_data(valid_obs_table, add_weights=True)
    assert supplementary_data == []
    assert all(data.weighing_scheme == "VFCC17" for data in tracking_data)


def test_read_astropy_optical_data_returns_tracking_data(valid_obs_table):
    _check_tracking_data_pair(
        read_astropy_optical_data(Table.from_pandas(valid_obs_table)), expected_sets=2
    )


# Fakes astroquery's MPC.get_observations(code) return value: the real call
# returns an astropy-Table-like object exposing .to_pandas().
def _fake_mpc_obs(df):
    class _Wrap:
        def to_pandas(self_inner):
            return df

    return _Wrap()


def _mpc_observations(df):
    return patch("astroquery.mpc.MPC.get_observations", return_value=_fake_mpc_obs(df))


# ---------------------------------------------------------------------------
# get_observations: identifier-resolution logic. This is the trickiest part
# of the class - for each MPC code it must decide the object's canonical
# "number" from whatever combination of number/desig/comet_type the MPC
# response contains. Each test below pins one specific branch of that logic.
# ---------------------------------------------------------------------------
class TestGetObservationsIdentifierResolution:
    def test_read_mpc_data_returns_tracking_data(self):
        df = pd.DataFrame(
            {
                "number": ["433"],
                "epoch": [2451545.0],
                "RA": [1.0],
                "DEC": [1.0],
                "band": ["V"],
                "observatory": ["500"],
                "note2": [""],
                "desig": [None],
            }
        )
        with _mpc_observations(df):
            _check_tracking_data_pair(read_mpc_data([433]), expected_sets=1)

    def test_numbered_asteroid_uses_number(self):
        # Simple case: a numbered asteroid reports its number directly.
        df = pd.DataFrame(
            {
                "number": ["433"],
                "epoch": [2451545.0],
                "RA": [1.0],
                "DEC": [1.0],
                "band": ["V"],
                "observatory": ["500"],
                "note2": [""],
                "desig": [None],
            }
        )
        with _mpc_observations(df):
            batch = BatchMPC()
            batch.get_observations([433])
        assert batch.table["number"].iloc[0] == "433"

    def test_comet_combines_number_and_comet_type(self):
        # Comets/interstellar objects report a "comet_type" alongside their
        # number; the two must be concatenated (e.g. "3" + "I" -> "3I").
        df = pd.DataFrame(
            {
                "number": ["3"],
                "epoch": [2451545.0],
                "RA": [1.0],
                "DEC": [1.0],
                "band": ["V"],
                "observatory": ["500"],
                "note2": [""],
                "desig": [None],
                "comet_type": ["I"],
            }
        )
        with _mpc_observations(df):
            batch = BatchMPC()
            batch.get_observations(["3I"])
        assert batch.table["number"].iloc[0] == "3I"

    def test_no_number_falls_back_to_designation(self):
        # Not-yet-numbered objects have no "number" but do have a
        # provisional designation, which should be used instead.
        df = pd.DataFrame(
            {
                "number": [pd.NA],
                "epoch": [2451545.0],
                "RA": [1.0],
                "DEC": [1.0],
                "band": ["V"],
                "observatory": ["500"],
                "note2": [""],
                "desig": ["2020 AB1"],
            }
        )
        with _mpc_observations(df):
            batch = BatchMPC()
            batch.get_observations(["2020 AB1"])
        assert batch.table["number"].iloc[0] == "2020 AB1"

    def test_no_identifier_found_raises(self):
        # If neither a number nor a designation can be found, the provided
        # code is almost certainly invalid - fail loudly rather than
        # silently dropping the object, so bad input surfaces immediately.
        df = pd.DataFrame(
            {
                "number": [pd.NA],
                "epoch": [2451545.0],
                "RA": [1.0],
                "DEC": [1.0],
                "band": ["V"],
                "observatory": ["500"],
                "note2": [""],
                "desig": [pd.NA],
            }
        )
        with _mpc_observations(df):
            batch = BatchMPC()
            with pytest.raises(ValueError, match="Could not find a valid identifier"):
                batch.get_observations(["bogus"])


# ---------------------------------------------------------------------------
# get_observations: input-validation guard clauses. These raise before any
# network call is made, so no mocking is needed here.
# ---------------------------------------------------------------------------
def test_rejects_non_list_mpccodes():
    with pytest.raises(ValueError):
        BatchMPC().get_observations("433")


def test_rejects_mismatched_id_types_length():
    with pytest.raises(ValueError):
        BatchMPC().get_observations([1, 2], id_types=["asteroid_number"])


def test_rejects_invalid_code_type():
    with pytest.raises(ValueError):
        BatchMPC().get_observations([1.5])


# ---------------------------------------------------------------------------
# read_80_column_data: the 80-column module entry point around parse_80cols_file()
# and the shared augmented-table-to-TrackingData path.
# ---------------------------------------------------------------------------
def test_read_80_column_data_reads_file(tmp_path):
    observation_file = tmp_path / "eros_obs.txt"
    observation_file.write_text(
        "00433         S2021 06 07.42640918 08 15.401-41 22 02.35         12.0 V      500\n"
        "00433         s2021 06 07.4264091 -198301.940 +198171.039 +56287.9850   ~6oMXC57"
    )

    tracking_data, supplementary_data = read_80_column_data(
        [str(observation_file)], custom_name="Eros"
    )

    assert len(tracking_data) == 1
    assert tracking_data[0].time_scale == "UTC"
    assert tracking_data[0].observable_type == "AngularPosition"
    assert len(supplementary_data) == 1
    assert supplementary_data[0].body_name == "500"
    assert supplementary_data[0].reference_point_name == ""


# ---------------------------------------------------------------------------
# _refresh_metadata / __add__ / copy: bookkeeping around the internal table
# populated by MPC retrieval.
# ---------------------------------------------------------------------------
def test_refresh_metadata_computes_epoch_bounds_and_lists(valid_obs_table):
    batch = _batch_from_optical_table(valid_obs_table)
    assert batch.epoch_start < batch.epoch_end
    assert set(batch.observatories) == {"500", "673"}


def test_add_merges_and_sorts_by_epoch(valid_obs_table):
    # __add__ concatenates two batches and re-sorts by epoch, regardless of
    # the order the original rows were in.
    a = _batch_from_optical_table(valid_obs_table.iloc[[1]])
    b = _batch_from_optical_table(valid_obs_table.iloc[[0]])
    combined = a + b
    assert combined.size == 2
    assert combined.table["epoch"].is_monotonic_increasing


def test_copy_is_independent_of_original(valid_obs_table):
    # copy() must deep-copy the underlying table so mutating the clone
    # (including via _refresh_metadata) never leaks back into the original.
    original = _batch_from_optical_table(valid_obs_table)
    clone = original.copy()
    clone._table = clone._table.iloc[0:0]
    clone._refresh_metadata()
    assert original.size == len(valid_obs_table)


def test_mpc_custom_name_metadata(valid_obs_table):
    # Provide a custom name
    custom_name = "Death_Star"
    original = _batch_from_optical_table(valid_obs_table, custom_name=custom_name)
    # Verify the custom_name column exists and has the correct value
    assert "custom_name" in original.table.columns
    assert (original.table["custom_name"] == custom_name).all()

    # Verify MPC_objects returns the custom name instead of mpc code
    assert custom_name in original.MPC_objects
    assert valid_obs_table["number"][0] not in original.MPC_objects
