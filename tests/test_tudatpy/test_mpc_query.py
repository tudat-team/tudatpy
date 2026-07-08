import pandas as pd
import pytest
from tudatpy.data.mpc import BatchMPC
import numpy as np
from astropy.table import Table


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
# BatchMPC.__init__ calls MPC.get_observatory_codes() over the network to
# build its list of observatories/space telescopes. Every test constructs a
# BatchMPC(), so this fixture is autouse: it stubs that call out with a small
# fake observatory table (C51 has no Longitude, mimicking a space telescope)
# so the whole suite runs offline and deterministically.
@pytest.fixture(autouse=True)
def mock_station_info(monkeypatch):
    """Prevent every BatchMPC() from hitting the network in __init__."""
    fake_obs_codes = pd.DataFrame(
        {
            "Code": ["500", "673", "089", "C51"],
            "Longitude": [0.0, 1.23, 4.56, float("nan")],  # C51 = space telescope
        }
    )

    def fake_get_observatory_codes():
        class _Wrap:
            def to_pandas(self_inner):
                return fake_obs_codes

        return _Wrap()

    monkeypatch.setattr("astroquery.mpc.MPC.get_observatory_codes", fake_get_observatory_codes)
    yield


# A minimal table that satisfies BatchMPC._req_cols, used as the baseline
# input for from_pandas/from_astropy-driven tests below. RA/DEC are in
# degrees (the default `in_degrees=True`), epoch is Julian Day.
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


# ---------------------------------------------------------------------------
# _standardize_dataframe: normalizes observatory codes and MPC numbers
# before they're stored, so downstream lookups/joins are consistent.
# ---------------------------------------------------------------------------
class TestStandardizeDataframe:
    def test_zero_pads_short_observatory_codes(self):
        # MPC observatory codes are always 3 chars (e.g. "89" -> "089");
        # alphanumeric codes like "C51" are already 3 chars and pass through.
        batch = BatchMPC()
        df = pd.DataFrame({"observatory": ["89", "673", "C51"], "number": [1, 2, 3]})
        out = batch._standardize_dataframe(df)
        assert out["observatory"].tolist() == ["089", "673", "C51"]

    def test_strips_whitespace_and_stringifies_number(self):
        batch = BatchMPC()
        df = pd.DataFrame({"number": [" 433 ", 1]})
        out = batch._standardize_dataframe(df)
        assert out["number"].tolist() == ["433", "1"]

    def test_does_not_mutate_input(self):
        # _standardize_dataframe must work on a copy; callers (e.g.
        # get_observations) rely on the original frame staying untouched.
        batch = BatchMPC()
        df = pd.DataFrame({"observatory": ["89"]})
        batch._standardize_dataframe(df)
        assert df["observatory"].iloc[0] == "89"  # original untouched


# ---------------------------------------------------------------------------
# _validate_table: the guard rail shared by from_pandas/from_astropy that
# checks reference frame, required columns, and non-emptiness before any
# data is ingested.
# ---------------------------------------------------------------------------
class TestValidateTable:
    def test_rejects_non_j2000_frame(self, valid_obs_table):
        # Only J2000 is currently supported (see class docstring notes).
        batch = BatchMPC()
        with pytest.raises(NotImplementedError):
            batch._validate_table(valid_obs_table, frame="FK5")

    def test_rejects_missing_required_columns(self):
        batch = BatchMPC()
        df = pd.DataFrame({"number": [1], "epoch": [2451545.0]})  # missing RA/DEC/band/observatory
        with pytest.raises(ValueError, match="mandatory columns"):
            batch._validate_table(df, frame="J2000")

    def test_rejects_empty_table(self, valid_obs_table):
        batch = BatchMPC()
        with pytest.raises(ValueError, match="zero rows"):
            batch._validate_table(valid_obs_table.iloc[0:0], frame="J2000")

    def test_rejects_unsupported_type(self):
        # _validate_table only understands (Q)Table and DataFrame inputs.
        batch = BatchMPC()
        with pytest.raises(TypeError):
            batch._validate_table([1, 2, 3], frame="J2000")


# ---------------------------------------------------------------------------
# RA/DEC unit handling in from_pandas/from_astropy: degrees must be
# converted to radians and wrapped into [-pi, pi]; radians must be passed
# through unchanged when in_degrees=False.
# ---------------------------------------------------------------------------
class TestAngleConversion:
    def test_degrees_converted_and_wrapped_to_minus_pi_pi(self, valid_obs_table):
        batch = BatchMPC()
        batch.from_pandas(valid_obs_table, in_degrees=True)
        ra = batch.table["RA"].to_numpy()
        assert (ra >= -np.pi).all() and (ra <= np.pi).all()
        # 200 deg -> should wrap to a negative radian value, not stay > pi
        assert ra[1] < 0

    def test_radians_left_untouched_when_in_degrees_false(self, valid_obs_table):
        table_rad = valid_obs_table.copy()
        table_rad["RA"] = np.radians(table_rad["RA"])
        table_rad["DEC"] = np.radians(table_rad["DEC"])
        batch = BatchMPC()
        batch.from_pandas(table_rad, in_degrees=False)
        assert np.isclose(batch.table["RA"].iloc[0], table_rad["RA"].iloc[0])


# ---------------------------------------------------------------------------
# _add_time_columns: converts the JD/UTC 'epoch' column into
# epoch_seconds_UTC and epoch_seconds_TDB using tudatpy's own time-scale
# converter. This is deliberately NOT mocked: the conversion is pure, local
# and deterministic, so exercising the real tudatpy code gives a genuine
# regression guard on that contract.
# ---------------------------------------------------------------------------
def test_tdb_is_ahead_of_utc_by_expected_offset(valid_obs_table):
    batch = BatchMPC()
    out = batch._add_time_columns(valid_obs_table)
    diff = out["epoch_seconds_TDB"] - out["epoch_seconds_UTC"]
    # TDB leads UTC by (TAI-UTC leap seconds) + 32.184s (TT-TAI); in the
    # current era that's roughly 69s, so bound it loosely to avoid a brittle
    # test if leap seconds change.
    assert (diff > 60).all() and (diff < 70).all()


# Fakes astroquery's MPC.get_observations(code) return value: the real call
# returns an astropy-Table-like object exposing .to_pandas().
def _fake_mpc_obs(df):
    class _Wrap:
        def to_pandas(self_inner):
            return df

    return _Wrap()


# ---------------------------------------------------------------------------
# get_observations: identifier-resolution logic. This is the trickiest part
# of the class - for each MPC code it must decide the object's canonical
# "number" from whatever combination of number/desig/comet_type the MPC
# response contains. Each test below pins one specific branch of that logic.
# ---------------------------------------------------------------------------
class TestGetObservationsIdentifierResolution:
    def test_numbered_asteroid_uses_number(self, monkeypatch):
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
        monkeypatch.setattr("astroquery.mpc.MPC.get_observations", lambda code: _fake_mpc_obs(df))
        batch = BatchMPC()
        batch.get_observations([433])
        assert batch.table["number"].iloc[0] == "433"

    def test_comet_combines_number_and_comet_type(self, monkeypatch):
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
        monkeypatch.setattr("astroquery.mpc.MPC.get_observations", lambda code: _fake_mpc_obs(df))
        batch = BatchMPC()
        batch.get_observations(["3I"])
        assert batch.table["number"].iloc[0] == "3I"

    def test_no_number_falls_back_to_designation(self, monkeypatch):
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
        monkeypatch.setattr("astroquery.mpc.MPC.get_observations", lambda code: _fake_mpc_obs(df))
        batch = BatchMPC()
        batch.get_observations(["2020 AB1"])
        assert batch.table["number"].iloc[0] == "2020 AB1"

    def test_no_identifier_found_raises(self, monkeypatch):
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
        monkeypatch.setattr("astroquery.mpc.MPC.get_observations", lambda code: _fake_mpc_obs(df))
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
# from_file: a thin wrapper around parse_80cols_file() + from_astropy(). We
# stub the parser (its own correctness is tudatpy's parser test suite's job)
# and just check the wiring: the parsed table flows through to the batch.
#
# Note: patch the name where it's *used* (tudatpy.data.mpc.mpc), not where
# it's *defined* (tudatpy.data.mpc.parser_80col) - mpc.py imports the
# function by name at module load time, so patching the original module
# doesn't affect the copy already bound in mpc.py's namespace.
# ---------------------------------------------------------------------------
def test_from_file_delegates_to_parser_and_from_astropy(monkeypatch, tmp_path):
    fake_table = Table(
        {
            "number": ["433"],
            "epoch": [2451545.0],
            "RA": [0.1],
            "DEC": [0.1],
            "band": ["V"],
            "observatory": ["500"],
        }
    )
    monkeypatch.setattr("tudatpy.data.mpc.mpc.parse_80cols_file", lambda filename: fake_table)
    batch = BatchMPC()
    batch.from_file(str(tmp_path / "fake.txt"), custom_name="Eros")
    assert batch.size == 1
    assert batch.MPC_objects == ["Eros"]


# ---------------------------------------------------------------------------
# _refresh_metadata / __add__ / copy: bookkeeping around the internal table
# that every ingestion path (from_pandas, from_astropy, get_observations,
# from_file) relies on to keep public properties (size, epoch_start/end,
# observatories, ...) in sync.
# ---------------------------------------------------------------------------
def test_refresh_metadata_computes_epoch_bounds_and_lists(valid_obs_table):
    batch = BatchMPC()
    batch.from_pandas(valid_obs_table)
    assert batch.epoch_start < batch.epoch_end
    assert set(batch.observatories) == {"500", "673"}


def test_add_merges_and_sorts_by_epoch(valid_obs_table):
    # __add__ concatenates two batches and re-sorts by epoch, regardless of
    # the order the original rows were in.
    a = BatchMPC()
    a.from_pandas(valid_obs_table.iloc[[1]])
    b = BatchMPC()
    b.from_pandas(valid_obs_table.iloc[[0]])
    combined = a + b
    assert combined.size == 2
    assert combined.table["epoch"].is_monotonic_increasing


def test_copy_is_independent_of_original(valid_obs_table):
    # copy() must deep-copy the underlying table so mutating the clone
    # (including via _refresh_metadata) never leaks back into the original.
    original = BatchMPC()
    original.from_pandas(valid_obs_table)
    clone = original.copy()
    clone._table = clone._table.iloc[0:0]
    clone._refresh_metadata()
    assert original.size == len(valid_obs_table)


def test_mpc_custom_name_metadata(valid_obs_table):
    # Provide a custom name
    custom_name = "Death_Star"
    original = BatchMPC()
    original.from_pandas(valid_obs_table, custom_name=custom_name)
    # Verify the custom_name column exists and has the correct value
    assert "custom_name" in original.table.columns
    assert (original.table["custom_name"] == custom_name).all()

    # Verify MPC_objects returns the custom name instead of mpc code
    assert custom_name in original.MPC_objects
    assert valid_obs_table["number"][0] not in original.MPC_objects
