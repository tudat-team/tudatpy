import os
import json
import pytest
from unittest.mock import patch, MagicMock
from tudatpy.data.spacetrack import SpaceTrackQuery, OMMUtils
from datetime import datetime, timedelta
import numpy as np


# --- Unit Tests (Offline) ---
@pytest.fixture
def temp_tle_folder(tmp_path):
    # Creates a temporary foler path called tle_data
    return str(tmp_path / "tle_data")


# Global patch for login to avoid API authentication
# This will run on all offline tests (autouse = True) and will prevent real API authentication
# It replaces _login() with a dummy static function (static = returns None)
# Every time a query is created through: query = SpaceTrackQuery(username="u", password="p", tle_data_folder=temp_tle_folder),
# its _init_ has _login, but this is patched so the actual _login method will never run.
@pytest.fixture(autouse=True)
def mock_login(request):
    # If the test is marked as remote_data, do NOT patch login
    if request.node.get_closest_marker("remote_data"):
        yield
    else:
        with patch("tudatpy.data.spacetrack.SpaceTrackQuery._login", return_value=None):
            yield


# --- 1. Testing URL/Query Construction ---
def test_descending_epoch_url_construction(temp_tle_folder):
    """Verifies that descending_epoch builds the correct URL and filename."""
    query = SpaceTrackQuery(username="u", password="p", tle_data_folder=temp_tle_folder)

    # _get_json_and_save is patched, because it hits the API. This replaces internal function with a mock.
    with patch.object(query, "_get_json_and_save") as mock_get_json_and_save:
        # Calls the function being tested
        query.descending_epoch(N=10, update_existing=True)

        # Verify URL parts: check if 'limit/10' and 'epoch desc' are in the call
        args, kwargs = mock_get_json_and_save.call_args
        url = args[0]
        filename = args[1]

        # Assert step to verify proper url construction
        assert "limit/10" in url
        assert "orderby/epoch%20desc" in url or "orderby/epoch desc" in url
        assert filename == "gp_descending_limit_10.json"
        assert kwargs["merge"] is True


def test_filtered_by_oe_dict_url(temp_tle_folder):
    """Verifies that orbital element filters are correctly appended to the query."""
    query = SpaceTrackQuery(username="u", password="p", tle_data_folder=temp_tle_folder)

    # _get_json_and_save is patched, because it hits the API. This replaces internal function with a mock.
    with patch.object(query, "_get_json_and_save") as mock_get_json_and_save:
        filters = {"INCLINATION": (97.0, 99.0), "SEMIMAJOR_AXIS": (None, 7000.0)}
        query.filtered_by_oe_dict(filters, limit=50)

        url = mock_get_json_and_save.call_args[0][0]
        # Check for range logic
        assert "INCLINATION/97.0--99.0" in url
        # Check for "less than" logic
        assert "SEMIMAJOR_AXIS/<7000.0" in url
        assert "limit/50" in url


# --- 2. Testing Caching & Cooldown Logic ---
def test_get_tles_date_range_cooldown_logic(temp_tle_folder, capsys):
    """Verifies the 1.5h cooldown logic specifically."""
    query = SpaceTrackQuery(username="u", password="p", tle_data_folder=temp_tle_folder)
    norad_id = "25544"
    filepath = os.path.join(temp_tle_folder, f"tle_{norad_id}.json")

    # Create a mock file that was "hit" 10 minutes ago
    recent_hit = (datetime.now() - timedelta(minutes=10)).isoformat()
    mock_content = {
        "last_api_hit": recent_hit,
        "data": [],  # Empty data to trigger a "needs_fetch" if not for cooldown
    }

    with open(filepath, "w") as f:
        json.dump(mock_content, f)

    # Attempt to fetch
    # we don't need a patch for get_tles_for_date_range, because
    # this function first goes to the local file and checks the latest API hit timestamp.
    query.get_tles_for_date_range(norad_id, "2025-01-01", "2025-01-02")

    # What we want to do here is to test the captured output of a recent API hit.
    # Users will know that was a recent API hit because of the printed out text.
    # We capture that text
    captured = capsys.readouterr()
    assert "Cache is recent" in captured.out

    # Now we force API hit, even if the last hit was recent, by setting override_last_api_hit = True
    # In a real scenario, override_last_api_hit will trigger _fetch_json, which requires network connection.
    # we patch _fetch_json to avoid that, but we assert whether mock_fetch_json gets triggered or not (must be).
    with patch.object(query, "_fetch_json", return_value=[]) as mock_fetch_json:
        query.get_tles_for_date_range(
            norad_id, "2025-01-01", "2025-01-02", override_last_api_hit=True
        )
        assert mock_fetch_json.called  # Fetch should be triggered


# --- 3. Testing Batch Processing ---
def test_save_batch_to_individual_files(temp_tle_folder):
    """Verifies that a batch list is correctly split into per-satellite files."""

    # Create query (this triggers mock_login)
    query = SpaceTrackQuery(username="u", password="p", tle_data_folder=temp_tle_folder)
    # Create mock batch data, very minimal.
    batch_data = [
        {"NORAD_CAT_ID": "111", "EPOCH": "2025-01-01T00:00:00"},
        {"NORAD_CAT_ID": "111", "EPOCH": "2025-01-01T01:00:00"},
        {"NORAD_CAT_ID": "222", "EPOCH": "2025-01-01T00:00:00"},
    ]

    # Call the function to test
    OMMUtils.save_batch_to_individual_files(batch_data, temp_tle_folder)

    # Verify that two separate files are created (gp_111_limit_1 111 will have 2 entries)
    assert os.path.exists(os.path.join(temp_tle_folder, "222_from_batch.json"))
    assert os.path.exists(os.path.join(temp_tle_folder, "111_from_batch.json"))

    # Open Files
    with open(os.path.join(temp_tle_folder, "111_from_batch.json"), "r") as first:
        first_dict = json.load(first)
    with open(os.path.join(temp_tle_folder, "222_from_batch.json"), "r") as second:
        second_dict = json.load(second)
    # Verify that the two files contain the dictionaries in batch data
    assert first_dict == [
        batch_data[0],
        batch_data[1],
    ]  # this also checks that both 111 entries are in the same json
    assert second_dict == [batch_data[2]]


def test_get_tles_skips_omm_records_without_tle_lines():
    records = [
        {
            "NORAD_CAT_ID": "111",
            "TLE_LINE1": "line 1",
            "TLE_LINE2": "line 2",
        },
        {"NORAD_CAT_ID": "222"},
        {"NORAD_CAT_ID": "333", "TLE_LINE1": "line 1"},
        {"NORAD_CAT_ID": "444", "TLE_LINE1": "", "TLE_LINE2": "line 2"},
    ]

    assert OMMUtils.get_tles(records) == {"111": [("line 1", "line 2")]}


# --- 4. Cleaning Logic ---
def test_clean_file_dict_format(temp_tle_folder):
    """Verifies cleaning works on the dictionary/metadata format."""
    os.makedirs(temp_tle_folder, exist_ok=True)
    filepath = os.path.join(temp_tle_folder, "to_clean.json")

    dirty_data = {
        "last_api_hit": "some_time",
        "data": [
            {"NORAD_CAT_ID": "1", "EPOCH": "E1", "CREATION_DATE": "2020"},
            {"NORAD_CAT_ID": "1", "EPOCH": "E1", "CREATION_DATE": "2021"},  # Newer duplicate
        ],
    }
    with open(filepath, "w") as f:
        json.dump(dirty_data, f)

    OMMUtils.clean_file(filepath)

    with open(filepath, "r") as f:
        cleaned = json.load(f)

    assert len(cleaned["data"]) == 1
    assert cleaned["data"][0]["CREATION_DATE"] == "2021"
    assert "last_api_hit" in cleaned  # Metadata preserved


def test_save_unique_sorted_logic(mock_login, temp_tle_folder):
    """Verifies deduplication logic without any network calls."""
    query = SpaceTrackQuery(username="u", password="p", tle_data_folder=temp_tle_folder)
    filepath = os.path.join(temp_tle_folder, "test_logic.json")

    # Simulate an existing file with one object and some entries
    existing_data = [
        {
            "NORAD_CAT_ID": "25544",
            "EPOCH": "2025-01-01T12:00:00",
            "CREATION_DATE": "2025-01-01T00:00:00",
        }
    ]
    query._save_unique_sorted(filepath, existing_data)

    # New data with same ID and Epoch but NEWER creation date
    new_data = [
        {
            "NORAD_CAT_ID": "25544",
            "EPOCH": "2025-01-01T12:00:00",
            "CREATION_DATE": "2025-01-01T10:00:00",
        }
    ]
    query._save_unique_sorted(filepath, new_data)

    with open(filepath, "r") as f:
        content = json.load(f)

    # Should still only have 1 entry, but with the newer creation date
    assert len(content["data"]) == 1
    assert content["data"][0]["CREATION_DATE"] == "2025-01-01T10:00:00"


def test_latest_on_orbit_custom_filename(temp_tle_folder):
    """Verifies that latest_on_orbit respects the custom filename parameter."""
    query = SpaceTrackQuery(username="u", password="p", tle_data_folder=temp_tle_folder)
    custom_name = "my_custom_orbit_file.json"

    # Mock data to return
    mock_data = [{"NORAD_CAT_ID": "1", "EPOCH": "2025-01-01T00:00:00"}]

    # We mock the session.get to return a mock response
    with patch.object(query.session, "get") as mock_get:
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = mock_data
        mock_get.return_value = mock_response

        query.latest_on_orbit(filename=custom_name)

    assert os.path.exists(os.path.join(temp_tle_folder, custom_name))


def test_descending_epoch_update_existing_logic(temp_tle_folder):
    """Verifies that update_existing=True merges data rather than overwriting."""
    query = SpaceTrackQuery(username="u", password="p", tle_data_folder=temp_tle_folder)
    master_file = "master_catalog.json"

    # Setup a reusable mock response helper
    def mock_api_call(data):
        m_get = MagicMock()
        m_get.status_code = 200
        m_get.json.return_value = data
        return m_get

    # 1. First Download
    data_1 = [{"NORAD_CAT_ID": "111", "EPOCH": "2025-01-01T00:00:00"}]
    with patch.object(query.session, "get", return_value=mock_api_call(data_1)):
        query.descending_epoch(N=1, filename=master_file, update_existing=False)

    # 2. Second Download (Merge)
    data_2 = [{"NORAD_CAT_ID": "222", "EPOCH": "2025-01-02T00:00:00"}]
    with patch.object(query.session, "get", return_value=mock_api_call(data_2)):
        query.descending_epoch(N=1, filename=master_file, update_existing=True)

    with open(os.path.join(temp_tle_folder, master_file), "r") as f:
        content = json.load(f)
        data_in_file = content["data"]

    ids = [item["NORAD_CAT_ID"] for item in data_in_file]
    assert "111" in ids and "222" in ids
    assert len(data_in_file) == 2


def test_descending_epoch_overwrite_logic(temp_tle_folder):
    """Verifies that update_existing=False (default) overwrites existing files."""
    query = SpaceTrackQuery(username="u", password="p", tle_data_folder=temp_tle_folder)
    target_file = "overwrite_test.json"

    def mock_api_call(data):
        m_get = MagicMock()
        m_get.status_code = 200
        m_get.json.return_value = data
        return m_get

    # 1. Initial save
    with patch.object(
        query.session, "get", return_value=mock_api_call([{"NORAD_CAT_ID": "1", "EPOCH": "A"}])
    ):
        query.descending_epoch(filename=target_file)

    # 2. Overwrite save
    with patch.object(
        query.session, "get", return_value=mock_api_call([{"NORAD_CAT_ID": "2", "EPOCH": "B"}])
    ):
        query.descending_epoch(filename=target_file, update_existing=False)

    with open(os.path.join(temp_tle_folder, target_file), "r") as f:
        content = json.load(f)
        assert content["data"][0]["NORAD_CAT_ID"] == "2"
        assert len(content["data"]) == 1


def test_get_tudat_keplerian_element_set_full_conversion(mock_login, temp_tle_folder):
    """
    Tests the full extraction and conversion process:
    - km to meters (SMA)
    - degrees to radians (i, Omega, raan, M)
    - Mean to True anomaly
    """
    query = SpaceTrackQuery(username="u", password="p", tle_data_folder=temp_tle_folder)

    # Define inputs in "Space-Track" units (km and deg)
    sma_km = 7000.0
    ecc = 0.1
    inc_deg = 51.6
    raan_deg = 120.0
    argp_deg = 45.0
    mean_anom_deg = 90.0

    mock_omm_list = [
        {
            "SEMIMAJOR_AXIS": str(sma_km),
            "ECCENTRICITY": str(ecc),
            "INCLINATION": str(inc_deg),
            "RA_OF_ASC_NODE": str(raan_deg),
            "ARG_OF_PERICENTER": str(argp_deg),
            "MEAN_ANOMALY": str(mean_anom_deg),
            "NORAD_CAT_ID": "25544",
        }
    ]

    # get_tudat_keplerian_element_set wants a single dictionary, not a list
    with pytest.raises(TypeError):
        OMMUtils.get_tudat_keplerian_element_set(mock_omm_list)

    mock_json_dict = mock_omm_list[0]
    # Execute conversion
    kep = OMMUtils.get_tudat_keplerian_element_set(mock_json_dict)

    # Extract results: (a, e, i, omega, raan, true_anomaly)
    a_res, e_res, i_res, omega_res, raan_res, nu_res = kep

    # 1. Check Semi-major axis (km -> m)
    assert a_res == sma_km * 1000.0

    # 2. Check Eccentricity (dimensionless, should remain same)
    assert e_res == ecc

    # 3. Check Angular elements (deg -> rad)
    assert i_res == pytest.approx(np.deg2rad(inc_deg))
    assert raan_res == pytest.approx(np.deg2rad(raan_deg))
    assert omega_res == pytest.approx(np.deg2rad(argp_deg))

    # 4. We do not check mean-true anomaly conversion,
    # since that is a tudapty function and we know it works well already.


# --- Integration Tests (Online) ---
@pytest.mark.remote_data(required_env=("SPACETRACK_USER", "SPACETRACK_PASS"), service="Space-Track")
def test_live_api_query(temp_tle_folder):
    """Actually hit the Space-Track API."""

    username = os.getenv("SPACETRACK_USER")
    password = os.getenv("SPACETRACK_PASS")

    query = SpaceTrackQuery(username=username, password=password, tle_data_folder=temp_tle_folder)

    # ISS ID: check if we get a valid list back
    res = query.get_tles_for_date_range("25544", "2025-01-01", "2026-01-02")

    assert isinstance(res, list)
    if len(res) > 0:
        assert res[0]["NORAD_CAT_ID"] == "25544"


@pytest.mark.remote_data(required_env=("SPACETRACK_USER", "SPACETRACK_PASS"), service="Space-Track")
def test_custom_query_from_url_live(temp_tle_folder):
    """
    Verifies that a full Space-Track GUI-style URL
    can be ingested and returns valid OMM data.
    """

    username = os.getenv("SPACETRACK_USER")
    password = os.getenv("SPACETRACK_PASS")

    query = SpaceTrackQuery(username=username, password=password, tle_data_folder=temp_tle_folder)

    full_url = (
        "https://www.space-track.org/"
        "basicspacedata/query/class/gp/"
        "NORAD_CAT_ID/25544/"
        "orderby/epoch desc/"
        "limit/1/format/json"
    )

    half_url = "basicspacedata/query/class/gp/NORAD_CAT_ID/25544/limit/1/format/json"

    data_full_url = query.query_from_query_builder_url(
        full_url, output_file="custom_test.json", update_existing=False
    )

    data_half_url = query.query_from_query_builder_url(
        half_url, output_file="custom_test_half_url.json", update_existing=False
    )

    # --- Assertions ---
    assert data_full_url is not None
    assert isinstance(data_full_url, list)
    assert len(data_full_url) == 1

    record = data_full_url[0]

    # Core OMM fields expected from Space-Track GP class
    assert "NORAD_CAT_ID" in record
    assert "EPOCH" in record
    assert "TLE_LINE1" in record
    assert "TLE_LINE2" in record

    # Verify compatibility with OMMUtils
    tles = OMMUtils.get_tles(data_full_url)
    assert "25544" in tles

    assert data_half_url is not None
    assert isinstance(data_half_url, list)
    assert len(data_half_url) == 1

    record = data_half_url[0]

    # Core OMM fields expected from Space-Track GP class
    assert "NORAD_CAT_ID" in record
    assert "EPOCH" in record
    assert "TLE_LINE1" in record
    assert "TLE_LINE2" in record

    # Verify compatibility with OMMUtils
    tles = OMMUtils.get_tles(data_half_url)
    assert "25544" in tles
