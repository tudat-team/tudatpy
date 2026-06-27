import json
import os
import tempfile
import pytest
import numpy as np
from unittest.mock import patch, MagicMock

from tudatpy.data.UDL.batch_utas import BatchUTAS, StationPairObservations, UTASMetadata

# =========================================================================
# Fixtures: create temporary UTAS JSON files for testing
# =========================================================================

SINGLE_STATION_PAIR_OBS = [
    {
        "frequency": 2260760000,
        "tdoa": 0.0049952785,
        "tdoaUnc": -0.00000629,
        "fdoa": 1004.84432812,
        "fdoaUnc": -0.01933109,
        "obTime": "2025-01-02 00:40:00.611999988555908Z",
        "senlat": -14.3754580415,
        "senlon": 132.152376453,
        "senalt": 0.1893,
        "sen2lat": -42.80358333,
        "sen2lon": 147.4405278,
        "sen2alt": 0.0,
        "origSensorId1": "KATHERINE",
        "origSensorId2": "HOBART",
        "ucts": 0,
        "sensor1Delay": 0.0,
        "sensor2Delay": 0.0,
        "bandwidth": 0.0,
        "source": "Unknown",
        "dataMode": "REAL",
        "origin": "Converted on Bingus",
        "satNo": 53365,
    },
    {
        "frequency": 2260760000,
        "tdoa": 0.0049948385,
        "tdoaUnc": 0.0000064,
        "fdoa": 1004.7833303,
        "fdoaUnc": 0.0133036,
        "obTime": "2025-01-02 00:40:01.611999988555908Z",
        "senlat": -14.3754580415,
        "senlon": 132.152376453,
        "senalt": 0.1893,
        "sen2lat": -42.80358333,
        "sen2lon": 147.4405278,
        "sen2alt": 0.0,
        "origSensorId1": "KATHERINE",
        "origSensorId2": "HOBART",
        "ucts": 0,
        "sensor1Delay": 0.0,
        "sensor2Delay": 0.0,
        "bandwidth": 0.0,
        "source": "Unknown",
        "dataMode": "REAL",
        "origin": "Converted on Bingus",
        "satNo": 53365,
    },
]


SECOND_STATION_PAIR_OBS = [
    {
        "frequency": 2260760000,
        "tdoa": 0.001234,
        "tdoaUnc": 0.000001,
        "fdoa": 500.0,
        "fdoaUnc": 0.01,
        "obTime": "2025-01-02 00:40:00.611999988555908Z",
        "senlat": 40.0,
        "senlon": -74.0,
        "senalt": 0.05,
        "sen2lat": 51.0,
        "sen2lon": -0.1,
        "sen2alt": 0.02,
        "origSensorId1": "STATION_A",
        "origSensorId2": "STATION_B",
        "ucts": 0,
        "sensor1Delay": 0.0,
        "sensor2Delay": 0.0,
        "bandwidth": 0.0,
        "source": "Unknown",
        "dataMode": "REAL",
        "origin": "Test",
        "satNo": 53365,
    },
]


def _write_json(file_path, data):
    """Helper to write data as JSON to a file."""
    with open(file_path, "w") as f:
        json.dump(data, f)


@pytest.fixture
def single_pair_file(tmp_path):
    """Create a temporary JSON file with a single station pair."""
    file_path = tmp_path / "single_pair.json"
    _write_json(str(file_path), SINGLE_STATION_PAIR_OBS)
    return str(file_path)


@pytest.fixture
def second_pair_file(tmp_path):
    """Create a temporary JSON file with a second station pair."""
    file_path = tmp_path / "second_pair.json"
    _write_json(str(file_path), SECOND_STATION_PAIR_OBS)
    return str(file_path)


@pytest.fixture
def two_files_single_pair(tmp_path):
    """Create two temporary JSON files, each with the same single station pair."""
    file1 = tmp_path / "file1.json"
    file2 = tmp_path / "file2.json"
    _write_json(str(file1), SINGLE_STATION_PAIR_OBS[:1])
    _write_json(str(file2), SINGLE_STATION_PAIR_OBS[1:])
    return [str(file1), str(file2)]


@pytest.fixture
def empty_file(tmp_path):
    """Create a temporary empty JSON file (no observations)."""
    file_path = tmp_path / "empty.json"
    _write_json(str(file_path), [])
    return str(file_path)


@pytest.fixture
def invalid_json_file(tmp_path):
    """Create a temporary file with invalid JSON content."""
    file_path = tmp_path / "invalid.json"
    with open(file_path, "w") as f:
        f.write("{invalid json content")
    return str(file_path)


@pytest.fixture
def wrong_structure_file(tmp_path):
    """Create a temporary JSON file with wrong structure (not a list, no 'observations' key)."""
    file_path = tmp_path / "wrong_structure.json"
    _write_json(str(file_path), {"key": "value"})
    return str(file_path)


# =========================================================================
# Mock-based tests: BatchUTAS with mocked Tudat dependencies
# =========================================================================


class TestBatchUTASProperties:
    """Test basic properties of BatchUTAS with mocked environment."""

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_target_id(self, mock_time_rep, mock_env_setup, single_pair_file):
        """Test that target_id is correctly extracted from data."""
        batch = BatchUTAS([single_pair_file])
        assert batch.target_id == "53365"

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_num_observations_single_file(self, mock_time_rep, mock_env_setup, single_pair_file):
        """Test num_observations with a single file."""
        batch = BatchUTAS([single_pair_file])
        assert batch.num_observations == 2

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_num_observations_two_files(self, mock_time_rep, mock_env_setup, two_files_single_pair):
        """Test num_observations sums across multiple files."""
        batch = BatchUTAS(two_files_single_pair)
        assert batch.num_observations == 2

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_station_pairs(self, mock_time_rep, mock_env_setup, single_pair_file):
        """Test that station pairs are correctly identified."""
        batch = BatchUTAS([single_pair_file])
        assert batch.station_pairs == [("KATHERINE", "HOBART")]

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_station_names(self, mock_time_rep, mock_env_setup, single_pair_file):
        """Test that station names are correctly extracted."""
        batch = BatchUTAS([single_pair_file])
        assert batch.station_names == {"KATHERINE", "HOBART"}

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_num_station_pairs(self, mock_time_rep, mock_env_setup, single_pair_file):
        """Test num_station_pairs for single pair."""
        batch = BatchUTAS([single_pair_file])
        assert batch.num_station_pairs == 1

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_num_station_pairs_multi(
        self, mock_time_rep, mock_env_setup, single_pair_file, second_pair_file
    ):
        """Test num_station_pairs for multiple pairs across files."""
        batch = BatchUTAS([single_pair_file, second_pair_file])
        assert batch.num_station_pairs == 2

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_station_names_multi(
        self, mock_time_rep, mock_env_setup, single_pair_file, second_pair_file
    ):
        """Test station_names with multiple station pairs across files."""
        batch = BatchUTAS([single_pair_file, second_pair_file])
        assert batch.station_names == {"KATHERINE", "HOBART", "STATION_A", "STATION_B"}


class TestBatchUTASMetadata:
    """Test metadata extraction."""

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_metadata_fields(self, mock_time_rep, mock_env_setup, single_pair_file):
        """Test that all metadata fields are correctly populated."""
        batch = BatchUTAS([single_pair_file])
        metadata = batch.get_metadata()

        assert isinstance(metadata, UTASMetadata)
        assert metadata.target_id == "53365"
        assert metadata.frequency == 2260760000
        assert metadata.bandwidth == 0.0
        assert metadata.sensor1_delay == 0.0
        assert metadata.sensor2_delay == 0.0
        assert metadata.data_mode == "REAL"
        assert metadata.origin == "Converted on Bingus"
        assert metadata.source == "Unknown"
        assert metadata.ucts == 0


class TestBatchUTASGetObservations:
    """Test get_observations_for_station_pair method."""

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_get_observations_success(self, mock_time_rep, mock_env_setup, single_pair_file):
        """Test getting observations for a valid station pair."""
        batch = BatchUTAS([single_pair_file])
        obs = batch.get_observations_for_station_pair(("KATHERINE", "HOBART"))

        assert isinstance(obs, StationPairObservations)
        assert len(obs) == 2
        assert len(obs.tdoa) == 2
        assert len(obs.fdoa) == 2
        assert len(obs.epochs) == 2

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_get_observations_invalid_pair(self, mock_time_rep, mock_env_setup, single_pair_file):
        """Test getting observations for a non-existent station pair raises error."""
        batch = BatchUTAS([single_pair_file])
        with pytest.raises(RuntimeError, match="Station pair .* not found"):
            batch.get_observations_for_station_pair(("NONEXISTENT", "ALSO_NONEXISTENT"))

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_get_all_observations(self, mock_time_rep, mock_env_setup, single_pair_file):
        """Test get_all_observations_by_station_pair returns correct dict."""
        batch = BatchUTAS([single_pair_file])
        all_obs = batch.get_all_observations_by_station_pair()

        assert isinstance(all_obs, dict)
        assert len(all_obs) == 1
        assert ("KATHERINE", "HOBART") in all_obs
        assert len(all_obs[("KATHERINE", "HOBART")]) == 2


class TestBatchUTASMultipleFiles:
    """Test BatchUTAS with multiple input files."""

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_two_files_same_station_pair(
        self, mock_time_rep, mock_env_setup, two_files_single_pair
    ):
        """Test loading two files with the same station pair."""
        batch = BatchUTAS(two_files_single_pair)
        assert batch.num_observations == 2
        assert batch.num_station_pairs == 1
        assert len(batch.get_observations_for_station_pair(("KATHERINE", "HOBART"))) == 2

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_multi_file_multi_pair(
        self, mock_time_rep, mock_env_setup, single_pair_file, second_pair_file
    ):
        """Test loading two files with different station pairs."""
        batch = BatchUTAS([single_pair_file, second_pair_file])
        assert batch.num_observations == 3
        assert batch.num_station_pairs == 2
        assert ("KATHERINE", "HOBART") in batch.station_pairs
        assert ("STATION_A", "STATION_B") in batch.station_pairs


class TestBatchUTASErrorHandling:
    """Test error handling in BatchUTAS."""

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_empty_file_raises_error(self, mock_time_rep, mock_env_setup, empty_file):
        """Test that an empty observations list raises RuntimeError."""
        with pytest.raises(RuntimeError, match="No observations found"):
            BatchUTAS([empty_file])

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_invalid_json_raises_error(self, mock_time_rep, mock_env_setup, invalid_json_file):
        """Test that invalid JSON raises RuntimeError."""
        with pytest.raises(RuntimeError, match="JSON parse error"):
            BatchUTAS([invalid_json_file])

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_wrong_structure_raises_error(
        self, mock_time_rep, mock_env_setup, wrong_structure_file
    ):
        """Test that wrong JSON structure raises RuntimeError."""
        with pytest.raises(RuntimeError, match="Unexpected JSON structure"):
            BatchUTAS([wrong_structure_file])

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_missing_required_field(self, mock_time_rep, mock_env_setup, tmp_path):
        """Test that missing required field raises RuntimeError."""
        bad_data = [
            {
                "frequency": 2260760000,
                # missing tdoa, fdoa, obTime, etc.
                "origSensorId1": "A",
                "origSensorId2": "B",
                "satNo": 12345,
            }
        ]
        file_path = tmp_path / "bad.json"
        _write_json(str(file_path), bad_data)
        with pytest.raises(RuntimeError):
            BatchUTAS([str(file_path)])

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_multiple_targets_raises_error(self, mock_time_rep, mock_env_setup, tmp_path):
        """Test that multiple targets in one file raises RuntimeError."""
        multi_target_data = [
            {
                "frequency": 2260760000,
                "tdoa": 0.004,
                "tdoaUnc": 0.0,
                "fdoa": 100.0,
                "fdoaUnc": 0.0,
                "obTime": "2025-01-02 00:40:00.611999988555908Z",
                "senlat": -14.375,
                "senlon": 132.152,
                "senalt": 0.1893,
                "sen2lat": -42.803,
                "sen2lon": 147.440,
                "sen2alt": 0.0,
                "origSensorId1": "A",
                "origSensorId2": "B",
                "ucts": 0,
                "sensor1Delay": 0.0,
                "sensor2Delay": 0.0,
                "bandwidth": 0.0,
                "source": "Unknown",
                "dataMode": "REAL",
                "origin": "Test",
                "satNo": 11111,
            },
            {
                "frequency": 2260760000,
                "tdoa": 0.005,
                "tdoaUnc": 0.0,
                "fdoa": 200.0,
                "fdoaUnc": 0.0,
                "obTime": "2025-01-02 00:40:01.611999988555908Z",
                "senlat": -14.375,
                "senlon": 132.152,
                "senalt": 0.1893,
                "sen2lat": -42.803,
                "sen2lon": 147.440,
                "sen2alt": 0.0,
                "origSensorId1": "A",
                "origSensorId2": "B",
                "ucts": 0,
                "sensor1Delay": 0.0,
                "sensor2Delay": 0.0,
                "bandwidth": 0.0,
                "source": "Unknown",
                "dataMode": "REAL",
                "origin": "Test",
                "satNo": 22222,
            },
        ]
        file_path = tmp_path / "multi_target.json"
        _write_json(str(file_path), multi_target_data)
        with pytest.raises(RuntimeError, match="Multiple targets detected"):
            BatchUTAS([str(file_path)])

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_multiple_targets_across_files(self, mock_time_rep, mock_env_setup, tmp_path):
        """Test that multiple targets across files raises RuntimeError."""
        data1 = [
            {
                "frequency": 2260760000,
                "tdoa": 0.004,
                "tdoaUnc": 0.0,
                "fdoa": 100.0,
                "fdoaUnc": 0.0,
                "obTime": "2025-01-02 00:40:00.611999988555908Z",
                "senlat": -14.375,
                "senlon": 132.152,
                "senalt": 0.1893,
                "sen2lat": -42.803,
                "sen2lon": 147.440,
                "sen2alt": 0.0,
                "origSensorId1": "A",
                "origSensorId2": "B",
                "ucts": 0,
                "sensor1Delay": 0.0,
                "sensor2Delay": 0.0,
                "bandwidth": 0.0,
                "source": "Unknown",
                "dataMode": "REAL",
                "origin": "Test",
                "satNo": 11111,
            }
        ]
        data2 = [
            {
                "frequency": 2260760000,
                "tdoa": 0.005,
                "tdoaUnc": 0.0,
                "fdoa": 200.0,
                "fdoaUnc": 0.0,
                "obTime": "2025-01-02 00:40:01.611999988555908Z",
                "senlat": -14.375,
                "senlon": 132.152,
                "senalt": 0.1893,
                "sen2lat": -42.803,
                "sen2lon": 147.440,
                "sen2alt": 0.0,
                "origSensorId1": "A",
                "origSensorId2": "B",
                "ucts": 0,
                "sensor1Delay": 0.0,
                "sensor2Delay": 0.0,
                "bandwidth": 0.0,
                "source": "Unknown",
                "dataMode": "REAL",
                "origin": "Test",
                "satNo": 22222,
            }
        ]
        file1 = tmp_path / "target1.json"
        file2 = tmp_path / "target2.json"
        _write_json(str(file1), data1)
        _write_json(str(file2), data2)
        with pytest.raises(RuntimeError, match="Multiple targets detected"):
            BatchUTAS([str(file1), str(file2)])


class TestBatchUTASToTudat:
    """Test to_tudat method with mocked Tudat objects."""

    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    def test_to_tudat_returns_observation_collection(
        self, mock_env_setup, mock_time_rep, single_pair_file
    ):
        """Test that to_tudat returns an ObservationCollection."""
        # Mock the environment setup
        mock_body = MagicMock()
        mock_body.get_shape_model.return_value = None
        mock_body.ground_station_list = []

        mock_bodies = MagicMock()
        mock_bodies.get_body.return_value = mock_body

        mock_env_setup.create_body_shape_model.return_value = MagicMock()
        mock_env_setup.from_spice_oblate_spherical_body_shape_settings.return_value = {}
        mock_env_setup.ground_station.basic_station.return_value = MagicMock()
        mock_env_setup.add_ground_station = MagicMock()

        # Mock time conversion
        mock_time_scale_converter = MagicMock()
        mock_time_scale_converter.convert_time.return_value = 1000.0
        mock_time_rep.default_time_scale_converter.return_value = mock_time_scale_converter
        mock_time_rep.utc_scale = "UTC"
        mock_time_rep.tdb_scale = "TDB"

        batch = BatchUTAS([single_pair_file])

        # Mock the observation collection creation
        with patch(
            "tudatpy.data.UDL.batch_utas.observations.ObservationCollection"
        ) as mock_obs_collection, patch(
            "tudatpy.data.UDL.batch_utas.observations.create_single_observation_set"
        ) as mock_create_set, patch(
            "tudatpy.data.UDL.batch_utas.links"
        ) as mock_links, patch(
            "tudatpy.data.UDL.batch_utas.model_settings"
        ) as mock_model_settings:
            mock_links.receiver = "receiver"
            mock_links.receiver2 = "receiver2"
            mock_links.transmitter = "transmitter"
            mock_links.body_reference_point_link_end_id.return_value = "link_end_id"
            mock_links.body_origin_link_end_id.return_value = "target_link_end"
            mock_links.link_definition.return_value = MagicMock()

            mock_model_settings.differenced_time_of_arrival_type = "TDOA"
            mock_model_settings.differenced_frequency_of_arrival_type = "FDOA"

            mock_create_set.return_value = MagicMock()
            mock_obs_collection.return_value = MagicMock()

            result = batch.to_tudat(mock_bodies)

            # Should create 2 observation sets (1 TDOA + 1 FDOA) for 1 station pair
            assert mock_create_set.call_count == 2
            mock_obs_collection.assert_called_once()

    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    def test_to_tudat_multi_pair_creates_correct_number_of_sets(
        self, mock_env_setup, mock_time_rep, single_pair_file, second_pair_file
    ):
        """Test that to_tudat creates correct number of observation sets for multiple pairs."""
        mock_body = MagicMock()
        mock_body.get_shape_model.return_value = None
        mock_body.ground_station_list = []

        mock_bodies = MagicMock()
        mock_bodies.get_body.return_value = mock_body

        mock_env_setup.create_body_shape_model.return_value = MagicMock()
        mock_env_setup.from_spice_oblate_spherical_body_shape_settings.return_value = {}
        mock_env_setup.ground_station.basic_station.return_value = MagicMock()
        mock_env_setup.add_ground_station = MagicMock()

        mock_time_scale_converter = MagicMock()
        mock_time_scale_converter.convert_time.return_value = 1000.0
        mock_time_rep.default_time_scale_converter.return_value = mock_time_scale_converter
        mock_time_rep.utc_scale = "UTC"
        mock_time_rep.tdb_scale = "TDB"

        batch = BatchUTAS([single_pair_file, second_pair_file])
        assert batch.num_station_pairs == 2

        with patch(
            "tudatpy.data.UDL.batch_utas.observations.ObservationCollection"
        ) as mock_obs_collection, patch(
            "tudatpy.data.UDL.batch_utas.observations.create_single_observation_set"
        ) as mock_create_set, patch(
            "tudatpy.data.UDL.batch_utas.links"
        ) as mock_links, patch(
            "tudatpy.data.UDL.batch_utas.model_settings"
        ) as mock_model_settings:
            mock_links.receiver = "receiver"
            mock_links.receiver2 = "receiver2"
            mock_links.transmitter = "transmitter"
            mock_links.body_reference_point_link_end_id.return_value = "link_end_id"
            mock_links.body_origin_link_end_id.return_value = "target_link_end"
            mock_links.link_definition.return_value = MagicMock()

            mock_model_settings.differenced_time_of_arrival_type = "TDOA"
            mock_model_settings.differenced_frequency_of_arrival_type = "FDOA"

            mock_create_set.return_value = MagicMock()
            mock_obs_collection.return_value = MagicMock()

            batch.to_tudat(mock_bodies)

            # 2 station pairs × 2 types (TDOA + FDOA) = 4 observation sets
            assert mock_create_set.call_count == 4


class TestBatchUTASLinkDefinitions:
    """Test link definition generation."""

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_get_link_definitions_default_target(
        self, mock_time_rep, mock_env_setup, single_pair_file
    ):
        """Test link definitions use target_id by default."""
        batch = BatchUTAS([single_pair_file])

        with patch("tudatpy.data.UDL.batch_utas.links") as mock_links:
            mock_links.receiver = "receiver"
            mock_links.receiver2 = "receiver2"
            mock_links.transmitter = "transmitter"
            mock_links.body_reference_point_link_end_id.return_value = "link_end_id"
            mock_links.body_origin_link_end_id.return_value = "target_link_end"
            mock_links.link_definition.return_value = MagicMock()

            link_defs = batch.get_link_definitions()

            assert len(link_defs) == 1
            # Verify body_origin_link_end_id was called with target_id
            mock_links.body_origin_link_end_id.assert_called_with("53365")

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_get_link_definitions_custom_target(
        self, mock_time_rep, mock_env_setup, single_pair_file
    ):
        """Test link definitions use custom target name when provided."""
        batch = BatchUTAS([single_pair_file])

        with patch("tudatpy.data.UDL.batch_utas.links") as mock_links:
            mock_links.receiver = "receiver"
            mock_links.receiver2 = "receiver2"
            mock_links.transmitter = "transmitter"
            mock_links.body_reference_point_link_end_id.return_value = "link_end_id"
            mock_links.body_origin_link_end_id.return_value = "target_link_end"
            mock_links.link_definition.return_value = MagicMock()

            link_defs = batch.get_link_definitions(target_name_override="MySatellite")

            assert len(link_defs) == 1
            mock_links.body_origin_link_end_id.assert_called_with("MySatellite")

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_get_link_definitions_multiple_pairs(
        self, mock_time_rep, mock_env_setup, single_pair_file, second_pair_file
    ):
        """Test link definitions for multiple station pairs."""
        batch = BatchUTAS([single_pair_file, second_pair_file])

        with patch("tudatpy.data.UDL.batch_utas.links") as mock_links:
            mock_links.receiver = "receiver"
            mock_links.receiver2 = "receiver2"
            mock_links.transmitter = "transmitter"
            mock_links.body_reference_point_link_end_id.return_value = "link_end_id"
            mock_links.body_origin_link_end_id.return_value = "target_link_end"
            mock_links.link_definition.return_value = MagicMock()

            link_defs = batch.get_link_definitions()

            assert len(link_defs) == 2


class TestBatchUTASGeodeticConversion:
    """Test geodetic position conversion."""

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_convert_to_tudat_geodetic(self, mock_time_rep, mock_env_setup, single_pair_file):
        """Test geodetic conversion from degrees/meters to radians/meters."""
        batch = BatchUTAS([single_pair_file])

        position = {"altitude": 100.0, "latitude": 45.0, "longitude": -90.0}
        result = batch._convert_to_tudat_geodetic(position)

        assert result[0] == pytest.approx(100.0)  # altitude unchanged (meters)
        assert result[1] == pytest.approx(45.0 * np.pi / 180.0)  # latitude in radians
        assert result[2] == pytest.approx(-90.0 * np.pi / 180.0)  # longitude in radians


class TestBatchUTASStationPairObservations:
    """Test StationPairObservations class."""

    def test_empty_station_pair_observations(self):
        """Test empty StationPairObservations."""
        obs = StationPairObservations()
        assert len(obs) == 0
        assert obs.epochs == []
        assert obs.tdoa == []
        assert obs.fdoa == []

    def test_station_pair_observations_after_append(self):
        """Test StationPairObservations after appending data."""
        obs = StationPairObservations()
        obs.epochs.append(1000.0)
        obs.tdoa.append(0.5)
        obs.tdoa_uncertainties.append(0.01)
        obs.fdoa.append(100.0)
        obs.fdoa_uncertainties.append(0.1)

        assert len(obs) == 1
        assert obs.epochs == [1000.0]
        assert obs.tdoa == [0.5]


class TestBatchUTASUTASMetadata:
    """Test UTASMetadata class."""

    def test_default_metadata(self):
        """Test default UTASMetadata values."""
        meta = UTASMetadata()
        assert meta.target_id == ""
        assert meta.frequency == 0.0
        assert meta.bandwidth == 0.0
        assert meta.sensor1_delay == 0.0
        assert meta.sensor2_delay == 0.0
        assert meta.data_mode == ""
        assert meta.origin == ""
        assert meta.source == ""
        assert meta.ucts == 0


class TestBatchUTASJSONStructures:
    """Test handling of different JSON structures."""

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_dict_with_observations_key(self, mock_time_rep, mock_env_setup, tmp_path):
        """Test loading from dict with 'observations' key."""
        file_path = tmp_path / "dict_structure.json"
        _write_json(str(file_path), {"observations": SINGLE_STATION_PAIR_OBS})
        batch = BatchUTAS([str(file_path)])
        assert batch.num_observations == 2

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_list_structure(self, mock_time_rep, mock_env_setup, single_pair_file):
        """Test loading from list structure (standard)."""
        batch = BatchUTAS([single_pair_file])
        assert batch.num_observations == 2


class TestBatchUTASSatNoAsString:
    """Test satNo handling when provided as string."""

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_satno_as_string(self, mock_time_rep, mock_env_setup, tmp_path):
        """Test that satNo as string is handled correctly."""
        data = [
            {
                "frequency": 2260760000,
                "tdoa": 0.004,
                "tdoaUnc": 0.0,
                "fdoa": 100.0,
                "fdoaUnc": 0.0,
                "obTime": "2025-01-02 00:40:00.611999988555908Z",
                "senlat": -14.375,
                "senlon": 132.152,
                "senalt": 0.1893,
                "sen2lat": -42.803,
                "sen2lon": 147.440,
                "sen2alt": 0.0,
                "origSensorId1": "A",
                "origSensorId2": "B",
                "ucts": 0,
                "sensor1Delay": 0.0,
                "sensor2Delay": 0.0,
                "bandwidth": 0.0,
                "source": "Unknown",
                "dataMode": "REAL",
                "origin": "Test",
                "satNo": "53365",
            }
        ]
        file_path = tmp_path / "string_satno.json"
        _write_json(str(file_path), data)
        batch = BatchUTAS([str(file_path)])
        assert batch.target_id == "53365"


class TestBatchUTASOptionalFields:
    """Test handling of optional fields with defaults."""

    @patch("tudatpy.data.UDL.batch_utas.environment_setup")
    @patch("tudatpy.data.UDL.batch_utas.time_representation")
    def test_optional_fields_use_defaults(self, mock_time_rep, mock_env_setup, tmp_path):
        """Test that missing optional fields use default values."""
        data = [
            {
                "frequency": 2260760000,
                "tdoa": 0.004,
                "fdoa": 100.0,
                "obTime": "2025-01-02 00:40:00.611999988555908Z",
                "senlat": -14.375,
                "senlon": 132.152,
                "senalt": 0.1893,
                "sen2lat": -42.803,
                "sen2lon": 147.440,
                "sen2alt": 0.0,
                "origSensorId1": "A",
                "origSensorId2": "B",
                "ucts": 0,
                "satNo": 53365,
            }
        ]
        file_path = tmp_path / "no_optional.json"
        _write_json(str(file_path), data)
        batch = BatchUTAS([str(file_path)])
        metadata = batch.get_metadata()

        assert metadata.bandwidth == 0.0
        assert metadata.sensor1_delay == 0.0
        assert metadata.sensor2_delay == 0.0
        assert metadata.data_mode == ""
        assert metadata.origin == ""
        assert metadata.source == ""
