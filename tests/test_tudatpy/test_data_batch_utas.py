import json
import pytest
import numpy as np

from tudatpy.data.UDL.batch_utas import BatchUTAS

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
    with open(file_path, "w") as f:
        json.dump(data, f)


@pytest.fixture
def single_pair_file(tmp_path):
    file_path = tmp_path / "single_pair.json"
    _write_json(str(file_path), SINGLE_STATION_PAIR_OBS)
    return str(file_path)


@pytest.fixture
def second_pair_file(tmp_path):
    file_path = tmp_path / "second_pair.json"
    _write_json(str(file_path), SECOND_STATION_PAIR_OBS)
    return str(file_path)


@pytest.fixture
def two_files_single_pair(tmp_path):
    file1 = tmp_path / "file1.json"
    file2 = tmp_path / "file2.json"
    _write_json(str(file1), SINGLE_STATION_PAIR_OBS[:1])
    _write_json(str(file2), SINGLE_STATION_PAIR_OBS[1:])
    return [str(file1), str(file2)]


# =========================================================================
# Parsing — core data-loading logic
# =========================================================================


class TestBatchUTASProperties:
    """Verify that parsed properties are correct."""

    def test_target_id(self, single_pair_file):
        batch = BatchUTAS([single_pair_file])
        assert batch.target_id == "53365"

    def test_num_observations_single_file(self, single_pair_file):
        batch = BatchUTAS([single_pair_file])
        assert batch.num_observations == 2

    def test_num_observations_two_files(self, two_files_single_pair):
        batch = BatchUTAS(two_files_single_pair)
        assert batch.num_observations == 2

    def test_station_pairs(self, single_pair_file):
        batch = BatchUTAS([single_pair_file])
        assert batch.station_pairs == [("KATHERINE", "HOBART")]

    def test_station_names(self, single_pair_file):
        batch = BatchUTAS([single_pair_file])
        assert batch.station_names == {"KATHERINE", "HOBART"}

    def test_num_station_pairs(self, single_pair_file):
        batch = BatchUTAS([single_pair_file])
        assert batch.num_station_pairs == 1

    def test_num_station_pairs_multi(self, single_pair_file, second_pair_file):
        batch = BatchUTAS([single_pair_file, second_pair_file])
        assert batch.num_station_pairs == 2

    def test_station_names_multi(self, single_pair_file, second_pair_file):
        batch = BatchUTAS([single_pair_file, second_pair_file])
        assert batch.station_names == {"KATHERINE", "HOBART", "STATION_A", "STATION_B"}

    def test_metadata_fields(self, single_pair_file):
        """All metadata fields populated correctly."""
        batch = BatchUTAS([single_pair_file])
        meta = batch.get_metadata()
        assert meta.target_id == "53365"
        assert meta.frequency == 2260760000
        assert meta.bandwidth == 0.0
        assert meta.sensor1_delay == 0.0
        assert meta.sensor2_delay == 0.0
        assert meta.data_mode == "REAL"
        assert meta.origin == "Converted on Bingus"
        assert meta.source == "Unknown"
        assert meta.ucts == 0

    def test_get_observations_success(self, single_pair_file):
        batch = BatchUTAS([single_pair_file])
        obs = batch.get_observations_for_station_pair(("KATHERINE", "HOBART"))
        assert len(obs) == 2
        assert len(obs.tdoa) == 2
        assert len(obs.fdoa) == 2
        assert len(obs.epochs) == 2

    def test_get_observations_invalid_pair(self, single_pair_file):
        batch = BatchUTAS([single_pair_file])
        with pytest.raises(RuntimeError, match="Station pair .* not found"):
            batch.get_observations_for_station_pair(("NONEXISTENT", "ALSO_NONEXISTENT"))

    def test_get_all_observations(self, single_pair_file):
        batch = BatchUTAS([single_pair_file])
        all_obs = batch.get_all_observations_by_station_pair()
        assert len(all_obs) == 1
        assert ("KATHERINE", "HOBART") in all_obs
        assert len(all_obs[("KATHERINE", "HOBART")]) == 2

    def test_two_files_same_station_pair(self, two_files_single_pair):
        batch = BatchUTAS(two_files_single_pair)
        assert batch.num_observations == 2
        assert batch.num_station_pairs == 1
        assert len(batch.get_observations_for_station_pair(("KATHERINE", "HOBART"))) == 2

    def test_multi_file_multi_pair(self, single_pair_file, second_pair_file):
        batch = BatchUTAS([single_pair_file, second_pair_file])
        assert batch.num_observations == 3
        assert batch.num_station_pairs == 2
        assert ("KATHERINE", "HOBART") in batch.station_pairs
        assert ("STATION_A", "STATION_B") in batch.station_pairs

    def test_dict_with_observations_key(self, tmp_path):
        """Loading from dict with 'observations' key (alternate JSON structure)."""
        file_path = tmp_path / "dict_structure.json"
        _write_json(str(file_path), {"observations": SINGLE_STATION_PAIR_OBS})
        batch = BatchUTAS([str(file_path)])
        assert batch.num_observations == 2

    def test_satno_as_string(self, tmp_path):
        """satNo can be a string."""
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

    def test_optional_fields_use_defaults(self, tmp_path):
        """Missing optional fields fall back to default values."""
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
        meta = batch.get_metadata()
        assert meta.bandwidth == 0.0
        assert meta.sensor1_delay == 0.0
        assert meta.sensor2_delay == 0.0
        assert meta.data_mode == ""
        assert meta.origin == ""
        assert meta.source == ""


# =========================================================================
# Error handling
# =========================================================================


class TestBatchUTASErrors:
    """BatchUTAS must raise RuntimeError for invalid input."""

    def test_empty_file_raises_error(self, tmp_path):
        file_path = tmp_path / "empty.json"
        _write_json(str(file_path), [])
        with pytest.raises(RuntimeError, match="No observations found"):
            BatchUTAS([str(file_path)])

    def test_invalid_json_raises_error(self, tmp_path):
        file_path = tmp_path / "invalid.json"
        with open(file_path, "w") as f:
            f.write("{invalid json content")
        with pytest.raises(RuntimeError, match="JSON parse error"):
            BatchUTAS([str(file_path)])

    def test_wrong_structure_raises_error(self, tmp_path):
        file_path = tmp_path / "wrong_structure.json"
        _write_json(str(file_path), {"key": "value"})
        with pytest.raises(RuntimeError, match="Unexpected JSON structure"):
            BatchUTAS([str(file_path)])

    def test_missing_required_field(self, tmp_path):
        bad_data = [
            {
                "frequency": 2260760000,
                "origSensorId1": "A",
                "origSensorId2": "B",
                "satNo": 12345,
            }
        ]
        file_path = tmp_path / "bad.json"
        _write_json(str(file_path), bad_data)
        with pytest.raises(RuntimeError):
            BatchUTAS([str(file_path)])

    def test_multiple_targets_in_one_file(self, tmp_path):
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
        _write_json(str(file_path), data)
        with pytest.raises(RuntimeError, match="Multiple targets detected"):
            BatchUTAS([str(file_path)])

    def test_multiple_targets_across_files(self, tmp_path):
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


# =========================================================================
# Geodetic conversion (pure function, no side effects)
# =========================================================================


def test_convert_to_tudat_geodetic(single_pair_file):
    """Degrees → radians, altitude in metres unchanged."""
    batch = BatchUTAS([single_pair_file])
    position = {"altitude": 100.0, "latitude": 45.0, "longitude": -90.0}
    result = batch._convert_to_tudat_geodetic(position)
    assert result[0] == pytest.approx(100.0)
    assert result[1] == pytest.approx(45.0 * np.pi / 180.0)
    assert result[2] == pytest.approx(-90.0 * np.pi / 180.0)


# =========================================================================
# Integration: to_tudat with real Tudat objects
# =========================================================================


@pytest.mark.slow
def test_to_tudat_returns_valid_collection(single_pair_file):
    """Round-trip: parse two observations, convert to Tudat ObservationCollection.

    Uses real SPICE kernels and body creation so we validate that the ground
    stations, link definitions, and observation sets are structurally correct.
    """
    from tudatpy.interface import spice
    from tudatpy.dynamics import environment_setup
    from tudatpy.estimation import observations
    from tudatpy.estimation.observable_models_setup import model_settings

    spice.load_standard_kernels()

    body_settings = environment_setup.get_default_body_settings(["Earth"], "SSB", "J2000")
    bodies = environment_setup.create_system_of_bodies(body_settings)

    batch = BatchUTAS([single_pair_file])
    obs_collection = batch.to_tudat(bodies, target_name_override="KPLO")

    # --- structural checks ---
    assert isinstance(obs_collection, observations.ObservationCollection)

    # 1 station pair × 2 observable types = 2 observation sets
    obs_sets = obs_collection.get_single_observation_sets()
    assert len(obs_sets) == 2

    # --- check that the right observable types are present ---
    types_in_collection = obs_collection.get_observable_types()
    assert model_settings.differenced_time_of_arrival_type in types_in_collection
    assert model_settings.differenced_frequency_of_arrival_type in types_in_collection

    # --- each set should have 2 observations (matching input) ---
    concatenated_times = obs_collection.concatenated_times
    assert len(concatenated_times) == 4  # 2 obs × 2 types

    # --- ground stations should have been created ---
    earth = bodies.get_body("Earth")
    assert "KATHERINE" in earth.ground_station_list
    assert "HOBART" in earth.ground_station_list
