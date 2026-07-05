import numpy as np
import pandas as pd
import pytest

from tudatpy.data.jpl_radar import JPLRadarQuery
from tudatpy.data.radar import (
    DOPPLER_OBSERVABLE,
    RANGE_OBSERVABLE,
    RadarTrackingData,
    radar_tracking_data_from_mpc_table,
)
from tudatpy.dynamics import environment_setup
from tudatpy.estimation.observable_models_setup import model_settings
from tudatpy.constants import SPEED_OF_LIGHT
from tudatpy.interface import spice

WGS84_EQUATORIAL_RADIUS = 6378137.0
WGS84_FLATTENING = 1.0 / 298.257223563

spice.load_standard_kernels()


def _jpl_radar_response():
    return {
        "fields": [
            "epoch",
            "value",
            "sigma",
            "units",
            "freq",
            "xmit",
            "rcvr",
            "bp",
        ],
        "data": [
            [
                "2005-01-29 00:00:00",
                "4000.0",
                "0.025",
                "us",
                "8560.0",
                "-14",
                "-14",
                "C",
            ],
            [
                "2005-01-29 00:00:01",
                "-3.5",
                "0.1",
                "Hz",
                "8560.0",
                "-14",
                "-14",
                "C",
            ],
            [
                "2005-01-29 00:00:02",
                "123.0",
                "1.0",
                "us",
                "8560.0",
                "-14",
                "-14",
                "S",
            ],
        ],
    }


def test_jpl_radar_query_returns_generic_radar_tracking_data(monkeypatch):
    query = JPLRadarQuery("99942")
    monkeypatch.setattr(query, "_fetch_json", _jpl_radar_response)

    radar_data = query.to_radar_tracking_data(target_body="99942")
    table = radar_data.table

    assert len(table) == 2
    assert "RA" not in table
    assert "DEC" not in table
    assert set(table["observable_type"]) == {RANGE_OBSERVABLE, DOPPLER_OBSERVABLE}
    assert set(table["transmitter"]) == {"JPL:-14"}
    assert set(table["receiver"]) == {"JPL:-14"}

    range_row = table.loc[table["observable_type"] == RANGE_OBSERVABLE].iloc[0]
    assert range_row["value"] == pytest.approx(SPEED_OF_LIGHT * 4000.0e-6)
    assert range_row["sigma"] == pytest.approx(SPEED_OF_LIGHT * 0.025e-6)

    doppler_row = table.loc[table["observable_type"] == DOPPLER_OBSERVABLE].iloc[0]
    assert doppler_row["value"] == pytest.approx(8560.0e6 - 3.5)
    assert doppler_row["sigma"] == pytest.approx(0.1)
    assert doppler_row["transmitter_frequency_hz"] == pytest.approx(8560.0e6)


def test_jpl_radar_query_can_use_mpc_station_compatibility_mode(monkeypatch):
    query = JPLRadarQuery("99942")
    monkeypatch.setattr(query, "_fetch_json", _jpl_radar_response)

    table = query.to_radar_tracking_data(station_id_mode="mpc").table

    assert set(table["transmitter"]) == {"253"}
    assert set(table["receiver"]) == {"253"}


def test_radar_tracking_data_rejects_angular_columns():
    table = pd.DataFrame(
        {
            "target_body": ["99942"],
            "epoch_seconds_UTC": [0.0],
            "epoch_seconds_TDB": [0.0],
            "observable_type": [RANGE_OBSERVABLE],
            "value": [1.0],
            "sigma": [1.0],
            "transmitter": ["JPL:-14"],
            "receiver": ["JPL:-14"],
            "target_point": ["C"],
            "transmitter_frequency_hz": [np.nan],
            "source": ["JPL"],
            "RA": [np.nan],
        }
    )

    with pytest.raises(ValueError, match="must not contain angular astrometry"):
        RadarTrackingData(table)


def test_mpc_radar_table_adapter_returns_generic_schema_without_angular_columns():
    mpc_table = pd.DataFrame(
        {
            "number": ["99942"],
            "epoch_seconds_UTC": [1.0],
            "epoch_seconds_TDB": [2.0],
            "radar_target_point": ["C"],
            "radar_transmitter": ["251"],
            "radar_receiver": ["251"],
            "radar_range": [123.0],
            "radar_range_sigma": [4.0],
            "RA": [np.nan],
            "DEC": [np.nan],
        }
    )

    radar_data = radar_tracking_data_from_mpc_table(mpc_table)
    table = radar_data.table

    assert len(table) == 1
    assert "RA" not in table
    assert "DEC" not in table
    assert table.iloc[0]["target_body"] == "99942"
    assert table.iloc[0]["observable_type"] == RANGE_OBSERVABLE
    assert table.iloc[0]["value"] == 123.0
    assert table.iloc[0]["sigma"] == 4.0


def test_radar_tracking_data_to_tudat_adds_range_and_doppler_sets():
    table = pd.DataFrame(
        {
            "target_body": ["99942", "99942"],
            "epoch_seconds_UTC": [0.0, 1.0],
            "epoch_seconds_TDB": [0.0, 1.0],
            "observable_type": [RANGE_OBSERVABLE, DOPPLER_OBSERVABLE],
            "value": [1.0e6, 8560.0e6 - 3.5],
            "sigma": [10.0, 0.1],
            "transmitter": ["JPL:-14", "JPL:-14"],
            "receiver": ["JPL:-14", "JPL:-14"],
            "target_point": ["C", "C"],
            "transmitter_frequency_hz": [np.nan, 8560.0e6],
            "source": ["JPL", "JPL"],
        }
    )
    radar_data = RadarTrackingData(table)

    body_settings = environment_setup.get_default_body_settings(["Earth"], "SSB", "J2000")
    body_settings.get("Earth").shape_settings = environment_setup.shape.oblate_spherical(
        WGS84_EQUATORIAL_RADIUS,
        WGS84_FLATTENING,
    )
    bodies = environment_setup.create_system_of_bodies(body_settings)

    observation_dataset = radar_data.to_tudat(bodies)
    observable_types = {
        metadata.observable_type for metadata in observation_dataset.observation_set_metadata
    }

    assert "JPL:-14" in bodies.get("Earth").ground_station_list
    assert bodies.does_body_exist("99942")
    assert model_settings.n_way_range_type in observable_types
    assert model_settings.doppler_measured_frequency_type in observable_types
