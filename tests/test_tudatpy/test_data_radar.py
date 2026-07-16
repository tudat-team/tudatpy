import numpy as np
import pandas as pd
import pytest
from astropy.table import Table

from tudatpy.constants import SPEED_OF_LIGHT
import tudatpy.data_input.tracking_data.jpl_radar.jpl_radar as jpl_radar_backend
from tudatpy.data_input.tracking_data.jpl_radar import (
    JPLRadarQuery,
    get_available_radar_targets,
)
from tudatpy.data_input.tracking_data.mpc import BatchMPC
from tudatpy.data_input.tracking_data.radar_utilities import (
    DOPPLER_OBSERVABLE,
    RADAR_COLUMNS,
    RANGE_OBSERVABLE,
    empty_radar_table,
    radar_data_to_tracking_data,
    radar_frequency_band_string_from_hz,
)
from tudatpy.data_input.tracking_data.radar_utilities.stations import (
    get_radar_station_geodetic_positions,
    radar_ground_station_settings,
)


def _jpl_radar_response():
    return {
        "fields": [
            "epoch",
            "des",
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
                "99942",
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
                "99942",
                "-3.5",
                "0.1",
                "Hz",
                "8560.0",
                "-14",
                "-14",
                "C",
            ],
        ],
    }


def test_empty_radar_table_uses_canonical_columns():
    table = empty_radar_table()

    assert tuple(table.columns) == RADAR_COLUMNS
    assert table.empty


def test_get_available_radar_targets_returns_sorted_unique_designations(monkeypatch):
    monkeypatch.setattr(
        jpl_radar_backend,
        "_fetch_json",
        lambda query_parameters, timeout, context: {
            "fields": ["des", "epoch"],
            "data": [
                ["99942", "2005-01-29 00:00:00"],
                ["433", "1975-01-01 00:00:00"],
                ["99942", "2005-01-29 00:00:01"],
            ],
        },
    )

    assert get_available_radar_targets() == ["433", "99942"]


def test_jpl_radar_query_returns_canonical_radar_data(monkeypatch):
    query = JPLRadarQuery("99942")
    monkeypatch.setattr(query, "_fetch_json", _jpl_radar_response)

    table = query.to_radar_data(target_body="99942")

    assert len(table) == 2
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


def test_radar_data_converts_to_tracking_and_supplementary_data():
    table = pd.DataFrame(
        {
            "target_body": ["99942", "99942"],
            "epoch_seconds_UTC": [0.0, 1.0],
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

    tracking_data, supplementary_data = radar_data_to_tracking_data(table)

    assert [data.observable_type for data in tracking_data] == [
        "NWayRange",
        "DopplerMeasuredFrequency",
    ]
    assert tracking_data[0].link_ends == [
        (("Earth", "JPL:-14"), "transmitter"),
        (("99942", ""), "reflector_1"),
        (("Earth", "JPL:-14"), "receiver"),
    ]
    assert len(supplementary_data) == 1
    assert supplementary_data[0].body_name == "Earth"
    assert supplementary_data[0].reference_point_name == "JPL:-14"


def test_radar_station_settings_can_include_all_known_stations():
    positions = get_radar_station_geodetic_positions()
    settings = radar_ground_station_settings()

    assert {"251", "253", "JPL:-14"}.issubset(positions)
    assert len(settings) == len(positions)


def test_radar_frequency_band_boundary_uses_ku_at_12_ghz():
    assert radar_frequency_band_string_from_hz(12.0e9) == "Ku-band"


def test_batchmpc_raw_mpc80_path_loads_space_astrometry_and_radar(monkeypatch):
    line_eros_valid = (
        "00433         S2021 06 07.42640918 08 15.401-41 22 02.35         12.0 V      500"
    )
    line_eros_parallax = (
        "00433         s2021 06 07.4264091 -198301.940 +198171.039 +56287.9850   ~6oMXC57"
    )

    first_radar_record = list(" " * 80)
    second_radar_record = list(" " * 80)

    def set_field(record, start, value):
        record[start : start + len(value)] = list(value)

    set_field(first_radar_record, 0, "00433")
    set_field(second_radar_record, 0, "00433")
    first_radar_record[14] = "R"
    second_radar_record[14] = "r"
    set_field(first_radar_record, 15, "2005")
    set_field(first_radar_record, 20, "01")
    set_field(first_radar_record, 23, "29.000000")
    set_field(first_radar_record, 32, "00000004000000")
    set_field(first_radar_record, 62, "085600")
    set_field(first_radar_record, 68, "253")
    set_field(first_radar_record, 77, "253")
    second_radar_record[32] = "C"
    set_field(second_radar_record, 33, "00000000000250")

    def fake_get_observations(code, **kwargs):
        assert kwargs["get_mpcformat"] is True
        return Table(
            {
                "obs": [
                    line_eros_valid + line_eros_parallax,
                    "".join(first_radar_record),
                    "".join(second_radar_record),
                ]
            }
        )

    monkeypatch.setattr("astroquery.mpc.MPC.get_observations", fake_get_observations)

    batch = BatchMPC()
    batch.get_observations([433], drop_misc_observations=False)

    assert len(batch.table) == 1
    assert len(batch.radar_table) == 1
    assert batch.size == 2
    np.testing.assert_allclose(batch.table["spacecraft_position_x"].iloc[0], -198301940.0)

    filtered_batch = batch.filter(observation_types=["R"], in_place=False)
    assert len(filtered_batch.table) == 0
    assert len(filtered_batch.radar_table) == 1
    assert filtered_batch.size == 1

    tracking_data, supplementary_data = batch.to_tracking_dataset()
    assert {data.observable_type for data in tracking_data} == {"AngularPosition", "NWayRange"}
    assert len(supplementary_data) == 1
    assert supplementary_data[0].body_name == "500"
