import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.table import Table

from tudatpy.astro import time_representation
from tudatpy.constants import SPEED_OF_LIGHT
import tudatpy.data_input.tracking_data.jpl_radar.jpl_radar as jpl_radar_backend
from tudatpy.data_input.tracking_data.jpl_radar import (
    JPLRadarQuery,
    get_available_radar_targets,
)
from tudatpy.data_input.tracking_data.mpc import BatchMPC
from tudatpy.data_input.tracking_data.obs_80_cols.parsers import parse_80cols_data
from tudatpy.data_input.tracking_data.radar_utilities import (
    DOPPLER_OBSERVABLE,
    RADAR_COLUMNS,
    RADAR_TABLE_META_KEY,
    RANGE_OBSERVABLE,
    empty_radar_table,
    radar_data_to_tracking_data,
)

_JPL_RESPONSE = {
    "signature": {"version": "1.1", "source": "NASA/JPL Small-Body Radar Astrometry API"},
    "count": "3",
    "fields": ["des", "epoch", "value", "sigma", "units", "freq", "rcvr", "xmit", "bp"],
    "data": [
        [
            "1997 WQ23",
            "2013-11-22 01:53:00",
            "111982347.7",
            "1.000",
            "us",
            "2380",
            "-1",
            "-1",
            "C",
        ],
        [
            "1997 WQ23",
            "2013-11-19 00:38:00",
            "96148022.94",
            "4.000",
            "us",
            "2380",
            "-1",
            "-1",
            "C",
        ],
        [
            "1997 WQ23",
            "2013-11-19 00:25:00",
            "-124458.2215",
            "0.300",
            "Hz",
            "2380",
            "-1",
            "-1",
            "C",
        ],
    ],
    "coords": {
        "-1": {
            "alt_units": "m",
            "altitude": "450.56",
            "latitude": "18.34420",
            "longitude": "293.2473",
            "name": "Arecibo",
        }
    },
}


def _radar_table(**columns):
    defaults = {
        "target_body": "433",
        "observable_type": DOPPLER_OBSERVABLE,
        "value": 2380.0e6,
        "sigma": 0.1,
        "transmitter": "251",
        "receiver": "251",
        "target_point": "C",
        "transmitter_frequency_hz": 2380.0e6,
        "source": "test",
    }
    table = pd.DataFrame(columns)
    for column, value in defaults.items():
        if column not in table:
            table[column] = value
    return table[RADAR_COLUMNS]


def _parsed_radar(lines):
    return parse_80cols_data(lines).meta[RADAR_TABLE_META_KEY]


def _roving_observer_pair():
    optical = "     K25F22A  C2025 10 13.24277 00 20 45.76 +25 53 06.1          18.3 RrET147718"
    first = optical[:14] + "V" + optical[15:77] + "247"
    second = (
        optical[:14]
        + "v"
        + optical[15:32]
        + f"1 {243.1105:10.4f} {35.4259:+10.4f} {1001:5d}"
        + " " * 16
        + "247"
    )
    assert len(first) == len(second) == 80
    return optical, first, second


def test_empty_radar_table_uses_canonical_columns():
    assert list(empty_radar_table().columns) == RADAR_COLUMNS


def test_mpc_radar_delay_record(mpc_radar_pair):
    row = _parsed_radar(mpc_radar_pair(delay_us=4000.0, delay_sigma_us=0.25)).iloc[0]
    assert row["observable_type"] == RANGE_OBSERVABLE
    assert row["value"] == pytest.approx(SPEED_OF_LIGHT * 4000.0e-6)
    assert row["sigma"] == pytest.approx(SPEED_OF_LIGHT * 0.25e-6)
    assert (row["target_body"], row["transmitter"], row["receiver"]) == (
        "433",
        "251",
        "251",
    )


def test_mpc_radar_doppler_record(mpc_radar_pair):
    row = _parsed_radar(mpc_radar_pair(doppler_hz=-124458.2215, doppler_sigma_hz=0.3)).iloc[0]
    assert row["observable_type"] == DOPPLER_OBSERVABLE
    assert row["value"] == pytest.approx(2380.0e6 - 124458.2215)
    assert row["sigma"] == pytest.approx(0.3)
    assert row["transmitter_frequency_hz"] == pytest.approx(2380.0e6)


def test_mpc_radar_record_with_delay_and_doppler_gives_two_rows(mpc_radar_pair):
    table = _parsed_radar(
        mpc_radar_pair(
            delay_us=4000.0,
            delay_sigma_us=0.25,
            doppler_hz=-3.5,
            doppler_sigma_hz=0.1,
        )
    )
    assert set(table["observable_type"]) == {RANGE_OBSERVABLE, DOPPLER_OBSERVABLE}


def test_mpc_radar_epoch_is_rounded_to_integral_second(mpc_radar_pair):
    table = _parsed_radar(
        mpc_radar_pair(
            date="1990 07 15.326389",
            delay_us=4000.0,
            delay_sigma_us=0.25,
        )
    )
    expected = time_representation.date_time_components_to_epoch(1990, 7, 15, 7, 50, 0.0)
    assert table["epoch_seconds_UTC"].iloc[0] == expected


def test_mpc_radar_keeps_surface_bounce_point_but_does_not_convert_it(mpc_radar_pair):
    table = _parsed_radar(
        mpc_radar_pair(
            delay_us=4000.0,
            delay_sigma_us=0.25,
            bounce_point="S",
        )
    )
    assert list(table["target_point"]) == ["S"]
    tracking_data, supplementary_data = radar_data_to_tracking_data(table)
    assert len(tracking_data) == len(supplementary_data) == 0


def test_roving_observer_pair_is_not_parsed_as_radar():
    optical, roving_first, roving_second = _roving_observer_pair()
    parsed = parse_80cols_data([optical, roving_first, roving_second])
    assert RADAR_TABLE_META_KEY not in parsed.meta
    assert len(parsed) == 1


def test_real_mpc_radar_lines_parse():
    # Published 80-column record from https://projectpluto.com/radar/99942.htm.
    lines = (Path(__file__).parent / "fixtures" / "mpc_radar_sample.txt").read_text().splitlines()
    table = _parsed_radar(lines)
    assert set(table["observable_type"]) == {RANGE_OBSERVABLE, DOPPLER_OBSERVABLE}


def test_jpl_radar_query_returns_canonical_radar_data(monkeypatch):
    monkeypatch.setattr(jpl_radar_backend, "_query", lambda params, timeout: _JPL_RESPONSE)
    table = JPLRadarQuery("1997 WQ23").to_radar_data()

    assert len(table) == 3
    assert set(table["transmitter"]) == set(table["receiver"]) == {"251"}
    delay = (
        table[table["observable_type"] == RANGE_OBSERVABLE].sort_values("epoch_seconds_UTC").iloc[0]
    )
    assert delay["value"] == pytest.approx(SPEED_OF_LIGHT * 96148022.94e-6)
    assert delay["sigma"] == pytest.approx(SPEED_OF_LIGHT * 4.0e-6)
    doppler = table[table["observable_type"] == DOPPLER_OBSERVABLE].iloc[0]
    assert doppler["value"] == pytest.approx(2380.0e6 - 124458.2215)


def test_jpl_radar_query_epoch_filter_and_station_positions(monkeypatch):
    monkeypatch.setattr(jpl_radar_backend, "_query", lambda params, timeout: _JPL_RESPONSE)
    query = JPLRadarQuery("1997 WQ23")
    table = query.to_radar_data(epoch_start=datetime.datetime(2013, 11, 20))
    assert len(table) == 1
    np.testing.assert_allclose(
        query.station_geodetic_positions()["251"],
        [450.56, np.deg2rad(18.3442), np.deg2rad(293.2473)],
    )


def test_jpl_empty_response_gives_empty_table(monkeypatch):
    empty = {"count": 0, "signature": _JPL_RESPONSE["signature"]}
    monkeypatch.setattr(jpl_radar_backend, "_query", lambda params, timeout: empty)
    assert JPLRadarQuery("0000").to_radar_data().empty
    assert get_available_radar_targets() == []


@pytest.mark.parametrize(
    "jpl_code, mpc_code",
    [
        ("-1", "251"),
        ("-2", "254"),
        ("-9", "256"),
        ("-13", "252"),
        ("-14", "253"),
        ("-25", "257"),
        ("-38", "255"),
        ("-73", "259"),
        ("-99", "JPL:-99"),
    ],
)
def test_jpl_station_codes_map_to_mpc_codes(jpl_code, mpc_code):
    assert jpl_radar_backend._station_id(jpl_code) == mpc_code


def test_radar_data_converts_with_weights_and_link_ends():
    table = _radar_table(
        epoch_seconds_UTC=[0.0, 1.0],
        observable_type=[RANGE_OBSERVABLE, DOPPLER_OBSERVABLE],
        value=[1.0e6, 8560.0e6 - 3.5],
        sigma=[10.0, 0.1],
        transmitter_frequency_hz=[np.nan, 8560.0e6],
    )
    tracking_data, supplementary_data = radar_data_to_tracking_data(table)

    assert [data.observable_type for data in tracking_data] == [
        RANGE_OBSERVABLE,
        DOPPLER_OBSERVABLE,
    ]
    assert tracking_data[0].link_ends == [
        (("Earth", "251"), "transmitter"),
        (("433", ""), "reflector_1"),
        (("Earth", "251"), "receiver"),
    ]
    np.testing.assert_allclose(np.concatenate(tracking_data[0].get_observation_weights()), [1.0e-2])
    np.testing.assert_allclose(np.concatenate(tracking_data[1].get_observation_weights()), [1.0e2])
    assert [(item.body_name, item.reference_point_name) for item in supplementary_data] == [
        ("Earth", "251")
    ]


def test_frequency_history_switches_between_tracks():
    table = _radar_table(
        epoch_seconds_UTC=[0.0, 60.0, 86400.0, 86460.0],
        transmitter_frequency_hz=[2380.0e6, 2380.0e6, 8560.0e6, 8560.0e6],
        value=[2380.0e6, 2380.0e6, 8560.0e6, 8560.0e6],
    )
    tracking_data, supplementary_data = radar_data_to_tracking_data(table)

    assert len(tracking_data) == 2
    history = supplementary_data[0].frequency_supplementary_data[0].frequency_history
    assert history == {-86400.0: 2380.0e6, 43230.0: 8560.0e6}


def test_batchmpc_get_satellite_state_history():
    batch = BatchMPC()
    batch._table = pd.DataFrame(
        {
            "observatory": ["C51", "C51", "C51"],
            "epoch_seconds_TDB": [0.0, 10.0, 20.0],
            "spacecraft_position_x": [0.0, 10.0, 20.0],
            "spacecraft_position_y": [10.0, 10.0, 10.0],
            "spacecraft_position_z": [0.0, -10.0, -20.0],
        }
    )

    state_history = batch.get_satellite_state_history("C51")
    epochs = sorted(state_history)
    states = np.array([state_history[epoch] for epoch in epochs])
    assert epochs == [0.0, 10.0, 20.0]
    np.testing.assert_allclose(
        states[:, :3],
        [[0.0, 10.0, 0.0], [10.0, 10.0, -10.0], [20.0, 10.0, -20.0]],
    )
    np.testing.assert_allclose(states[:, 3:], [[1.0, 0.0, -1.0]] * 3)


def test_batchmpc_mpc80_path_loads_space_astrometry_and_radar(monkeypatch, mpc_radar_pair):
    eros_observation = (
        "00433         S2021 06 07.42640918 08 15.401-41 22 02.35         12.0 V      500"
    )
    eros_parallax = (
        "00433         s2021 06 07.4264091 -198301.940 +198171.039 +56287.9850   ~6oMXC57"
    )
    radar_lines = mpc_radar_pair(
        delay_us=4000.0,
        delay_sigma_us=0.25,
        doppler_hz=-3.5,
        doppler_sigma_hz=0.1,
        transmitter="253",
        receiver="253",
    )

    def fake_get_observations(code, **kwargs):
        assert kwargs["get_mpcformat"] is True
        return Table({"obs": [eros_observation + eros_parallax, *radar_lines]})

    monkeypatch.setattr("astroquery.mpc.MPC.get_observations", fake_get_observations)
    batch = BatchMPC()
    batch.get_observations([433], use_mpc80_format=True)

    assert len(batch.table) == 1 and len(batch.radar_table) == 2 and batch.size == 3
    assert batch.epoch_start == batch.radar_table["epoch_seconds_UTC"].min()
    assert batch.epoch_end == batch.table["epoch_seconds_UTC"].max()
    np.testing.assert_allclose(batch.table["spacecraft_position_x"].iloc[0], -198301940.0)

    radar_only = batch.filter(observation_types=["R"], in_place=False)
    assert len(radar_only.table) == 0 and len(radar_only.radar_table) == 2

    tracking_data, supplementary_data = batch.to_tracking_dataset()
    assert {data.observable_type for data in tracking_data} == {
        "AngularPosition",
        RANGE_OBSERVABLE,
        DOPPLER_OBSERVABLE,
    }
    assert {item.body_name for item in supplementary_data} == {"500", "Earth"}
