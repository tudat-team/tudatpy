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
            "alt_units": "km",
            "altitude": "0.45334",
            "latitude": "18.3442199",
            "longitude": "293.2473068",
            "name": "Arecibo (305-m, 1963 to 2020)",
        }
    },
}


def _implicit_decimal_field(value, integer_width, fraction_width=4, signed=False):
    """Return a fixed-width MPC field with an implied decimal point."""
    if value is None:
        return " " * (integer_width + fraction_width)
    width = integer_width + fraction_width + 1
    return f"{value:{'+' if signed else ''}0{width}.{fraction_width}f}".replace(".", "")


def _mpc_radar_pair(
    number="00433",
    date="1990 07 15.326389",
    delay_us=None,
    delay_sigma_us=None,
    doppler_hz=None,
    doppler_sigma_hz=None,
    frequency_mhz=2380.0,
    transmitter="251",
    receiver="251",
    bounce_point="C",
):
    """Return two 80-column MPC radar records."""
    head = f"{number:<5}{'':7}  "
    tail = f"{transmitter:>3}{'':6}{receiver:>3}"
    first = (
        head
        + "R"
        + date
        + _implicit_decimal_field(delay_us, 11)
        + _implicit_decimal_field(doppler_hz, 11, signed=True)
        + _implicit_decimal_field(frequency_mhz, 5, 1)
        + tail
    )
    second = (
        head
        + "r"
        + date
        + bounce_point
        + _implicit_decimal_field(delay_sigma_us, 10)
        + _implicit_decimal_field(doppler_sigma_hz, 11)
        + " " * 6
        + tail
    )
    # The factory must preserve the MPC fixed-width record size.
    assert len(first) == len(second) == 80
    return [first, second]


@pytest.fixture
def mpc_radar_pair():
    """Return a factory for MPC radar record pairs."""
    return _mpc_radar_pair


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
    # The helper must produce valid fixed-width MPC records.
    assert len(first) == len(second) == 80
    return optical, first, second


def test_empty_radar_table_uses_canonical_columns():
    # Check that even an empty radar result has the schema expected by all readers.
    assert list(empty_radar_table().columns) == RADAR_COLUMNS


def test_mpc_radar_delay_record(mpc_radar_pair):
    # Check the MPC delay-field interpretation, uncertainty conversion and station IDs.
    row = _parsed_radar(mpc_radar_pair(delay_us=4000.0, delay_sigma_us=0.25)).iloc[0]

    # Delay and its sigma are round-trip light times in microseconds and must
    # become an N-way range and uncertainty in metres.
    assert row["observable_type"] == RANGE_OBSERVABLE
    assert row["value"] == pytest.approx(SPEED_OF_LIGHT * 4000.0e-6)
    assert row["sigma"] == pytest.approx(SPEED_OF_LIGHT * 0.25e-6)

    # Object and station identifiers must be unpacked and normalized to MPC codes.
    assert (row["target_body"], row["transmitter"], row["receiver"]) == (
        "433",
        "251",
        "251",
    )


def test_mpc_radar_doppler_record(mpc_radar_pair):
    # Check the MPC Doppler columns and conversion from shift to received frequency.
    row = _parsed_radar(mpc_radar_pair(doppler_hz=-124458.2215, doppler_sigma_hz=0.3)).iloc[0]

    # The parser must preserve the source frequency and add the signed Doppler
    # shift while leaving the Doppler uncertainty in hertz.
    assert row["observable_type"] == DOPPLER_OBSERVABLE
    assert row["value"] == pytest.approx(2380.0e6 - 124458.2215)
    assert row["sigma"] == pytest.approx(0.3)
    assert row["transmitter_frequency_hz"] == pytest.approx(2380.0e6)


def test_mpc_radar_record_with_delay_and_doppler_gives_two_rows(mpc_radar_pair):
    # Check that one MPC R/r pair carrying both measurements emits two observations.
    table = _parsed_radar(
        mpc_radar_pair(
            delay_us=4000.0,
            delay_sigma_us=0.25,
            doppler_hz=-3.5,
            doppler_sigma_hz=0.1,
        )
    )

    # Both Tudat observable types must be present exactly as separate table rows.
    assert set(table["observable_type"]) == {RANGE_OBSERVABLE, DOPPLER_OBSERVABLE}


def test_mpc_radar_epoch_is_rounded_to_integral_second(mpc_radar_pair):
    # Check the MPC requirement that radar reception epochs denote integral seconds.
    table = _parsed_radar(
        mpc_radar_pair(
            date="1990 07 15.326389",
            delay_us=4000.0,
            delay_sigma_us=0.25,
        )
    )
    expected = time_representation.date_time_components_to_epoch(1990, 7, 15, 7, 50, 0.0)

    # The rounded parser epoch must equal the exact calendar epoch, not the raw
    # fractional-day value before rounding.
    assert table["epoch_seconds_UTC"].iloc[0] == expected


def test_mpc_radar_keeps_surface_bounce_point_but_does_not_convert_it(mpc_radar_pair):
    # Check that unsupported surface-bounce data remain inspectable but are not converted.
    table = _parsed_radar(
        mpc_radar_pair(
            delay_us=4000.0,
            delay_sigma_us=0.25,
            bounce_point="S",
        )
    )

    # Parsing must retain the source bounce-point label.
    assert list(table["target_point"]) == ["S"]
    tracking_data, supplementary_data = radar_data_to_tracking_data(table)

    # The default conversion accepts centre-of-mass observations only.
    assert len(tracking_data) == len(supplementary_data) == 0


def test_roving_observer_pair_is_not_parsed_as_radar():
    # Check that optical V/v roving-observer records are not mistaken for radar data.
    optical, roving_first, roving_second = _roving_observer_pair()
    parsed = parse_80cols_data([optical, roving_first, roving_second])

    # No radar metadata may be created, while the ordinary optical row remains.
    assert RADAR_TABLE_META_KEY not in parsed.meta
    assert len(parsed) == 1


def test_real_mpc_radar_lines_parse():
    # Check the column layout against a published, non-synthetic MPC radar pair.
    # Published 80-column record from https://projectpluto.com/radar/99942.htm.
    lines = (Path(__file__).parent / "fixtures" / "mpc_radar_sample.txt").read_text().splitlines()
    table = _parsed_radar(lines)

    # The real pair contains both a delay and a Doppler measurement.
    assert set(table["observable_type"]) == {RANGE_OBSERVABLE, DOPPLER_OBSERVABLE}


def test_jpl_radar_query_returns_canonical_radar_data(monkeypatch):
    # Check conversion of a representative JPL API response to the shared schema.
    monkeypatch.setattr(jpl_radar_backend, "_query", lambda params, timeout: _JPL_RESPONSE)
    table = JPLRadarQuery("1997 WQ23").to_radar_data()

    # All API measurements must survive, with known JPL stations mapped to MPC codes.
    assert len(table) == 3
    assert set(table["transmitter"]) == set(table["receiver"]) == {"251"}

    # Delay values and sigmas must be converted from microseconds to metres.
    delay = (
        table[table["observable_type"] == RANGE_OBSERVABLE].sort_values("epoch_seconds_UTC").iloc[0]
    )
    assert delay["value"] == pytest.approx(SPEED_OF_LIGHT * 96148022.94e-6)
    assert delay["sigma"] == pytest.approx(SPEED_OF_LIGHT * 4.0e-6)

    # Doppler shifts must be added to the transmitted frequency.
    doppler = table[table["observable_type"] == DOPPLER_OBSERVABLE].iloc[0]
    assert doppler["value"] == pytest.approx(2380.0e6 - 124458.2215)


def test_jpl_radar_query_epoch_filter_and_station_positions(monkeypatch):
    # Check datetime filtering and conversion of JPL geodetic station coordinates.
    monkeypatch.setattr(jpl_radar_backend, "_query", lambda params, timeout: _JPL_RESPONSE)
    query = JPLRadarQuery("1997 WQ23")
    table = query.to_radar_data(epoch_start=datetime.datetime(2013, 11, 20))

    # Only the measurement inside the requested epoch interval must remain.
    assert len(table) == 1

    # Altitude stays in metres and angular coordinates are converted to radians.
    np.testing.assert_allclose(
        query.station_geodetic_positions()["251"],
        [453.34, np.deg2rad(18.3442199), np.deg2rad(293.2473068)],
    )


def test_jpl_empty_response_gives_empty_table(monkeypatch):
    # Check that valid empty API responses are handled by both public query paths.
    empty = {"count": 0, "signature": _JPL_RESPONSE["signature"]}
    monkeypatch.setattr(jpl_radar_backend, "_query", lambda params, timeout: empty)

    # Neither a target query nor the catalogue query may fail or invent rows.
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
    # Check every supported JPL-to-MPC station mapping and the unknown-code fallback.
    assert jpl_radar_backend._station_id(jpl_code) == mpc_code


def test_radar_data_converts_with_weights_and_link_ends():
    # Check conversion of canonical range and Doppler rows to weighted TrackingData.
    table = _radar_table(
        epoch_seconds_UTC=[0.0, 1.0],
        observable_type=[RANGE_OBSERVABLE, DOPPLER_OBSERVABLE],
        value=[1.0e6, 8560.0e6 - 3.5],
        sigma=[10.0, 0.1],
        transmitter_frequency_hz=[np.nan, 8560.0e6],
    )
    tracking_data, supplementary_data = radar_data_to_tracking_data(table)

    # Range and Doppler rows must become separate observation sets in input order.
    assert [data.observable_type for data in tracking_data] == [
        RANGE_OBSERVABLE,
        DOPPLER_OBSERVABLE,
    ]

    # The N-way link must connect transmitter, passive target and receiver.
    assert tracking_data[0].link_ends == [
        (("Earth", "251"), "transmitter"),
        (("433", ""), "reflector_1"),
        (("Earth", "251"), "receiver"),
    ]

    # Each source sigma must be stored as the inverse-variance observation weight.
    np.testing.assert_allclose(np.concatenate(tracking_data[0].get_observation_weights()), [1.0e-2])
    np.testing.assert_allclose(np.concatenate(tracking_data[1].get_observation_weights()), [1.0e2])

    # Doppler conversion must produce frequency data for the transmitting station.
    assert [(item.body_name, item.reference_point_name) for item in supplementary_data] == [
        ("Earth", "251")
    ]


def test_frequency_history_switches_between_tracks():
    # Check grouping by frequency band and placement of frequency-history change points.
    table = _radar_table(
        epoch_seconds_UTC=[0.0, 60.0, 86400.0, 86460.0],
        transmitter_frequency_hz=[2380.0e6, 2380.0e6, 8560.0e6, 8560.0e6],
        value=[2380.0e6, 2380.0e6, 8560.0e6, 8560.0e6],
    )
    tracking_data, supplementary_data = radar_data_to_tracking_data(table)

    # S- and X-band observations must be split into separate tracking data sets.
    assert len(tracking_data) == 2
    history = supplementary_data[0].frequency_supplementary_data[0].frequency_history

    # The history starts one day early and switches midway between the two tracks.
    assert history == {-86400.0: 2380.0e6, 43230.0: 8560.0e6}


def test_batchmpc_get_satellite_state_history():
    # Check conversion of MPC spacecraft positions to a TDB Cartesian state history.
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

    # The state history must retain all source epochs in sorted order.
    assert epochs == [0.0, 10.0, 20.0]

    # Position components must be preserved without changing their frame or units.
    np.testing.assert_allclose(
        states[:, :3],
        [[0.0, 10.0, 0.0], [10.0, 10.0, -10.0], [20.0, 10.0, -20.0]],
    )

    # Finite differences must recover the constant velocity of this linear track.
    np.testing.assert_allclose(states[:, 3:], [[1.0, 0.0, -1.0]] * 3)


def test_batchmpc_mpc80_path_loads_space_astrometry_and_radar(monkeypatch, mpc_radar_pair):
    # Check the raw MPC80 backend with space astrometry and both radar observables.
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
        # The raw backend must explicitly request MPC-format records from astroquery.
        assert kwargs["get_mpcformat"] is True
        return Table({"obs": [eros_observation + eros_parallax, *radar_lines]})

    monkeypatch.setattr("astroquery.mpc.MPC.get_observations", fake_get_observations)
    batch = BatchMPC()
    batch.get_observations([433], use_mpc80_format=True)

    # Optical and radar rows must all contribute to the batch metadata and size.
    assert len(batch.table) == 1 and len(batch.radar_table) == 2 and batch.size == 3
    assert batch.epoch_start == batch.radar_table["epoch_seconds_UTC"].min()
    assert batch.epoch_end == batch.table["epoch_seconds_UTC"].max()

    # The paired parallax line must supply the optical receiver position in metres.
    np.testing.assert_allclose(batch.table["spacecraft_position_x"].iloc[0], -198301940.0)

    radar_only = batch.filter(observation_types=["R"], in_place=False)

    # Filtering by MPC radar type must keep both delay and Doppler, but no optical row.
    assert len(radar_only.table) == 0 and len(radar_only.radar_table) == 2

    tracking_data, supplementary_data = batch.to_tracking_dataset()

    # Conversion must create all three observable types and supplementary data
    # for both the spacecraft receiver and terrestrial transmitter.
    assert {data.observable_type for data in tracking_data} == {
        "AngularPosition",
        RANGE_OBSERVABLE,
        DOPPLER_OBSERVABLE,
    }
    assert {item.body_name for item in supplementary_data} == {"500", "Earth"}
