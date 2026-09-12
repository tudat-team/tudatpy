"""Exercise deprecated APIs using deterministic data and the real Tudat containers."""

import importlib
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy.table import Table
from astroquery.mpc import MPC

from tudatpy.dynamics import environment, environment_setup
from tudatpy.estimation.observations import ObservationCollection


def legacy_symbol(module, name):
    with pytest.warns(DeprecationWarning, match="is deprecated"):
        return getattr(importlib.import_module(module), name)


def earth_bodies():
    settings = environment_setup.BodyListSettings("SSB", "J2000")
    settings.add_empty_settings("Earth")
    settings.get("Earth").shape_settings = environment_setup.shape.spherical(6378137.0)
    return environment_setup.create_system_of_bodies(settings)


@pytest.fixture
def legacy_mpc_table(monkeypatch):
    stations = Table({"Code": ["500"], "Longitude": [0.0], "cos": [1.0], "sin": [0.0]})
    monkeypatch.setattr(MPC, "get_observatory_codes", lambda: stations)
    return pd.DataFrame(
        {
            "number": ["433"],
            "epoch": [2460310.5],
            "RA": [30.0],
            "DEC": [0.0],
            "band": ["V"],
            "observatory": ["500"],
            "note2": ["C"],
            "catalog": ["U"],
        }
    )


@pytest.mark.parametrize("module", ["tudatpy.data.mpc", "tudatpy.data.mpc.mpc"])
@pytest.mark.parametrize("source", ["pandas", "astropy", "file"])
def test_legacy_mpc_ingestion_weights_and_environment(
    legacy_mpc_table, monkeypatch, module, source
):
    BatchMPC = legacy_symbol(module, "BatchMPC")
    batch = BatchMPC()
    if source == "pandas":
        batch.from_pandas(legacy_mpc_table, custom_name="Eros")
    elif source == "astropy":
        batch.from_astropy(Table.from_pandas(legacy_mpc_table), custom_name="Eros")
    else:
        parsed = legacy_mpc_table.copy()
        parsed[["RA", "DEC"]] = np.deg2rad(parsed[["RA", "DEC"]])
        legacy = importlib.import_module("tudatpy.data.mpc._legacy")
        monkeypatch.setattr(legacy, "parse_80cols_file", lambda filename: Table.from_pandas(parsed))
        batch.from_file("fixture.obs", custom_name="Eros")
    batch.set_weights([7.0])
    bodies = earth_bodies()
    collection = batch.to_tudat(
        bodies,
        included_satellites=None,
        apply_weights_VFCC17=False,
        apply_star_catalog_debias=False,
    )
    assert isinstance(collection, ObservationCollection)
    np.testing.assert_allclose(
        np.array(collection.concatenated_observations).reshape(-1), [np.pi / 6, 0.0]
    )
    np.testing.assert_array_equal(np.array(collection.concatenated_weights).reshape(-1), [7.0, 7.0])
    assert bodies.does_body_exist("Eros")
    assert "500" in bodies.get("Earth").ground_station_list
    assert batch.bodies_created["Eros"] == "empty body"


def test_legacy_mpc_default_weighting_is_numerical(legacy_mpc_table):
    batch = legacy_symbol("tudatpy.data.mpc", "BatchMPC")()
    batch.from_pandas(legacy_mpc_table)
    bodies = earth_bodies()
    with pytest.warns(DeprecationWarning, match="get_weights_VFCC17"):
        collection = batch.to_tudat(bodies, None, apply_star_catalog_debias=False)
    np.testing.assert_allclose(collection.concatenated_weights, 1.0 / np.deg2rad(1.0 / 3600.0) ** 2)


def test_legacy_horizons_ephemeris_methods():
    Query = legacy_symbol("tudatpy.data.horizons", "HorizonsQuery")
    query = Query.__new__(Query)
    states = np.array([[0.0, 1, 2, 3, 4, 5, 6], [60.0, 7, 8, 9, 10, 11, 12]])
    query.cartesian = lambda **kwargs: states
    with pytest.warns(DeprecationWarning, match="create_ephemeris_tabulated"):
        ephemeris = query.create_ephemeris_tabulated("SSB", "J2000")
    assert isinstance(ephemeris, environment_setup.ephemeris.TabulatedEphemerisSettings)
    Batch = legacy_symbol("tudatpy.data.horizons.horizons", "HorizonsBatch")
    batch = Batch.__new__(Batch)
    batch._query_objects = {"433": SimpleNamespace(name="Eros", cartesian=query.cartesian)}
    settings = environment_setup.BodyListSettings("SSB", "J2000")
    with pytest.warns(DeprecationWarning, match="add_batch_ephemerides"):
        result = batch.add_batch_ephemerides(settings, "SSB", "J2000")
    assert result is None
    assert isinstance(
        settings.get("Eros").ephemeris_settings,
        environment_setup.ephemeris.TabulatedEphemerisSettings,
    )
    assert batch._names == ["Eros"]


def test_legacy_omm_tle_conversion_methods():
    Utils = legacy_symbol("tudatpy.data.spacetrack", "OMMUtils")
    first = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753"
    second = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"
    assert isinstance(Utils.tle_to_Tle_object(first, second), environment.Tle)
    assert isinstance(Utils.tle_to_TleEphemeris_object(first, second), environment.TleEphemeris)


@pytest.mark.parametrize("module", ["tudatpy.data", "tudatpy.data.processTrk234"])
def test_legacy_tnf_process_returns_collection(monkeypatch, module):
    Processor = legacy_symbol(module, "Trk234Processor")
    legacy = importlib.import_module("tudatpy.data.processTrk234._legacy")
    monkeypatch.setattr(
        legacy.trk234,
        "Reader",
        lambda path: SimpleNamespace(decode=lambda **kwargs: None, sfdu_list=[]),
    )
    processor = Processor(["fixture.tnf"], ["doppler"], "Probe")
    monkeypatch.setattr(processor.converters["doppler"], "extract", lambda records: tnf_records())
    collection = processor.process()
    assert isinstance(collection, ObservationCollection)
    np.testing.assert_array_equal(
        np.array(collection.concatenated_observations).reshape(-1), [123.0]
    )
    assert isinstance(Processor([], []).process(), ObservationCollection)


def test_legacy_tnf_environment_method_installs_ramps(monkeypatch):
    Processor = legacy_symbol("tudatpy.data.processTrk234.processor", "Trk234Processor")
    legacy = importlib.import_module("tudatpy.data.processTrk234._legacy")
    monkeypatch.setattr(
        legacy.trk234,
        "Reader",
        lambda path: SimpleNamespace(decode=lambda **kwargs: None, sfdu_list=[]),
    )
    start = pd.Timestamp("2000-01-01T12:00:00").to_pydatetime()
    end = pd.Timestamp("2000-01-01T12:01:00").to_pydatetime()
    ramps = pd.DataFrame(
        {
            "epoch": [start],
            "station": ["DSS-14"],
            "start_time": [start],
            "end_time": [end],
            "freq": [8.4e9],
            "rate": [0.1],
        }
    )
    processor = Processor(["fixture.tnf"], [])
    processor.ramp_converter = SimpleNamespace(
        extract=lambda records: ramps.copy(),
        process=lambda frame: frame,
        handle_open_ramps=lambda frame, handling: frame,
    )
    bodies = earth_bodies()
    environment_setup.add_ground_station(
        bodies.get("Earth"),
        environment_setup.ground_station.basic_station("DSS-14", [6378137.0, 0.0, 0.0]),
    )
    processor.set_tnf_information_in_bodies(bodies)
    calculator = bodies.get("Earth").get_ground_station("DSS-14").transmitting_frequency_calculator
    assert isinstance(calculator, environment.PiecewiseLinearFrequencyInterpolator)
    np.testing.assert_allclose(calculator.start_frequencies, [8.4e9])


def tnf_records():
    return pd.DataFrame(
        {
            "epoch": [pd.Timestamp("2024-01-01T12:00:00").to_pydatetime()],
            "link_ends": [("DSS-14", "123", "DSS-14")],
            "band": [("X", "X")],
            "link_delays": [(0.0, 0.0, 0.0)],
            "count_time": [60.0],
            "lowest_ranging_component": [1],
            "obs": [123.0],
        }
    )


@pytest.mark.parametrize("name", ["DerivedDopplerConverter", "DerivedSraRangeConverter"])
def test_legacy_tnf_converters_keep_return_types_and_station_times(name):
    from tudatpy.estimation.observations import SingleObservationSet

    Converter = legacy_symbol("tudatpy.data.processTrk234.converters", name)
    converter = Converter()
    records = tnf_records()
    result = converter.process(records, "Probe")
    assert len(result) == 1
    assert isinstance(result[0], SingleObservationSet)
    np.testing.assert_array_equal(
        np.array(result[0].concatenated_observations).reshape(-1), [123.0]
    )
    from tudatpy.astro import time_representation
    from tudatpy.dynamics.environment_setup.ground_station import (
        get_approximate_dsn_ground_station_positions,
    )

    expected_time = time_representation.default_time_scale_converter().convert_time(
        input_scale=time_representation.utc_scale,
        output_scale=time_representation.tdb_scale,
        input_value=time_representation.DateTime.from_python_datetime(
            records.iloc[0]["epoch"]
        ).to_epoch(),
        earth_fixed_position=get_approximate_dsn_ground_station_positions()["DSS-14"],
    )
    assert abs(float(result[0].observation_times[0]) - float(expected_time)) < 1.0e-9
    assert isinstance(converter.build_link_ends_dict(("DSS-14", "123", "DSS-14"), "Probe"), dict)
