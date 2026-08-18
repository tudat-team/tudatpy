import numpy as np
import pandas as pd
import pytest
from tudatpy.data_input.environment_data import spice
from tudatpy.dynamics.environment_setup import (
    get_default_body_settings,
    ground_station,
    create_system_of_bodies,
)
from tudatpy.estimation.observations_setup import ancillary_settings
from tudatpy.estimation.observable_models_setup import links
from tudatpy.estimation.observable_models_setup.model_settings import ObservableType
from tudatpy.estimation.observations import create_observation_collection_from_tracking_data
from tudatpy.data_input.tracking_data.atdf._converters import (
    AtdfNwayDopplerConverter,
    AtdfRampConverter,
    AtdfNwayRangeConverter,
)

_SPACECRAFT_NAME = "-202"


# -----------------------------------------------------------------------------
# AtdfConverter helper tests
# -----------------------------------------------------------------------------
def test_atdf_station_to_tudat():
    converter = AtdfNwayRangeConverter()
    assert converter.atdf_station_to_tudat(" DSS 65 ") == "DSS-65"


def test_get_link_ends_two_way():
    converter = AtdfNwayRangeConverter()
    row = pd.Series({"Xmtr": "DSS 65", "Rcvr": "DSS 65"})
    assert converter.get_link_ends(row, _SPACECRAFT_NAME) == (
        "DSS-65",
        _SPACECRAFT_NAME,
        "DSS-65",
    )


def test_get_link_ends_one_way():
    converter = AtdfNwayRangeConverter()
    row = pd.Series({"Xmtr": "S/C", "Rcvr": "DSS 65"})
    assert converter.get_link_ends(row, _SPACECRAFT_NAME) == (
        "",
        _SPACECRAFT_NAME,
        "DSS-65",
    )


def test_get_link_end_delays_two_way():
    converter = AtdfNwayRangeConverter()
    row = pd.Series(
        {"Xmtr": "DSS 65", "XmtrDly (nsec)": 100.0, "ScDly (nsec)": 200.0, "RcvrDly (nsec)": 50.0}
    )
    assert converter.get_link_end_delays(row) == pytest.approx((1e-7, 2e-7, 5e-8))


def test_get_link_end_delays_one_way_folds_transmitter_into_spacecraft():
    """When the spacecraft is the transmitter, its delay is folded into the spacecraft delay."""
    converter = AtdfNwayRangeConverter()
    row = pd.Series(
        {"Xmtr": "S/C", "XmtrDly (nsec)": 100.0, "ScDly (nsec)": 200.0, "RcvrDly (nsec)": 50.0}
    )
    assert converter.get_link_end_delays(row) == pytest.approx((0.0, 3e-7, 5e-8))


def test_build_link_ends_two_way():
    converter = AtdfNwayRangeConverter()
    link_ends = converter.build_link_ends(("DSS-65", _SPACECRAFT_NAME, "DSS-65"), "Earth")
    assert link_ends == [
        (("Earth", "DSS-65"), "transmitter"),
        ((_SPACECRAFT_NAME, ""), "reflector_1"),
        (("Earth", "DSS-65"), "receiver"),
    ]


def test_build_link_ends_one_way():
    converter = AtdfNwayRangeConverter()
    link_ends = converter.build_link_ends(("", _SPACECRAFT_NAME, "DSS-65"), "Earth")
    assert link_ends == [
        ((_SPACECRAFT_NAME, ""), "transmitter"),
        (("Earth", "DSS-65"), "receiver"),
    ]


# -----------------------------------------------------------------------------
# Converter extract/process tests
# -----------------------------------------------------------------------------
def _msr_row(**overrides):
    row = {
        "Data Type": "2-Way-Range",
        "time_tag (DT)": pd.Timestamp("2021-01-01 10:00:00"),
        "Xmtr": "DSS 65",
        "Rcvr": "DSS 65",
        "UL": "X",
        "DL": "X",
        "XmtrDly (nsec)": 100.0,
        "ScDly (nsec)": 200.0,
        "RcvrDly (nsec)": 50.0,
        "Rng-LC": 1.0,
        "CT (sec)": 1.0,
        "Ref-Freq (Hz)": 22.0,
        "Ex": "X",
        "Observed": 1000.0,
    }
    row.update(overrides)
    return row


def test_range_converter_extract_and_process():
    df = pd.DataFrame([_msr_row()])
    converter = AtdfNwayRangeConverter()

    extracted = converter.extract(df, _SPACECRAFT_NAME)
    assert len(extracted) == 1

    tracking_data_list = converter.process(extracted)
    assert len(tracking_data_list) == 1

    tracking_data = tracking_data_list[0]
    assert tracking_data.observable_type == "DsnNWayRange"
    assert tracking_data.get_ancillary_settings_string_vector()["frequency bands"] == [
        "X-band",
        "X-band",
    ]
    assert tracking_data.get_ancillary_settings_double()[
        "DSN sequential range lowest ranging component"
    ] == pytest.approx(1.0)
    assert tracking_data.get_ancillary_settings_double_vector()[
        "link ends time delays"
    ] == pytest.approx([1e-7, 2e-7, 5e-8])


def test_doppler_converter_extract_and_process_groups_matching_observations():
    """Observations sharing link ends, band, ref. frequency, delays, count time, and exciter
    band are merged into a single TrackingData object."""
    df = pd.DataFrame(
        [
            _msr_row(
                **{
                    "Data Type": "3-Way-Doppler",
                    "time_tag (DT)": pd.Timestamp("2021-01-01 10:00:00"),
                    "Observed": -1000.0,
                }
            ),
            _msr_row(
                **{
                    "Data Type": "3-Way-Doppler",
                    "time_tag (DT)": pd.Timestamp("2021-01-01 10:01:00"),
                    "Observed": -1005.0,
                }
            ),
        ]
    )
    converter = AtdfNwayDopplerConverter()

    extracted = converter.extract(df, _SPACECRAFT_NAME)
    assert len(extracted) == 2

    tracking_data_list = converter.process(extracted)
    assert len(tracking_data_list) == 1

    tracking_data = tracking_data_list[0]
    assert tracking_data.observable_type == "DsnNWayAveragedDoppler"
    assert np.array(tracking_data.observations).flatten() == pytest.approx([-1000.0, -1005.0])
    double_settings = tracking_data.get_ancillary_settings_double()
    assert double_settings["DSN Doppler reference frequency"] == pytest.approx(22.0)
    assert double_settings["Doppler observable integration time"] == pytest.approx(1.0)
    assert tracking_data.get_ancillary_settings_string_vector()["frequency bands"] == [
        "X-band",
        "X-band",
    ]


def test_ramp_converter_extract():
    df = pd.DataFrame(
        [
            {
                "Start-Time (DT)": pd.Timestamp("2021-01-01 10:00:00"),
                "End-Time (DT)": pd.Timestamp("2021-01-01 10:05:00"),
                "Station": " DSS 65 ",
                "Frequency (Hz)": 7.2e9,
                "Rate (Hz/sec)": 0.1,
            }
        ]
    )
    result = AtdfRampConverter().extract(df)
    assert result["station"].iloc[0] == "DSS-65"
    assert result["freq"].iloc[0] == pytest.approx(7.2e9)
    assert result["rate"].iloc[0] == pytest.approx(0.1)


# -----------------------------------------------------------------------------
# Pipeline tests: TrackingData -> create_observation_collection_from_tracking_data
# -----------------------------------------------------------------------------
def _dsn_bodies():
    """System of bodies with DSN stations and a dummy spacecraft."""
    spice.load_standard_kernels()
    body_settings = get_default_body_settings(["Earth"], "SSB", "J2000")
    body_settings.get("Earth").ground_station_settings = ground_station.dsn_stations()
    body_settings.add_empty_settings(_SPACECRAFT_NAME)
    return create_system_of_bodies(body_settings)


def test_pipeline_range_synthetic():
    bodies = _dsn_bodies()

    df = pd.DataFrame([_msr_row()])
    extracted = AtdfNwayRangeConverter().extract(df, _SPACECRAFT_NAME)
    tracking_data = AtdfNwayRangeConverter().process(extracted)

    observation_collection = create_observation_collection_from_tracking_data(tracking_data, bodies)
    obs_sets = observation_collection.get_single_observation_sets()

    assert len(obs_sets) == 1
    obs_set = obs_sets[0]
    assert obs_set.observable_type == ObservableType.dsn_n_way_range_type

    lrc = obs_set.ancillary_settings.get_float_settings(
        ancillary_settings.ObservationAncillarySimulationVariable.sequential_range_lowest_ranging_component
    )
    assert lrc == pytest.approx(1.0)

    link_def = obs_set.link_definition
    assert link_def.link_end_id(links.LinkEndType.transmitter).reference_point == "DSS-65"
    assert link_def.link_end_id(links.LinkEndType.reflector1).body_name == _SPACECRAFT_NAME
    assert link_def.link_end_id(links.LinkEndType.receiver).reference_point == "DSS-65"


def test_pipeline_doppler_synthetic():
    bodies = _dsn_bodies()

    df = pd.DataFrame(
        [
            _msr_row(**{"Data Type": "2-Way-Doppler", "Observed": -8.4e9}),
            _msr_row(
                **{
                    "Data Type": "2-Way-Doppler",
                    "time_tag (DT)": pd.Timestamp("2021-01-01 10:01:00"),
                    "Observed": -8.4e9,
                }
            ),
        ]
    )
    extracted = AtdfNwayDopplerConverter().extract(df, _SPACECRAFT_NAME)
    tracking_data = AtdfNwayDopplerConverter().process(extracted)

    observation_collection = create_observation_collection_from_tracking_data(tracking_data, bodies)
    obs_sets = observation_collection.get_single_observation_sets()

    assert len(obs_sets) == 1
    obs_set = obs_sets[0]
    assert obs_set.observable_type == ObservableType.dsn_n_way_averaged_doppler_type

    doppler_count = obs_set.ancillary_settings.get_float_settings(
        ancillary_settings.ObservationAncillarySimulationVariable.doppler_integration_time
    )
    assert doppler_count == pytest.approx(1.0)

    link_def = obs_set.link_definition
    assert link_def.link_end_id(links.LinkEndType.transmitter).reference_point == "DSS-65"
    assert link_def.link_end_id(links.LinkEndType.reflector1).body_name == _SPACECRAFT_NAME
    assert link_def.link_end_id(links.LinkEndType.receiver).reference_point == "DSS-65"


if __name__ == "__main__":
    pytest.main([__file__])
