# %%
import numpy as np
import pandas as pd
import requests
import os
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
from tudatpy.estimation.observations import create_observation_dataset_from_tracking_data
from tudatpy.data_input.tracking_data import TrackingData, tnf


# -----------------------------------------------------------------------------
# Ramp converter tests
# -----------------------------------------------------------------------------
def test_single_ramp():
    """A simple ramp with a clear start and end event should produce one merged interval."""
    data = [
        {
            "station": "A",
            "epoch": pd.Timestamp("2021-01-01 10:00:00"),
            "type": 1,
            "freq": 50.0,
            "rate": 0.0,
        },
        {
            "station": "A",
            "epoch": pd.Timestamp("2021-01-01 10:05:00"),
            "type": 4,
            "freq": 50.0,
            "rate": 0.0,
        },
    ]
    df = pd.DataFrame(data)
    converter = tnf._RampConverter()
    result = converter.process(df)
    expected = pd.DataFrame(
        [
            {
                "start_time": pd.Timestamp("2021-01-01 10:00:00"),
                "end_time": pd.Timestamp("2021-01-01 10:05:00"),
                "station": "A",
                "type": 1,
                "freq": 50.0,
                "rate": 0.0,
            }
        ]
    )
    pd.testing.assert_frame_equal(result.reset_index(drop=True), expected)


def test_continuation():
    """
    If a new start event is a continuation (its frequency matches the extrapolated value
    and the rate is equal), the previous ramp's end time is updated without creating a new record.
    """
    data = [
        {
            "station": "A",
            "epoch": pd.Timestamp("2021-01-01 10:00:00"),
            "type": 1,
            "freq": 50.0,
            "rate": 0.1,
        },
        {
            "station": "A",
            "epoch": pd.Timestamp("2021-01-01 10:05:00"),
            "type": 1,
            "freq": 80.0,
            "rate": 0.1,
        },
    ]
    df = pd.DataFrame(data)
    converter = tnf._RampConverter()
    result = converter.process(df)
    expected = pd.DataFrame(
        [
            {
                "start_time": pd.Timestamp("2021-01-01 10:00:00"),
                "end_time": pd.Timestamp("2021-01-01 10:05:00"),
                "station": "A",
                "type": 1,
                "freq": 50.0,
                "rate": 0.1,
            },
        ]
    )
    pd.testing.assert_frame_equal(result.reset_index(drop=True), expected)


def test_new_ramp_non_continuation():
    """
    When a new start event does not qualify as a continuation,
    the previous ramp's end time is finalized and a new ramp is started.
    """
    data = [
        {
            "station": "A",
            "epoch": pd.Timestamp("2021-01-01 10:00:00"),
            "type": 1,
            "freq": 50.0,
            "rate": 0.1,
        },
        # This start event does not match the expected continuation frequency.
        {
            "station": "A",
            "epoch": pd.Timestamp("2021-01-01 10:05:00"),
            "type": 1,
            "freq": 75.0,
            "rate": 0.1,
        },
        {
            "station": "A",
            "epoch": pd.Timestamp("2021-01-01 10:10:00"),
            "type": 4,
            "freq": 75.0,
            "rate": 0.1,
        },
    ]
    df = pd.DataFrame(data)
    converter = tnf._RampConverter()
    result = converter.process(df)
    expected = pd.DataFrame(
        [
            {
                "start_time": pd.Timestamp("2021-01-01 10:00:00"),
                "end_time": pd.Timestamp("2021-01-01 10:05:00"),
                "station": "A",
                "type": 1,
                "freq": 50.0,
                "rate": 0.1,
            },
            {
                "start_time": pd.Timestamp("2021-01-01 10:05:00"),
                "end_time": pd.Timestamp("2021-01-01 10:10:00"),
                "station": "A",
                "type": 1,
                "freq": 75.0,
                "rate": 0.1,
            },
        ]
    )
    pd.testing.assert_frame_equal(result.reset_index(drop=True), expected)


def test_end_without_open_ramp():
    """
    An end event (type 4 or 5) with no corresponding open ramp should be ignored,
    resulting in an empty merged DataFrame.
    """
    data = [
        {
            "station": "A",
            "epoch": pd.Timestamp("2021-01-01 10:00:00"),
            "type": 4,
            "freq": 50.0,
            "rate": 0.0,
        },
    ]
    df = pd.DataFrame(data)
    converter = tnf._RampConverter()
    result = converter.process(df)
    # Simply assert that the resulting DataFrame is empty.
    assert result.empty, "Expected an empty DataFrame when no open ramp exists."


def test_multiple_stations():
    """
    Test that ramps are merged independently per station.
    """
    data_A = [
        {
            "station": "A",
            "epoch": pd.Timestamp("2021-01-01 10:00:00"),
            "type": 1,
            "freq": 50.0,
            "rate": 0.0,
        },
        {
            "station": "A",
            "epoch": pd.Timestamp("2021-01-01 10:10:00"),
            "type": 4,
            "freq": 50.0,
            "rate": 0.0,
        },
    ]
    data_B = [
        {
            "station": "B",
            "epoch": pd.Timestamp("2021-01-01 11:00:00"),
            "type": 1,
            "freq": 60.0,
            "rate": 0.0,
        },
        {
            "station": "B",
            "epoch": pd.Timestamp("2021-01-01 11:05:00"),
            "type": 4,
            "freq": 60.0,
            "rate": 0.0,
        },
    ]
    df_A = pd.DataFrame(data_A)
    df_B = pd.DataFrame(data_B)

    all_ramps = pd.concat([df_A, df_B], ignore_index=True)
    all_ramps.sort_values("epoch", inplace=True)
    all_ramps.reset_index(drop=True, inplace=True)
    converter = tnf._RampConverter()
    result = converter.process(all_ramps)

    expected_A = pd.DataFrame(
        [
            {
                "start_time": pd.Timestamp("2021-01-01 10:00:00"),
                "end_time": pd.Timestamp("2021-01-01 10:10:00"),
                "station": "A",
                "type": 1,
                "freq": 50.0,
                "rate": 0.0,
            },
        ]
    )
    expected_B = pd.DataFrame(
        [
            {
                "start_time": pd.Timestamp("2021-01-01 11:00:00"),
                "end_time": pd.Timestamp("2021-01-01 11:05:00"),
                "station": "B",
                "type": 1,
                "freq": 60.0,
                "rate": 0.0,
            },
        ]
    )
    expected = pd.concat([expected_A, expected_B], ignore_index=True)
    # Sorting for consistency.
    result = result.sort_values(["station", "start_time"]).reset_index(drop=True)
    expected = expected.sort_values(["station", "start_time"]).reset_index(drop=True)
    pd.testing.assert_frame_equal(result, expected)


def test_open_final_ramp_left_unbounded():
    """A final ramp with no end event (no following start, no type 4/5) is returned by process() with end_time NaT."""
    data = [
        {
            "station": "A",
            "epoch": pd.Timestamp("2021-01-01 10:00:00"),
            "type": 1,
            "freq": 50.0,
            "rate": 0.0,
        },
    ]
    result = tnf._RampConverter().process(pd.DataFrame(data))
    assert len(result) == 1
    assert pd.isna(result["end_time"].iloc[0])


def test_handle_open_ramps_raise_exception():
    """raise_exception raises ValueError when an open ramp is present."""
    ramp_df = pd.DataFrame(
        [
            {
                "start_time": pd.Timestamp("2021-01-01 10:00:00"),
                "end_time": pd.NaT,
                "station": "A",
                "type": 1,
                "freq": 50.0,
                "rate": 0.0,
            }
        ]
    )
    with pytest.raises(ValueError):
        tnf._RampConverter().handle_open_ramps(ramp_df, tnf.OpenRampHandling.raise_exception)


def test_handle_open_ramps_raise_exception_no_open_ramps():
    """raise_exception does not raise when all ramps are already closed."""
    ramp_df = pd.DataFrame(
        [
            {
                "start_time": pd.Timestamp("2021-01-01 10:00:00"),
                "end_time": pd.Timestamp("2021-01-01 10:05:00"),
                "station": "A",
                "type": 1,
                "freq": 50.0,
                "rate": 0.0,
            }
        ]
    )
    result = tnf._RampConverter().handle_open_ramps(ramp_df, tnf.OpenRampHandling.raise_exception)
    assert result["end_time"].iloc[0] == pd.Timestamp("2021-01-01 10:05:00")


def test_handle_open_ramps_close_silently():
    """close_silently closes an open ramp with end_time = start_time + 1 s."""
    start = pd.Timestamp("2021-01-01 10:00:00")
    ramp_df = pd.DataFrame(
        [
            {
                "start_time": start,
                "end_time": pd.NaT,
                "station": "A",
                "type": 1,
                "freq": 50.0,
                "rate": 0.0,
            }
        ]
    )
    result = tnf._RampConverter().handle_open_ramps(ramp_df, tnf.OpenRampHandling.close_silently)
    assert result["end_time"].iloc[0] == start + pd.Timedelta(seconds=1)


def test_handle_open_ramps_close_silently_leaves_closed_untouched():
    """close_silently only modifies open ramps; closed intervals are unchanged."""
    closed_end = pd.Timestamp("2021-01-01 10:05:00")
    open_start = pd.Timestamp("2021-01-01 10:10:00")
    ramp_df = pd.DataFrame(
        [
            {
                "start_time": pd.Timestamp("2021-01-01 10:00:00"),
                "end_time": closed_end,
                "station": "A",
                "type": 1,
                "freq": 50.0,
                "rate": 0.0,
            },
            {
                "start_time": open_start,
                "end_time": pd.NaT,
                "station": "A",
                "type": 1,
                "freq": 50.0,
                "rate": 0.0,
            },
        ]
    )
    result = tnf._RampConverter().handle_open_ramps(ramp_df, tnf.OpenRampHandling.close_silently)
    assert result["end_time"].iloc[0] == closed_end
    assert result["end_time"].iloc[1] == open_start + pd.Timedelta(seconds=1)


# -----------------------------------------------------------------------------
# Pipeline tests: TrackingData -> create_observation_dataset_from_tracking_data
# -----------------------------------------------------------------------------
def _dsn_bodies():
    """System of bodies with DSN stations and a dummy spacecraft."""
    spice.load_standard_kernels()
    body_settings = get_default_body_settings(["Earth"], "SSB", "J2000")
    body_settings.get("Earth").ground_station_settings = ground_station.dsn_stations()
    body_settings.add_empty_settings("-202")
    return create_system_of_bodies(body_settings)


_DSN_LINK_ENDS = [
    (("Earth", "DSS-65"), "transmitter"),
    (("-202", ""), "reflector_1"),
    (("Earth", "DSS-65"), "receiver"),
]


def test_pipeline_doppler_synthetic():
    bodies = _dsn_bodies()

    tracking_data = TrackingData(
        "DsnNWayAveragedDoppler",
        _DSN_LINK_ENDS,
        np.array([[-8.4e9], [-8.4e9]]),
        [617245672.68, 617245673.68],
        "receiver",
        "UTC",
    )
    tracking_data.add_string_vector_ancillary_setting("frequency bands", ["X-band", "X-band"])
    tracking_data.add_string_ancillary_setting(
        "DSN reference frequency band at reception", "X-band"
    )
    tracking_data.add_double_ancillary_setting("DSN Doppler reference frequency", 0.0)
    tracking_data.add_double_ancillary_setting("Doppler observable integration time", 1.0)
    tracking_data.add_double_vector_ancillary_setting("link ends time delays", [0.0, 0.0, 0.0])

    observation_dataset = create_observation_dataset_from_tracking_data([tracking_data], bodies)
    assert observation_dataset.number_of_observation_sets == 1
    metadata = observation_dataset.get_observation_set_metadata(0)
    assert metadata.observable_type == ObservableType.dsn_n_way_averaged_doppler_type

    dataset_ancillary_settings = observation_dataset.ancillary_settings(
        metadata.ancillary_settings_id
    )
    doppler_count = dataset_ancillary_settings.get_float_settings(
        ancillary_settings.doppler_integration_time
    )
    assert doppler_count == 1.0

    link_def = observation_dataset.link_definition(metadata.link_definition_id)
    assert link_def.link_end_id(links.transmitter).reference_point == "DSS-65"
    assert link_def.link_end_id(links.reflector1).body_name == "-202"
    assert link_def.link_end_id(links.receiver).reference_point == "DSS-65"


def test_pipeline_range_synthetic():
    bodies = _dsn_bodies()

    tracking_data = TrackingData(
        "DsnNWayRange",
        _DSN_LINK_ENDS,
        np.array([[1000.0]]),
        [617245672.68],
        "receiver",
        "UTC",
    )
    tracking_data.add_string_vector_ancillary_setting("frequency bands", ["X-band", "X-band"])
    tracking_data.add_double_ancillary_setting("DSN sequential range lowest ranging component", 7.0)
    tracking_data.add_double_vector_ancillary_setting(
        "link ends time delays", [4.9151e-08, 0.0, -1.837e-07]
    )

    observation_dataset = create_observation_dataset_from_tracking_data([tracking_data], bodies)
    assert observation_dataset.number_of_observation_sets == 1
    metadata = observation_dataset.get_observation_set_metadata(0)
    assert metadata.observable_type == ObservableType.dsn_n_way_range_type

    dataset_ancillary_settings = observation_dataset.ancillary_settings(
        metadata.ancillary_settings_id
    )
    lrc = dataset_ancillary_settings.get_float_settings(
        ancillary_settings.sequential_range_lowest_ranging_component
    )
    assert lrc == 7.0

    expected_delays = [4.9151e-08, 0.0, -1.837e-07]
    link_end_delays = dataset_ancillary_settings.get_float_list_settings(
        ancillary_settings.link_ends_delays
    )
    assert link_end_delays == pytest.approx(expected_delays)

    link_def = observation_dataset.link_definition(metadata.link_definition_id)
    assert link_def.link_end_id(links.transmitter).reference_point == "DSS-65"
    assert link_def.link_end_id(links.reflector1).body_name == "-202"
    assert link_def.link_end_id(links.receiver).reference_point == "DSS-65"


# -----------------------------------------------------------------------------
# End-to-end reader test
# -----------------------------------------------------------------------------
@pytest.mark.remote_data
def test_reader(tmp_path):
    local_filename = str(tmp_path / "tnfp_tracking_data.dat")

    # Download the TNF file if not already present.
    url_tnf = "https://pds-geosciences.wustl.edu/radiosciencedocs/urn-nasa-pds-radiosci_documentation/dsn_trk-2-34/tnfp.dat"
    if not os.path.exists(local_filename) or os.path.getsize(local_filename) == 0:
        try:
            response = requests.get(url_tnf, timeout=30)
            response.raise_for_status()
        except requests.RequestException as error:
            pytest.skip(f"Required remote TNF test data is unavailable: {error}")
        with open(local_filename, "wb") as f:
            f.write(response.content)

    # Create system of bodies.

    spice.load_standard_kernels()
    global_frame_origin = "SSB"
    global_frame_orientation = "J2000"
    body_settings = get_default_body_settings(
        ["Earth"], global_frame_origin, global_frame_orientation
    )
    body_settings.get("Earth").ground_station_settings = ground_station.dsn_stations()
    body_settings.add_empty_settings("-202")

    bodies = create_system_of_bodies(body_settings)

    # Create tracking data from the TNF file.
    tracking_data, supplementary_data = tnf.read_tnf_data(
        [local_filename],
        ["doppler"],
        spacecraft_name="-202",
    )
    assert tracking_data, "No tracking data found."

    observation_dataset = create_observation_dataset_from_tracking_data(tracking_data, bodies)
    assert observation_dataset.number_of_observation_sets > 0, "No observation sets found."
    metadata = observation_dataset.get_observation_set_metadata(0)
    dataset_ancillary_settings = observation_dataset.ancillary_settings(
        metadata.ancillary_settings_id
    )

    # Check doppler integration time.
    dopplerCount = dataset_ancillary_settings.get_float_settings(
        ancillary_settings.doppler_integration_time
    )
    assert dopplerCount == 1.0, f"Expected doppler integration time 1.0, got {dopplerCount}"

    # Check link end delays.
    linkEndDelays = dataset_ancillary_settings.get_float_list_settings(
        ancillary_settings.link_ends_delays
    )
    expected_delays = [4.915100149105456e-08, 0.0, -1.8370300836068054e-07]
    assert linkEndDelays == pytest.approx(
        expected_delays
    ), f"Expected link end delays {expected_delays}, got {linkEndDelays}"

    # Check link definition.
    linkEndType = observation_dataset.link_definition(metadata.link_definition_id)
    transmitter = linkEndType.link_end_id(links.transmitter).reference_point
    sc = linkEndType.link_end_id(links.reflector1).body_name
    rcv = linkEndType.link_end_id(links.receiver).reference_point
    assert transmitter == "DSS-65", f"Expected transmitter 'DSS-65', got {transmitter}"
    assert sc == "-202", f"Expected spacecraft '-202', got {sc}"
    assert rcv == "DSS-65", f"Expected receiver 'DSS-65', got {rcv}"

    # Check observation times and values.
    obsTimes = observation_dataset.observation_times_for_set(0)

    # This requires tudatpy to be compiled with time scalar type tudat::Time
    # assert obsTimes[0].to_float() == pytest.approx(
    #     617245672.6834568
    # ), f"Unexpected observation time: {obsTimes[0].to_float()}"
    assert float(obsTimes[0]) == pytest.approx(
        617245672.6834568
    ), f"Unexpected observation time: {obsTimes[0]}"
    obsValues = observation_dataset.observation_vector_for_set(0)
    assert float(obsValues[0]) == pytest.approx(
        -8445646929.490659
    ), f"Unexpected observation value: {obsValues[0]}"


if __name__ == "__main__":
    pytest.main([__file__])
