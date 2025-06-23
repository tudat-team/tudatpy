# %%
import pandas as pd
import requests
import os
import pytest
from tudatpy.interface import spice
from tudatpy.numerical_simulation.environment_setup import (
    get_default_body_settings,
    ground_station,
    create_system_of_bodies,
)
from tudatpy.numerical_simulation.estimation_setup import observation
from tudatpy.data.processTrk234.processor import Trk234Processor
from tudatpy.data.processTrk234 import converters as cnv


# -----------------------------------------------------------------------------
# Test function
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
    converter = cnv.RampConverter()
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
    converter = cnv.RampConverter()
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
    converter = cnv.RampConverter()
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
    converter = cnv.RampConverter()
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
    converter = cnv.RampConverter()
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


def test_reader():
    # Use the current directory as temporary path.
    tmp_dir = os.getcwd()
    local_filename = os.path.join(tmp_dir, "tnfp.dat")

    # Download the TNF file if not already present.
    url_tnf = "https://pds-geosciences.wustl.edu/radiosciencedocs/urn-nasa-pds-radiosci_documentation/dsn_trk-2-34/tnfp.dat"
    response = requests.get(url_tnf)
    assert response.status_code == 200, f"Failed to download TNF file from {url_tnf}"
    with open(local_filename, "wb") as f:
        f.write(response.content)

    # Create system of bodies.

    spice.load_standard_kernels()
    global_frame_origin = "SSB"
    global_frame_orientation = "J2000"
    body_settings = get_default_body_settings(
        ["Earth"], global_frame_origin, global_frame_orientation
    )
    # body_settings.get("Earth").shape_settings = shape.oblate_spherical_spice()
    body_settings.get("Earth").ground_station_settings = ground_station.dsn_stations()
    body_settings.add_empty_settings("-202")

    bodies = create_system_of_bodies(body_settings)

    # Create observation collection from the TNF file.
    trkProcessor = Trk234Processor(
        [local_filename],
        ["doppler"],
        spacecraft_name="-202",
    )
    observationCollection = trkProcessor.process()
    # trkProcessor.set_tnf_information_in_bodies(bodies) This requires tudatpy to be compiled with time scalr type tudat::Time

    single_obs_sets = observationCollection.get_single_observation_sets()
    assert single_obs_sets, "No observation sets found in the observation collection."

    obs_set = single_obs_sets[0]

    # Check doppler integration time.
    dopplerCount = obs_set.ancilliary_settings.get_float_settings(
        observation.doppler_integration_time
    )
    assert (
        dopplerCount == 1.0
    ), f"Expected doppler integration time 1.0, got {dopplerCount}"

    # Check link end delays.
    linkEndDelays = obs_set.ancilliary_settings.get_float_list_settings(
        observation.link_ends_delays
    )
    expected_delays = [4.915100149105456e-08, 0.0, -1.8370300836068054e-07]
    assert linkEndDelays == pytest.approx(
        expected_delays
    ), f"Expected link end delays {expected_delays}, got {linkEndDelays}"

    # Check link definition.
    linkEndType = obs_set.link_definition
    transmitter = linkEndType.link_end_id(observation.transmitter).reference_point
    sc = linkEndType.link_end_id(observation.reflector1).body_name
    rcv = linkEndType.link_end_id(observation.receiver).reference_point
    assert transmitter == "DSS-65", f"Expected transmitter 'DSS-65', got {transmitter}"
    assert sc == "-202", f"Expected spacecraft '-202', got {sc}"
    assert rcv == "DSS-65", f"Expected receiver 'DSS-65', got {rcv}"

    # Check observation times and values.
    obsTimes = obs_set.observation_times
    
    # Check observation times and values.
    obsTimes = obs_set.observation_times
    # This requires tudatpy to be compiled with time scalr type tudat::Time
    # assert obsTimes[0].to_float() == pytest.approx(
    #     617245672.6834568
    # ), f"Unexpected observation time: {obsTimes[0].to_float()}"
    assert obsTimes[0] == pytest.approx(
        617245672.6834568
    ), f"Unexpected observation time: {obsTimes[0]}"
    obsValues = obs_set.concatenated_observations
    assert obsValues[0] == pytest.approx(
        -8445646929.490659
    ), f"Unexpected observation value: {obsValues[0]}"

    os.remove(local_filename)


if __name__ == "__main__":
    pytest.main([__file__])
