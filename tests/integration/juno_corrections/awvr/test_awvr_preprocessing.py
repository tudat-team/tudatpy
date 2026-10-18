"""Exercise AMC screening and station-event masking without mission data."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


@pytest.fixture
def helpers():
    path = Path(__file__).with_name("juno_utils.py")
    spec = importlib.util.spec_from_file_location("juno_awvr_integration_utils", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def record(**changes):
    fields = dict(
        utc_time=548445540.0,
        zenith_wet_delay=0.04261,
        zenith_dry_delay=2.06209,
        algorithm_id="010",
        cloud_flag=0,
        wvr_flag=0,
        mtp_flag=0,
        surface_flag=0,
        wind_flag=0,
        path_delay_flag=0,
    )
    fields.update(changes)
    return SimpleNamespace(**fields)


@pytest.mark.parametrize(
    "flag", ["cloud_flag", "wvr_flag", "mtp_flag", "surface_flag", "wind_flag", "path_delay_flag"]
)
@pytest.mark.parametrize("value", [1, 2, 3, 9])
def test_every_nonzero_quality_flag_is_rejected(helpers, flag, value):
    assert flag in helpers.awvr_rejection_reasons(record(**{flag: value}))


@pytest.mark.parametrize("component", ["wet", "dry"])
@pytest.mark.parametrize(
    "value, reason",
    [(np.nan, "nonfinite"), (np.inf, "nonfinite"), (-0.1, "negative"), (9.99, "missing")],
)
def test_invalid_and_missing_path_lengths_are_rejected(helpers, component, value, reason):
    reasons = helpers.awvr_rejection_reasons(record(**{f"zenith_{component}_delay": value}))
    assert f"{reason}_{component}" in reasons


def test_algorithm_metadata_and_zero_paths_are_preserved(helpers):
    assert helpers.awvr_rejection_reasons(record()) == []
    assert helpers.awvr_rejection_reasons(record(zenith_wet_delay=0.0)) == []
    assert "nonfinite_epoch" in helpers.awvr_rejection_reasons(record(utc_time=np.nan))


def test_single_unit_selection_preserves_quality_gaps(helpers, monkeypatch, tmp_path, capsys):
    unit01 = tmp_path / "GRV_JUGR_2017139_2017139_01_V03.AMC.txt"
    unit02 = tmp_path / "GRV_JUGR_2017139_2017139_02_V03.AMC.txt"
    unit01.touch()
    unit02.touch()
    read_paths = []

    def reader(path):
        read_paths.append(Path(path))
        return SimpleNamespace(
            records=[
                record(utc_time=1000.0),
                record(utc_time=1030.0, surface_flag=3),
                record(utc_time=1060.0),
            ]
        )

    monkeypatch.setattr(helpers, "read_dsn_amc_file", reader)
    delays = helpers.load_awvr_zenith_delays(tmp_path, "01")
    assert read_paths == [unit01]
    assert delays == {1000.0: [2.06209, 0.04261], 1060.0: [2.06209, 0.04261]}
    summary = capsys.readouterr().out
    assert "2/3 samples retained" in summary
    assert "010" in summary
    assert "surface_flag" in summary
    np.testing.assert_array_equal(helpers.awvr_segment_boundaries(delays), [1000.0, 1060.0])


def test_source_boundaries_include_endpoints_and_preserve_tsac_gap_interiors(helpers):
    boundaries = helpers.awvr_segment_boundaries(
        {
            1000.0: [2.0, 0.04],
            1030.0: [2.0, 0.1],
            1100.0: [2.0, 0.2],
            1130.0: [2.0, 0.3],
        }
    )
    np.testing.assert_array_equal(boundaries, [1000.0, 1030.0, 1100.0, 1130.0])
    assert helpers.crosses_awvr_transition(999.0, 1000.0, boundaries)
    assert helpers.crosses_awvr_transition(1029.0, 1031.0, boundaries)
    assert not helpers.crosses_awvr_transition(1040.0, 1090.0, boundaries)
    assert not helpers.crosses_awvr_transition(1001.0, 1029.0, boundaries)
    assert not helpers.crosses_awvr_transition(1200.0, 1260.0, boundaries)


def test_transmission_and_reception_events_both_control_the_mask(helpers, monkeypatch):
    # Fixed 600 s one-way light time and a synthetic 70 s TDB-UTC offset.
    # Receive epoch 2270 therefore has Earth transmission UTC=1000.
    speed = 299792458.0
    station = SimpleNamespace(
        station_state=SimpleNamespace(cartesian_position_at_reference_epoch=np.zeros(3))
    )
    earth = SimpleNamespace(
        get_ground_station=lambda name: station,
        rotation_model=SimpleNamespace(body_fixed_to_inertial_rotation=lambda t: np.eye(3)),
    )

    def state(body, *args):
        return np.array([600.0 * speed, 0, 0, 0, 0, 0]) if body == "JUNO" else np.zeros(6)

    monkeypatch.setattr(helpers, "constants", SimpleNamespace(SPEED_OF_LIGHT=speed))
    monkeypatch.setattr(helpers, "spice", SimpleNamespace(get_body_cartesian_state_at_epoch=state))
    monkeypatch.setattr(
        helpers,
        "time_representation",
        SimpleNamespace(
            default_time_scale_converter=lambda: SimpleNamespace(
                convert_time=lambda a, b, t, p: t - 70.0
            ),
            tdb_scale=1,
            utc_scale=2,
        ),
    )
    monkeypatch.setattr(
        helpers,
        "observable_models_setup",
        SimpleNamespace(
            links=SimpleNamespace(LinkEndType=SimpleNamespace(transmitter="tx", receiver="rx"))
        ),
    )
    monkeypatch.setattr(
        helpers,
        "estimation",
        SimpleNamespace(
            observations=SimpleNamespace(
                create_single_observation_set=lambda *args: args,
                ObservationCollection=lambda sets: sets,
            )
        ),
    )

    def observation_set(tx, rx, times):
        return SimpleNamespace(
            link_definition=SimpleNamespace(
                link_ends={
                    "tx": SimpleNamespace(reference_point=tx),
                    "rx": SimpleNamespace(reference_point=rx),
                }
            ),
            observation_times=times,
            list_of_observations=list(range(len(times))),
            observable_type="Doppler",
            reference_link_end="rx",
            ancillary_settings=None,
        )

    collection = SimpleNamespace(
        get_single_observation_sets=lambda: [
            observation_set("DSS-25", "DSS-43", [2270.0, 2390.0]),
            observation_set("DSS-43", "DSS-25", [1085.0, 1300.0]),
        ]
    )
    result = helpers.exclude_awvr_source_transition_counts(
        collection,
        SimpleNamespace(get=lambda name: earth),
        {1000.0: [2.0, 0.04], 1030.0: [2.0, 0.1]},
        [(0.0, 10000.0)],
    )
    assert [item[3] for item in result] == [[2390.0], [1300.0]]
