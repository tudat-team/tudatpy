from tudatpy.data.mpc import BatchMPC
from tudatpy.data.horizons import HorizonsQuery
from tudatpy.astro import time_representation
from tudatpy.dynamics import environment_setup
from tudatpy.estimation import observable_models_setup, observations
from tudatpy.interface import spice
from astropy.table import Table
import numpy as np
import datetime
import json
import pandas as pd
import pytest
import warnings
from pathlib import Path
import tudatpy.data.mpc.mpc as mpc_module
from tudatpy.data.mpc.parser_80col import parse_80cols_data
from test_tudatpy.shared_constants import WGS84_EQUATORIAL_RADIUS, WGS84_FLATTENING

spice.load_standard_kernels()

# coverage = 88%
# TESTS DO NOT CHECK/VALIDATE:
# positions of observatories.

# Parameterised inputs
mpc_codes_test = [222, 999]
mpc_codes_test2 = [3]

get_observations_input = [
    ([999, 222], {"999", "222"}),
    ([222, "C/2012 S1"], {"222", "2012 S1"}),
]
get_observations_input2 = [
    (222, ValueError, "MPCcodes parameter must be list of integers/strings"),
    (
        [222, 1.0],
        ValueError,
        "All codes in the MPCcodes parameter must be integers or string",
    ),
]

filter_test_input = [
    (
        999,
        datetime.datetime(2022, 1, 1),
        datetime.datetime(2023, 1, 1),
        ["C51"],
        ["T08", "T05", "U55"],
        684,
        264,
        241,
        141,
    ),
    (
        222,
        datetime.datetime(2022, 1, 1),
        datetime.datetime(2023, 1, 1),
        ["C51"],
        ["T08", "T05", "U55"],
        575,
        214,
        209,
        7,
    ),
]


# for the weights tests
observatory_set_single = ["M22"]
observatory_set_multi = ["K19", "D67", "089", "706"]
weights_test_combinations = [
    (observatory_set_single, True),  # just one obs
    (observatory_set_single, False),
    (observatory_set_multi, False),
    (None, False),  # all data
]

APOPHIS_FIGURE2_RADAR_FIXTURE_PATH = (
    Path(__file__).resolve().parent / "data" / "mpc_radar_apophis_figure2_horizons_fixtures.json"
)
APOPHIS_FIGURE2_MAX_NORMALIZED_RADAR_RESIDUAL = 6.0


def _flattened_radec_and_times(observation_dataset):
    """Return flattened angular-position data in the same shape used by the MPC source table."""
    flattened_data = observation_dataset.ordered_flattened_observation_data()
    obs_radec = np.array(flattened_data.observation_vector).reshape(2, -1, order="F")
    obs_times = np.array(flattened_data.times).reshape(2, -1, order="F")
    return obs_radec, obs_times


# @pytest.mark.parametrize("inp,expected", get_observations_input)
# def test_BatchMPC_getobservations(inp, expected):
#    query = BatchMPC()
#    query.get_observations(inp)
#    assert set(query.MPC_objects) == expected


# @pytest.mark.parametrize("inp,errtype,errvalue", get_observations_input2)
# def test_BatchMPC_getobservations2(inp, errtype, errvalue):
#    query = BatchMPC()
#    with pytest.raises(Exception) as exc_info:
#        query.get_observations(inp)
#
#    assert exc_info.type is errtype
#    assert str(exc_info.value) == errvalue


@pytest.mark.parametrize("mpc_code", mpc_codes_test)
def test_create_observation_dataset_from_astropy_table(mpc_code):
    """Check if observatory table matches the primary ObservationDataset output."""
    query = BatchMPC()
    query.get_observations([mpc_code])
    query.filter(observatories=["T05", "T08"])
    # table values are sorted for easier comparison
    query._table = query._table.sort_values(["observatory", "epoch_seconds_TDB"])

    RADEC = query.table.loc[:, ["RA", "DEC"]].to_numpy().T
    times = query.table.loc[:, ["epoch_seconds_TDB"]].to_numpy().T[0]
    times = np.array([times, times])  # concat times are doubled due to RA + DEC

    # we created a table by using get_observations.
    # This yields observations in radians, so we have to set
    # in_degrees = False
    observation_dataset = query.create_observation_dataset_from_astropy_table(
        query._table, apply_weights_VFCC17=True, apply_star_catalog_debias=False, in_degrees=False
    )

    dataset_RADEC, dataset_times = _flattened_radec_and_times(observation_dataset)

    # Full-array comparisons catch ordering regressions that max/sum checks can hide.
    np.testing.assert_allclose(dataset_times, times)
    np.testing.assert_allclose(dataset_RADEC, RADEC)


@pytest.mark.parametrize("mpc_code", mpc_codes_test)
def test_mpc_custom_name_metadata(mpc_code):
    batch = BatchMPC()
    # Use a real MPC code but provide a custom name
    custom_name = "Death_Star"
    batch.get_observations([mpc_code], custom_name=custom_name)

    # Verify the custom_name column exists and has the correct value
    assert "custom_name" in batch.table.columns
    assert (batch.table["custom_name"] == custom_name).all()

    # Verify MPC_objects returns the custom name instead of mpc code
    assert custom_name in batch.MPC_objects
    assert mpc_code not in batch.MPC_objects


@pytest.mark.parametrize("mpc_code", mpc_codes_test)
def test_BatchMPC_to_tudat(mpc_code):
    """Check if observatory table matches the ObservationDataset returned by to_tudat."""
    query = BatchMPC()
    query.get_observations([mpc_code])
    query.filter(observatories=["T05", "T08"])

    # table values are sorted for easier comparison
    query._table = query._table.sort_values(["observatory", "epoch_seconds_TDB"])

    RADEC = query.table.loc[:, ["RA", "DEC"]].to_numpy().T
    times = query.table.loc[:, ["epoch_seconds_TDB"]].to_numpy().T[0]
    times = np.array([times, times])  # concat times are doubled due to RA + DEC

    # to_tudat needs a system of bodies with earth in it as input
    bodies_to_create = [
        "Earth",
    ]
    global_frame_origin = "SSB"
    global_frame_orientation = "J2000"
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation
    )
    bodies = environment_setup.create_system_of_bodies(body_settings)

    observation_dataset = query.to_tudat(
        bodies=bodies, included_satellites=None, apply_star_catalog_debias=False
    )

    dataset_RADEC, dataset_times = _flattened_radec_and_times(observation_dataset)

    # Full-array comparisons catch ordering regressions that max/sum checks can hide.
    np.testing.assert_allclose(dataset_times, times)
    np.testing.assert_allclose(dataset_RADEC, RADEC)


@pytest.mark.parametrize("mpc_code", mpc_codes_test)
def test_BatchMPC_to_tudat_with_satelite(mpc_code):
    """Check if space-telescope observations match the ObservationDataset returned by to_tudat."""
    query = BatchMPC()
    query.get_observations([mpc_code])
    query.filter(observatories=["C51"])

    # table values are sorted for easier comparison
    query._table = query._table.sort_values(["observatory", "epoch_seconds_TDB"])

    RADEC = query.table.loc[:, ["RA", "DEC"]].to_numpy().T
    times = query.table.loc[:, ["epoch_seconds_TDB"]].to_numpy().T[0]
    times = np.array([times, times])  # concat times are doubled due to RA + DEC

    # to_tudat needs a system of bodies with earth in it as input
    bodies_to_create = [
        "Earth",
    ]
    global_frame_origin = "SSB"
    global_frame_orientation = "J2000"
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation
    )
    bodies = environment_setup.create_system_of_bodies(body_settings)
    bodies.create_empty_body("Wise")

    observation_dataset = query.to_tudat(
        bodies=bodies, included_satellites={"C51": "Wise"}, apply_star_catalog_debias=False
    )

    dataset_RADEC, dataset_times = _flattened_radec_and_times(observation_dataset)

    # Full-array comparisons catch ordering regressions that max/sum checks can hide.
    np.testing.assert_allclose(dataset_times, times)
    np.testing.assert_allclose(dataset_RADEC, RADEC)


def test_BatchMPC_get_satellite_state_history():
    batch = BatchMPC()
    batch._MPC_space_telescopes = []
    batch._table = pd.DataFrame(
        {
            "number": ["433", "433", "433"],
            "epoch_seconds_TDB": [0.0, 10.0, 20.0],
            "epoch_seconds_UTC": [0.0, 10.0, 20.0],
            "epoch": [2451545.0, 2451545.0 + 10.0 / 86400.0, 2451545.0 + 20.0 / 86400.0],
            "RA": [0.0, 0.0, 0.0],
            "DEC": [0.0, 0.0, 0.0],
            "band": ["V", "V", "V"],
            "observatory": ["C51", "C51", "C51"],
            "spacecraft_position_x": [0.0, 10.0, 20.0],
            "spacecraft_position_y": [10.0, 10.0, 10.0],
            "spacecraft_position_z": [0.0, -10.0, -20.0],
        }
    )
    batch._refresh_metadata()

    state_history = batch.get_satellite_state_history("C51")
    epochs = sorted(state_history)
    states = np.array([state_history[epoch] for epoch in epochs])

    assert epochs == [0.0, 10.0, 20.0]
    np.testing.assert_allclose(
        states[:, :3], [[0.0, 10.0, 0.0], [10.0, 10.0, -10.0], [20.0, 10.0, -20.0]]
    )
    np.testing.assert_allclose(
        states[:, 3:], [[1.0, 0.0, -1.0], [1.0, 0.0, -1.0], [1.0, 0.0, -1.0]]
    )


def test_BatchMPC_get_observations_skips_raw_query_without_satellite_data(monkeypatch):
    calls = []

    def fake_get_observations(*args, **kwargs):
        calls.append(kwargs)
        if kwargs.get("get_mpcformat", False):
            raise AssertionError("Raw MPC query should be skipped for ground-based data.")
        return Table(
            {
                "number": [433],
                "epoch": [2451545.0],
                "RA": [0.0],
                "DEC": [0.0],
                "band": ["V"],
                "observatory": ["T08"],
                "note2": ["C"],
                "catalog": ["G"],
                "desig": [None],
            }
        )

    monkeypatch.setattr(mpc_module.MPC, "get_observations", fake_get_observations)
    batch = BatchMPC()
    batch._MPC_space_telescopes = ["C51"]

    batch.get_observations([433])

    assert len(calls) == 1
    assert "spacecraft_position_x" in batch.table


def test_BatchMPC_merges_split_raw_satellite_records(monkeypatch):
    raw_ground = "00433         C2021 06 07.42640918 08 15.401-41 22 02.35         12.0 V      500"
    raw_satellite = (
        "00433         S2021 06 07.42640918 08 15.401-41 22 02.35         12.0 V      500"
    )
    raw_parallax = (
        "00433         s2021 06 07.4264091 -198301.940 +198171.039 +56287.9850   ~6oMXC57"
    )

    def fake_get_observations(*args, **kwargs):
        return Table({"obs": [raw_ground, raw_satellite, raw_parallax]})

    monkeypatch.setattr(mpc_module.MPC, "get_observations", fake_get_observations)
    batch = BatchMPC()
    parsed_table = pd.DataFrame(
        {
            "observatory": ["500", "500"],
            "epoch_seconds_TDB": [0.0, 1.0],
            "RA": [0.0, 0.0],
            "DEC": [0.0, 0.0],
        }
    )

    table = batch._add_spacecraft_positions_from_mpc80(parsed_table, 433, None)

    assert np.isnan(table.loc[0, "spacecraft_position_x"])
    np.testing.assert_allclose(table.loc[1, "spacecraft_position_x"], -198301940.0)
    np.testing.assert_allclose(table.loc[1, "spacecraft_position_y"], 198171039.0)
    np.testing.assert_allclose(table.loc[1, "spacecraft_position_z"], 56287985.0)


def _horizons_cartesian_states(horizons_id, epochs):
    query = HorizonsQuery(
        query_id=horizons_id,
        location="500@SSB",
        epoch_list=list(sorted(set(float(epoch) for epoch in epochs))),
        extended_query=True,
    )
    return query.cartesian(frame_orientation="J2000")


def _local_target_ephemeris_settings(batch, horizons_id, observatory_code):
    speed_of_light = 299792458.0
    epochs = batch.table.sort_values("epoch_seconds_TDB", kind="stable")[
        "epoch_seconds_TDB"
    ].to_numpy(dtype=float)
    spacecraft_history = batch.get_satellite_state_history(observatory_code)
    spacecraft_positions = np.array([spacecraft_history[float(epoch)][:3] for epoch in epochs])
    earth_positions = np.array(
        [
            spice.get_body_cartesian_state_at_epoch(
                "Earth",
                "SSB",
                "J2000",
                "NONE",
                float(epoch),
            )[:3]
            for epoch in epochs
        ]
    )
    target_states_at_receive = _horizons_cartesian_states(horizons_id, epochs)[:, 1:7]
    receiver_positions = earth_positions + spacecraft_positions
    transmit_epochs = epochs - (
        np.linalg.norm(target_states_at_receive[:, :3] - receiver_positions, axis=1)
        / speed_of_light
    )

    sample_epochs = []
    for epoch_array in (epochs, transmit_epochs):
        for offset in (-600.0, -300.0, 0.0, 300.0, 600.0):
            sample_epochs.extend(epoch_array + offset)

    target_table = _horizons_cartesian_states(horizons_id, sample_epochs)
    return environment_setup.ephemeris.tabulated(
        {float(row[0]): row[1:7] for row in target_table},
        "SSB",
        "J2000",
    )


def _angular_position_observation_settings(observation_dataset):
    observation_settings = []
    links_done = []
    for set_metadata in observation_dataset.observation_set_metadata:
        if (
            set_metadata.observable_type
            != observable_models_setup.model_settings.angular_position_type
        ):
            continue
        link = observation_dataset.link_definition(set_metadata.link_definition_id)
        if link not in links_done:
            links_done.append(link)
            observation_settings.append(
                observable_models_setup.model_settings.angular_position(
                    link,
                    bias_settings=None,
                )
            )
    return observation_settings


def _radar_observation_settings(observation_dataset):
    observation_settings = []
    links_done = []
    light_time_corrections = [
        observable_models_setup.light_time_corrections.first_order_relativistic_light_time_correction(
            ["Sun"]
        )
    ]
    for set_metadata in observation_dataset.observation_set_metadata:
        link = observation_dataset.link_definition(set_metadata.link_definition_id)
        link_key = (set_metadata.observable_type, link)
        if link_key in links_done:
            continue
        links_done.append(link_key)
        if set_metadata.observable_type == observable_models_setup.model_settings.n_way_range_type:
            observation_settings.append(
                observable_models_setup.model_settings.n_way_range(
                    link,
                    light_time_corrections,
                    bias_settings=None,
                    time_scale_for_observable=time_representation.utc_scale,
                )
            )
        elif (
            set_metadata.observable_type
            == observable_models_setup.model_settings.doppler_measured_frequency_type
        ):
            observation_settings.append(
                observable_models_setup.model_settings.two_way_doppler_instantaneous_frequency(
                    link,
                    light_time_corrections,
                    bias_settings=None,
                )
            )
    return observation_settings


def _load_apophis_figure2_radar_fixture_cases():
    with APOPHIS_FIGURE2_RADAR_FIXTURE_PATH.open() as fixture_file:
        fixture = json.load(fixture_file)
    return fixture["cases"]


def _target_horizons_ephemeris_from_fixture(case):
    ephemeris_data = case["target_ephemeris"]
    return environment_setup.ephemeris.tabulated(
        {float(row[0]): np.asarray(row[1:7], dtype=float) for row in ephemeris_data["states"]},
        ephemeris_data["frame_origin"],
        ephemeris_data["frame_orientation"],
    )


def _compute_fixture_radar_prefit_residual(case):
    radar_table = parse_80cols_data(case["mpc_80col_records"])
    batch = BatchMPC()
    batch.from_astropy(radar_table, in_degrees=False)

    target_body = batch.MPC_objects[0]
    body_settings = environment_setup.get_default_body_settings(
        ["Sun", "Earth", "Moon"],
        "SSB",
        "J2000",
    )
    body_settings.get("Earth").shape_settings = environment_setup.shape.oblate_spherical(
        WGS84_EQUATORIAL_RADIUS,
        WGS84_FLATTENING,
    )
    body_settings.get("Earth").rotation_model_settings = (
        environment_setup.rotation_model.gcrs_to_itrs(
            environment_setup.rotation_model.iau_2006,
            "J2000",
        )
    )
    body_settings.add_empty_settings(target_body)
    body_settings.get(target_body).ephemeris_settings = _target_horizons_ephemeris_from_fixture(
        case
    )

    bodies = environment_setup.create_system_of_bodies(body_settings)
    observation_dataset = batch.to_tudat(
        bodies=bodies,
        included_satellites=None,
        apply_star_catalog_debias=False,
        apply_weights_VFCC17=False,
    )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        from tudatpy.numerical_simulation import estimation_setup

        observation_simulators = estimation_setup.create_observation_simulators(
            _radar_observation_settings(observation_dataset),
            bodies,
        )
    observations.compute_residuals_and_dependent_variables_for_dataset(
        observation_dataset,
        observation_simulators,
        bodies,
    )

    flattened = observation_dataset.ordered_flattened_observation_data()
    residual = np.asarray(flattened.residual_vector)[0]
    sigma = batch.radar_table["sigma"].to_numpy(dtype=float)[0]
    return residual, sigma, observation_dataset


def _compute_hst_prefit_residuals(target, horizons_id):
    hst_observatory_code = "250"
    batch = BatchMPC()
    batch.get_observations([target])
    batch.filter(
        epoch_start=datetime.datetime(1990, 1, 1),
        epoch_end=datetime.datetime(2026, 1, 1),
        observatories=[hst_observatory_code],
    )

    target_body = batch.MPC_objects[0]
    hst_body = f"HST_{target.replace(' ', '_')}"
    body_settings = environment_setup.get_default_body_settings(
        ["Sun", "Earth"],
        "SSB",
        "J2000",
    )
    body_settings.add_empty_settings(target_body)
    body_settings.get(target_body).ephemeris_settings = _local_target_ephemeris_settings(
        batch,
        horizons_id,
        hst_observatory_code,
    )
    body_settings.add_empty_settings(hst_body)
    body_settings.get(hst_body).ephemeris_settings = batch.get_satellite_ephemeris_settings(
        hst_observatory_code,
        frame_origin="Earth",
        frame_orientation="J2000",
    )

    bodies = environment_setup.create_system_of_bodies(body_settings)
    observation_dataset = batch.to_tudat(
        bodies=bodies,
        included_satellites={hst_observatory_code: hst_body},
        apply_star_catalog_debias=False,
        apply_weights_VFCC17=False,
    )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        from tudatpy.numerical_simulation import estimation_setup

        observation_simulators = estimation_setup.create_observation_simulators(
            _angular_position_observation_settings(observation_dataset),
            bodies,
        )
    observations.compute_residuals_and_dependent_variables_for_dataset(
        observation_dataset,
        observation_simulators,
        bodies,
    )
    flattened = observation_dataset.ordered_flattened_observation_data()
    return np.array(flattened.residual_vector).reshape(2, -1, order="F").T


def test_mpc_hst_space_astrometry_prefit_residuals():
    """Check HST MPC spacecraft astrometry against Tudat-computed residuals."""
    rad_to_arcsec = 180.0 / np.pi * 3600.0
    cases = [
        ("2003 BF91", "2003 BF91;"),
        ("2003 BH91", "2003 BH91;"),
        ("2020 KP11", "2020 KP11;"),
    ]

    for target, horizons_id in cases:
        residuals_arcsec = _compute_hst_prefit_residuals(target, horizons_id) * rad_to_arcsec
        assert np.max(np.abs(residuals_arcsec)) < 0.1, target
        assert np.max(np.sqrt(np.mean(residuals_arcsec**2, axis=0))) < 0.05, target


def test_mpc_apophis_figure2_radar_prefit_residuals_against_horizons():
    """Check Apophis Figure 2 radar residuals against frozen Horizons target states."""
    fixture_cases = _load_apophis_figure2_radar_fixture_cases()
    assert len(fixture_cases) == 50
    print(f"\nComparing {len(fixture_cases)} Apophis Figure 2 radar observations.")

    for case in fixture_cases:
        residual, sigma, _ = _compute_fixture_radar_prefit_residual(case)
        unit = "m" if case["jpl_radar_row"]["units"] == "us" else "Hz"
        normalized_residual = residual / sigma
        diagnostic = (
            f"{case['epoch_utc']}  {case['jpl_radar_row']['units']:>2s}  "
            f"station = {case['mpc_80col_records'][0][68:71]}->"
            f"{case['mpc_80col_records'][1][68:71]}  "
            f"residual = {residual:12.3f} {unit}  sigma = {sigma:10.3f} {unit}  "
            f"residual/sigma = {normalized_residual:10.3f}"
        )
        if abs(normalized_residual) >= APOPHIS_FIGURE2_MAX_NORMALIZED_RADAR_RESIDUAL:
            print(diagnostic)
        assert abs(normalized_residual) < APOPHIS_FIGURE2_MAX_NORMALIZED_RADAR_RESIDUAL, diagnostic


def test_compare_mpc_horizons_eph():
    """Compares true observations from BatchMPC to interpolated simulated RA/DEC from JPL Horizons"""
    batch = BatchMPC()
    batch.get_observations([433])

    # batch.filter takes python datetimes in UTC!
    batch.filter(
        epoch_start=datetime.datetime(2017, 1, 1),
        epoch_end=datetime.datetime(2022, 1, 1),
        observatories=["T08"],
    )

    # Horizons Query wants batch_times (or start_epoch, end_epoch) in UTC!!!
    batch_times = batch.table.epoch_seconds_UTC.to_numpy()
    eros = HorizonsQuery(
        query_id="433;", location="T08@399", epoch_list=batch_times, extended_query=True
    )

    # interpolated_observations returns times in TDB!!!
    radec_horizons = eros.interpolated_observations(degrees=False)

    # the retrieved batch.table has time columns: epoch [julian days in UTC], epoch_seconds_UTC [UTC datetime], epoch_seconds_TDB [TDB seconds]
    radec_mpc = batch.table.loc[:, ["epoch_seconds_TDB", "RA", "DEC"]].reset_index(drop=True)

    diff = (radec_horizons - radec_mpc).to_numpy()
    diff = np.abs(diff).max(axis=0)
    time_diff = diff[0]
    RA_diff = diff[1]
    DEC_diff = diff[2]

    assert time_diff < 1e-3
    assert RA_diff < 1e-5
    assert DEC_diff < 1e-5


# COMMENTED DUE TO REGULAR TIMEOUT AND FAILURE ON AZURE
# @pytest.mark.parametrize("mpc_code", mpc_codes_test2)
# def test_mpc_coverage(mpc_code):
#     batch_base = BatchMPC()
#     batch_base.get_observations([mpc_code])
#     batch_base.filter(
#         epoch_start=datetime.datetime(2021, 1, 1),
#         epoch_end=datetime.datetime(2022, 1, 1),
#     )
#
#     # properties
#     batch_base.table
#     batch_base.observatories
#     batch_base.space_telescopes
#     batch_base.MPC_objects
#     batch_base.size
#     batch_base.bands
#     batch_base.epoch_start
#     batch_base.epoch_end
#     len(batch_base)
#
#     # addition
#     batch2 = BatchMPC()
#     batch2.get_observations([1])
#     batch2.filter(
#         epoch_start=datetime.datetime(2021, 1, 1),
#         epoch_end=datetime.datetime(2022, 1, 1),
#     )
#     batch3 = batch_base + batch2
#
#     # copy
#     batch3copy = batch3.copy()
#
#     # from_pandas + from_astropy
#     batch4 = BatchMPC()
#     batch5 = BatchMPC()
#
#     batch4.from_astropy(astroquery_MPC.get_observations(mpc_code))
#     batch5.from_pandas(batch_base._table)  # type: ignore
#
#     # plotting
#     batch_base.plot_observations_temporal()
#     batch_base.plot_observations_sky()
#     batch_base.plot_observations_sky(projection="hammer")
#     batch_base.plot_observations_sky(projection="mollweide")
#     batch_base.plot_observations_sky(projection="lambert")
#
#     # obs_table
#     batch_base.observatories_table(only_in_batch=False)
#     batch_base.observatories_table(only_space_telescopes=True)
#     batch_base.observatories_table(exclude_space_telescopes=True)
#     batch_base.observatories_table(include_positions=True)
#
#     # summary
#     batch_base.summary()
