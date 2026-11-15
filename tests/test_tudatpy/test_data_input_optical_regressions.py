"""Check names, catalogue corrections, and weights when optical observations are read."""

import importlib
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from astropy.table import Table

from tudatpy.dynamics import environment_setup
from tudatpy.estimation.observations import (
    create_observation_collection_from_tracking_data,
    create_observation_dataset_from_tracking_data,
)
from tudatpy.data_input.tracking_data.mpc import BatchMPC
from tudatpy.data_input.tracking_data.obs_80_cols import read_80_column_data
from tudatpy.data_input.tracking_data.optical_utilities import read_optical_data

optical = importlib.import_module(
    "tudatpy.data_input.tracking_data.optical_utilities.optical_utilities"
)
ades_reader = importlib.import_module("tudatpy.data_input.tracking_data.ades.reader")


@pytest.fixture
def optical_table():
    """Provide one asteroid position together with its catalogue and observing details."""
    return pd.DataFrame(
        {
            "number": ["433"],
            "observatory": ["500"],
            "epoch": [2460310.5],
            "RA": [30.0],
            "DEC": [0.0],
            "note2": ["C"],
            "catalog": ["U"],
            "band": ["V"],
            "mag": [18.0],
            "phottype": ["V"],
            "discovery": ["*"],
        }
    )


def empty_bodies():
    """Create the otherwise empty set of bodies needed to collect the observation."""
    return environment_setup.create_system_of_bodies(
        environment_setup.BodyListSettings("SSB", "J2000")
    )


def test_catalog_bias_is_subtracted_in_final_collection(optical_table, monkeypatch):
    """Apply known one- and two-arcsecond errors and verify both are subtracted."""
    arcsec = np.deg2rad(1.0 / 3600.0)
    # Supply the same known catalogue errors at every possible sky location.
    index = pd.MultiIndex.from_product([range(12), ["U", "unknown"]])
    bias_map = pd.DataFrame(0.0, index=index, columns=["RA", "DEC", "PMRA", "PMDEC"])
    bias_map.loc[(slice(None), "U"), "RA"] = 1.0
    bias_map.loc[(slice(None), "U"), "DEC"] = 2.0
    monkeypatch.setattr(optical, "load_bias_file", lambda **kwargs: (bias_map, 1))
    data, _ = read_optical_data(optical_table, add_star_catalog_corrections=True)
    raw = np.array(data[0].observations).reshape(-1)
    collection = create_observation_collection_from_tracking_data(
        data, empty_bodies(), apply_corrections=True
    )
    # The final two angles must equal the measured angles minus the supplied errors.
    np.testing.assert_allclose(
        np.array(collection.concatenated_observations).reshape(-1) - raw,
        -arcsec * np.array([1.0, 2.0]),
        rtol=1.0e-10,
        atol=1.0e-16,
    )


def test_adopted_biases_are_applied(optical_table):
    """Subtract adopted RA cos(Dec) and Dec biases instead of the star-catalog corrections, and use
    the negative adopted time bias."""
    arcsec = np.deg2rad(1.0 / 3600.0)
    index = pd.MultiIndex.from_product([range(12), ["U", "unknown"]])
    bias_map = pd.DataFrame(1.0, index=index, columns=["RA", "DEC", "PMRA", "PMDEC"])
    table = optical_table.assign(DEC=60.0, biasRA=0.5, biasDec=0.25)
    with patch.object(optical, "load_bias_file", lambda **kwargs: (bias_map, 1)):
        data, _ = read_optical_data(table, add_star_catalog_corrections=True)
    raw = np.array(data[0].observations).reshape(-1)
    collection = create_observation_collection_from_tracking_data(
        data, empty_bodies(), apply_corrections=True
    )
    np.testing.assert_allclose(
        np.array(collection.concatenated_observations).reshape(-1) - raw,
        -arcsec * np.array([0.5 / np.cos(np.deg2rad(60.0)), 0.25]),
        rtol=1.0e-10,
    )
    data, _ = read_optical_data(table.assign(biasTime=2.0))
    assert data[0].correction_settings.time_biases == [-2.0]


@pytest.mark.parametrize("source", ["file", "mpc"])
@pytest.mark.parametrize(
    "catalog, ra_bias, dec_bias",
    [("U", 1.0, 2.0), ("u", 3.0, 4.0), (" ", 0.0, 0.0), ("?", 0.0, 0.0)],
)
def test_mpc80_catalog_reaches_weights_and_corrections(
    tmp_path, monkeypatch, source, catalog, ra_bias, dec_bias
):
    """File and MPC-batch loading must apply the bias for the record's catalogue."""
    line = "00433         C2021 06 07.42640918 08 15.401-41 22 02.35         12.0 V      500"
    line = line[:71] + catalog + line[72:]
    index = pd.MultiIndex.from_product([range(12), ["U", "u", "unknown"]])
    bias_map = pd.DataFrame(0.0, index=index, columns=["RA", "DEC", "PMRA", "PMDEC"])
    bias_map.loc[(slice(None), "U"), ["RA", "DEC"]] = [1.0, 2.0]
    bias_map.loc[(slice(None), "u"), ["RA", "DEC"]] = [3.0, 4.0]
    monkeypatch.setattr(optical, "load_bias_file", lambda **kwargs: (bias_map, 1))

    if source == "file":
        path = tmp_path / "observations.txt"
        path.write_text(line + "\n")
        data, _ = read_80_column_data(
            [str(path)], add_weights=True, add_star_catalog_corrections=True
        )
    else:
        monkeypatch.setattr(
            BatchMPC, "_fetch_mpc80_records", staticmethod(lambda code, id_type: [line])
        )
        batch = BatchMPC()
        batch.get_observations(["433"], use_mpc80_format=True)
        data, _ = batch.to_tracking_dataset(add_weights=True, add_star_catalog_corrections=True)

    assert data[0].get_ancillary_settings_string_vector()["catalog"] == [catalog.strip()]
    raw = np.array(data[0].observations).reshape(-1)
    collection = create_observation_collection_from_tracking_data(
        data, weighting_bodies(), apply_corrections=True
    )
    np.testing.assert_allclose(
        np.array(collection.concatenated_observations).reshape(-1) - raw,
        -np.deg2rad(1.0 / 3600.0) * np.array([ra_bias / np.cos(raw[1]), dec_bias]),
        rtol=1.0e-10,
        atol=1.0e-16,
    )


@pytest.mark.parametrize("ancillary", [False, True])
def test_optical_metadata_does_not_create_simulation_settings(optical_table, ancillary):
    """Read with and without extra asteroid details; neither may alter the angle calculation."""
    data, _ = read_optical_data(optical_table, add_ancillary_data=ancillary)
    collection = create_observation_collection_from_tracking_data(data, empty_bodies())

    # The asteroid number remains available, but it creates no extra calculation instructions.
    assert data[0].get_ancillary_settings_string_vector()["number"] == ["433"]
    sets = collection.get_single_observation_sets()
    assert len(sets) == 1
    assert sets[0].ancillary_settings is None


def weighting_bodies():
    """Provide the Earth station needed to assign a weight to the measurement."""
    bodies = empty_bodies()
    bodies.create_empty_body("Earth")
    environment_setup.add_ground_station(
        bodies.get("Earth"),
        environment_setup.ground_station.basic_station("500", [6378137.0, 0.0, 0.0]),
    )
    return bodies


@pytest.mark.parametrize(
    "technique,catalog,sigma",
    [("C", "U", 1.0), ("P", "U", 2.5), (" ", "U", 2.5), ("C", " ", 1.5)],
)
def test_vfcc17_assigns_weight_without_requesting_extra_details(
    optical_table, technique, catalog, sigma
):
    """Without extra details, assign the weight for CCD and photographic data; blank codes count as unknown."""
    optical_table["note2"] = technique
    optical_table["catalog"] = catalog
    data, _ = read_optical_data(optical_table, vfcc17_weights=True, add_ancillary_data=False)
    # The observing method and catalogue must remain available to choose the weight.
    metadata = data[0].get_ancillary_settings_string_vector()
    assert metadata["note2"] == [technique]
    assert metadata["catalog"] == [catalog]
    # The final weight must match the published value for the selected observing method.
    collection = create_observation_collection_from_tracking_data(data, weighting_bodies())
    expected = 1.0 / np.deg2rad(sigma / 3600.0) ** 2
    np.testing.assert_allclose(collection.concatenated_weights, expected, rtol=1.0e-13)


@pytest.mark.parametrize("correlation", [0.4, np.nan])
@pytest.mark.parametrize("nightly_deweighting,observations_per_night", [(False, 1), (True, 3)])
def test_ades_weights_from_uncertainties_and_local_night(
    optical_table, nightly_deweighting, observations_per_night, correlation
):
    """Weight positions by their uncertainties and correlation, or by VFCC17 without; all three share a night, and only the
    ADES weights are divided by the number of observations in the night."""
    table = pd.concat([optical_table] * 3, ignore_index=True)
    # Station 500 is at zero longitude, so the local night is the integer part of the Julian date.
    table["epoch"] = [2460310.55, 2460310.65, 2460310.75]
    table["DEC"] = [0.0, 30.0, 60.0]
    table["rmsRA"] = [0.15, 0.2, np.nan]
    table["rmsDec"] = [0.13, 0.1, np.nan]
    table["rmsCorr"] = [correlation, np.nan, np.nan]
    data, _ = read_optical_data(
        table, ades_weights=True, vfcc17_weights=True, nightly_deweighting=nightly_deweighting
    )
    dataset = create_observation_dataset_from_tracking_data(data, weighting_bodies())

    # The weights are the inverse of the covariance of each observation, with the RA uncertainty scaled by 1/cos(Dec);
    # with nightly deweighting, weights are divided by N.
    arcsec = np.deg2rad(1.0 / 3600.0)
    expected = np.zeros((6, 6))
    for i in range(2):
        ra_sigma = table["rmsRA"][i] * arcsec / np.cos(np.deg2rad(table["DEC"][i]))
        dec_sigma = table["rmsDec"][i] * arcsec
        covariance = np.nan_to_num(table["rmsCorr"][i]) * ra_sigma * dec_sigma
        expected[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] = (
            np.linalg.inv([[ra_sigma**2, covariance], [covariance, dec_sigma**2]])
            / observations_per_night
        )
    # The last observation has no uncertainties, and gets the VFCC17 weight of CCD data with catalogue U (not reduced, as
    # there are fewer than five observations in the night).
    expected[4:, 4:] = np.eye(2) / arcsec**2
    np.testing.assert_allclose(dataset.get_weight_matrix().toarray(), expected, rtol=1.0e-12)


def ades_bodies():
    """Provide a rotating Earth with station 500, a fixed asteroid 433 and the Sun, for the rates of right ascension and
    declination and for the light deflection."""
    body_settings = environment_setup.BodyListSettings("SSB", "J2000")
    body_settings.add_empty_settings("Sun")
    body_settings.get("Sun").ephemeris_settings = environment_setup.ephemeris.constant(
        np.zeros(6), "SSB", "J2000"
    )
    body_settings.get("Sun").gravity_field_settings = environment_setup.gravity_field.central(
        1.32712440018e20
    )
    body_settings.add_empty_settings("Earth")
    body_settings.get("Earth").ephemeris_settings = environment_setup.ephemeris.constant(
        np.array([1.5e11, 0.0, 0.0, 0.0, 3.0e4, 0.0]), "SSB", "J2000"
    )
    body_settings.get("Earth").rotation_model_settings = environment_setup.rotation_model.simple(
        "J2000", "IAU_Earth", np.eye(3), 0.0, 7.292115e-5
    )
    body_settings.add_empty_settings("433")
    body_settings.get("433").ephemeris_settings = environment_setup.ephemeris.constant(
        np.array([3.0e11, 1.0e11, 2.0e10, 0.0, 0.0, 0.0]), "SSB", "J2000"
    )
    bodies = environment_setup.create_system_of_bodies(body_settings)
    environment_setup.add_ground_station(
        bodies.get("Earth"),
        environment_setup.ground_station.basic_station("500", [6378137.0, 0.0, 0.0]),
    )
    return bodies


def test_ades_minimum_uncertainty_time_uncertainties_and_adopted_biases():
    """Read four ADES observations at the same epoch (so with the same rates of right ascension and declination):
    0. CCD with rmsRA, rmsDec and rmsTime: ADES weights, with the timing uncertainty;
    1. CCD with rmsRA below minimum_ades_uncertainty and rmsTime: VFCC17 weights, with the timing uncertainty;
    2. stellar occultation with rmsRA and rmsDec below minimum_ades_uncertainty: ADES weights (not removed);
    3. only residual fields: weights from sigRA and sigDec, and the adopted biases replace the light deflection.
    """
    arcsec = np.deg2rad(1.0 / 3600.0)
    declination = 0.2
    parsed_table = Table(
        {
            "number": ["433"] * 4,
            "observatory": ["500"] * 4,
            "epoch": [2460310.55] * 4,
            "RA": [0.5] * 4,
            "DEC": [declination] * 4,
            "mode": ["CCD", "CCD", "OCC", "CCD"],
            "note2": ["C"] * 4,
            "catalog": ["U"] * 4,
            "band": ["V"] * 4,
            "mag": [18.0] * 4,
            "rmsRA": [0.5, 0.05, 0.01, np.nan],
            "rmsDec": [0.4, 0.4, 0.01, np.nan],
            "rmsTime": [1.0, 2.0, np.nan, np.nan],
            "sigRA": [np.nan, np.nan, np.nan, 0.3],
            "sigDec": [np.nan, np.nan, np.nan, 0.3],
            "biasRA": [np.nan, np.nan, np.nan, 0.2],
            "biasDec": [np.nan, np.nan, np.nan, -0.1],
        }
    )

    def read(time_uncertainties, light_deflection_bodies):
        with patch.object(ades_reader, "parse_ades_file", lambda file_path: parsed_table.copy()):
            data, _ = ades_reader.read_ades_data(
                "unused.psv",
                ades_weights=True,
                vfcc17_weights=True,
                time_uncertainties=time_uncertainties,
                minimum_ades_uncertainty=0.1,
                light_deflection_bodies=light_deflection_bodies,
            )
        return data

    bodies = ades_bodies()
    data = read(time_uncertainties=True, light_deflection_bodies=["Sun"])
    data_without_time_uncertainties = read(time_uncertainties=False, light_deflection_bodies=None)

    # Positional weights: ADES (with RA uncertainty scaled by 1/cos(Dec)), VFCC17 for CCD data with catalogue U, ADES
    # for the occultation, and the weights from the residual fields
    covariances = []
    for tracking_data in [data_without_time_uncertainties, data]:
        weights = (
            create_observation_dataset_from_tracking_data(tracking_data, bodies)
            .get_weight_matrix()
            .toarray()
        )
        covariances.append(
            [np.linalg.inv(weights[2 * i : 2 * i + 2, 2 * i : 2 * i + 2]) for i in range(4)]
        )
    expected_sigmas = [(0.5, 0.4), (1.0 * np.cos(declination), 1.0), (0.01, 0.01), (0.3, 0.3)]
    for covariance, (ra_cos_dec_sigma, dec_sigma) in zip(covariances[0], expected_sigmas):
        expected = np.diag(
            [(ra_cos_dec_sigma * arcsec / np.cos(declination)) ** 2, (dec_sigma * arcsec) ** 2]
        )
        np.testing.assert_allclose(covariance, expected, rtol=1.0e-12)

    # The timing uncertainties add sigma_T^2 r r^T, with the same rates r for all observations: rank one for the ADES
    # weights, four times as large for the VFCC17 weights (sigma_T = 2 s), and zero without rmsTime
    timing_covariances = [
        with_time - without for with_time, without in zip(covariances[1], covariances[0])
    ]
    scale = np.abs(timing_covariances[0]).max()
    assert scale > 0.0
    assert abs(np.linalg.det(timing_covariances[0])) < 1.0e-8 * scale**2
    np.testing.assert_allclose(
        timing_covariances[1], 4.0 * timing_covariances[0], rtol=1.0e-6, atol=1.0e-8 * scale
    )
    for timing_covariance in timing_covariances[2:]:
        np.testing.assert_allclose(timing_covariance, 0.0, atol=1.0e-8 * scale)

    # The light deflection is only applied to the observations without adopted biases, for which the adopted biases are
    # subtracted instead
    assert data[0].correction_settings.observation_correction_settings.observation_indices == [
        0,
        1,
        2,
    ]
    corrections = (
        np.array(
            create_observation_collection_from_tracking_data(
                data, bodies, apply_corrections=True
            ).concatenated_observations
        )
        - np.array(
            create_observation_collection_from_tracking_data(data, bodies).concatenated_observations
        )
    ).reshape(-1, 2)
    assert np.abs(corrections[0]).min() > 0.0
    np.testing.assert_allclose(corrections[1:3], corrections[[0, 0]], rtol=1.0e-12)
    np.testing.assert_allclose(
        corrections[3], -arcsec * np.array([0.2 / np.cos(declination), -0.1]), rtol=1.0e-12
    )


@pytest.mark.parametrize("source", ["pandas", "mpc"])
def test_custom_target_name_reaches_optical_link(optical_table, source):
    """Read asteroid 433 directly and through an MPC batch; both must name it Eros."""
    if source == "mpc":
        from tudatpy.data_input.tracking_data.mpc import BatchMPC

        batch = BatchMPC()
        batch._table = optical.create_augmented_optical_table(optical_table, custom_name="Eros")
        batch._refresh_metadata()
        data, _ = batch.to_tracking_dataset()
    else:
        data, _ = read_optical_data(optical_table, custom_name="Eros")
    # The calculation uses the name Eros while the original asteroid number remains available.
    assert (("Eros", ""), "transmitter") in data[0].link_ends
    assert data[0].get_ancillary_settings_string_vector()["number"] == ["433"]


def test_conflicting_optical_target_names_are_rejected(optical_table):
    """Reject two names for asteroid 433, but share its one supplied name between stations."""
    table = pd.concat([optical_table, optical_table], ignore_index=True)
    table["custom_name"] = ["Eros", "Other"]
    table["observatory"] = ["500", "501"]
    # Two different names for the same asteroid must be rejected, even at different stations.
    with pytest.raises(ValueError, match="Conflicting custom names"):
        read_optical_data(table)
    # When only one station supplies the name Eros, both stations must use it.
    table["custom_name"] = ["Eros", None]
    data, _ = read_optical_data(table)
    assert all((("Eros", ""), "transmitter") in item.link_ends for item in data)


def test_batch_metadata_preserves_target_order_and_resolves_names(optical_table):
    """Keep asteroid 433 before asteroid 1 and use one supplied Eros name for both 433 rows."""
    from tudatpy.data_input.tracking_data.mpc import BatchMPC

    table = pd.concat([optical_table, optical_table, optical_table], ignore_index=True)
    table["number"] = ["433", "433", "1"]
    table["custom_name"] = ["Eros", None, None]
    batch = BatchMPC()
    batch._table = optical.create_augmented_optical_table(table)
    batch._refresh_metadata()
    # The summary must keep first-seen order and share Eros's supplied name with its unnamed row.
    assert batch.MPC_objects == ["Eros", "1"]
