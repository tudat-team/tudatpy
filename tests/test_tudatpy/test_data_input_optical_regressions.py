import importlib

import numpy as np
import pandas as pd
import pytest
from astropy.table import Table

from tudatpy.dynamics import environment_setup
from tudatpy.estimation.observations import create_observation_collection_from_tracking_data
from tudatpy.data_input.tracking_data.optical_utilities import read_optical_data

optical = importlib.import_module(
    "tudatpy.data_input.tracking_data.optical_utilities.optical_utilities"
)


@pytest.fixture
def optical_table():
    """Provide one CCD observation with nonempty optical metadata."""
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
    """Create an empty environment using the supported factory."""
    return environment_setup.create_system_of_bodies(
        environment_setup.BodyListSettings("SSB", "J2000")
    )


def test_catalog_bias_is_subtracted_in_final_collection(optical_table, monkeypatch):
    """Subtract known catalog biases when the observation collection applies corrections."""
    arcsec = np.deg2rad(1.0 / 3600.0)
    # A deterministic all-sky map exercises the actual bias calculation too.
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
    # Compare the final measurements with raw minus the known RA/DEC biases.
    np.testing.assert_allclose(
        np.array(collection.concatenated_observations).reshape(-1) - raw,
        -arcsec * np.array([1.0, 2.0]),
        rtol=1.0e-10,
        atol=1.0e-16,
    )


@pytest.mark.parametrize("ancillary", [False, True])
def test_optical_metadata_does_not_create_simulation_settings(optical_table, ancillary):
    """Metadata-only optical input leaves no ancillary settings for the angular model to reject."""
    data, _ = read_optical_data(optical_table, add_ancillary_data=ancillary)
    collection = create_observation_collection_from_tracking_data(data, empty_bodies())

    # Both the mandatory target identifier and optional catalogue fields remain metadata.
    assert data[0].get_ancillary_settings_string_vector()["number"] == ["433"]
    sets = collection.get_single_observation_sets()
    assert len(sets) == 1
    assert sets[0].ancillary_settings is None


def weighting_bodies():
    """Provide the Earth station needed to calculate optical weights."""
    bodies = empty_bodies()
    bodies.create_empty_body("Earth")
    environment_setup.add_ground_station(
        bodies.get("Earth"),
        environment_setup.ground_station.basic_station("500", [6378137.0, 0.0, 0.0]),
    )
    return bodies


@pytest.mark.parametrize("ancillary", [False, True])
@pytest.mark.parametrize("technique,sigma", [("C", 1.0), ("P", 2.5)])
def test_vfcc17_uses_required_metadata_without_optional_ancillary(
    optical_table, ancillary, technique, sigma
):
    """Retain technique/catalog inputs and recover expected CCD and photographic weights."""
    optical_table["note2"] = technique
    data, _ = read_optical_data(optical_table, add_weights=True, add_ancillary_data=ancillary)
    # Required metadata survives independently of the optional ancillary switch.
    metadata = data[0].get_ancillary_settings_string_vector()
    assert metadata["note2"] == [technique]
    assert metadata["catalog"] == ["U"]
    # The actual collection must use the expected technique-dependent uncertainty.
    collection = create_observation_collection_from_tracking_data(data, weighting_bodies())
    expected = 1.0 / np.deg2rad(sigma / 3600.0) ** 2
    np.testing.assert_allclose(collection.concatenated_weights, expected, rtol=1.0e-13)


@pytest.mark.parametrize("source", ["pandas", "astropy", "mpc"])
def test_custom_target_name_reaches_optical_link(optical_table, source):
    """Propagate custom target names through pandas, astropy, and MPC ingestion."""
    if source == "mpc":
        from tudatpy.data_input.tracking_data.mpc import BatchMPC

        batch = BatchMPC()
        batch._table = optical.create_augmented_optical_table(optical_table, custom_name="Eros")
        batch._refresh_metadata()
        data, _ = batch.to_tracking_dataset()
    else:
        table = Table.from_pandas(optical_table) if source == "astropy" else optical_table
        data, _ = read_optical_data(table, custom_name="Eros")
    # Name the observation link with the custom label but retain the MPC identifier.
    assert (("Eros", ""), "transmitter") in data[0].link_ends
    assert data[0].get_ancillary_settings_string_vector()["number"] == ["433"]


def test_conflicting_optical_target_names_are_rejected(optical_table):
    """Reject conflicting names across stations and inherit a unique name for missing rows."""
    table = pd.concat([optical_table, optical_table], ignore_index=True)
    table["custom_name"] = ["Eros", "Other"]
    table["observatory"] = ["500", "501"]
    # Distinct names for one identifier are ambiguous even at different stations.
    with pytest.raises(ValueError, match="Conflicting custom names"):
        read_optical_data(table)
    # A missing name is unambiguous and inherits the identifier's supplied name.
    table["custom_name"] = ["Eros", None]
    data, _ = read_optical_data(table)
    assert all((("Eros", ""), "transmitter") in item.link_ends for item in data)


@pytest.mark.parametrize(
    "custom_names,expected_names",
    [([None, None, None], ["433", "1"]), (["Eros", None, None], ["Eros", "1"])],
)
def test_batch_metadata_preserves_target_order_and_resolves_names(
    optical_table, custom_names, expected_names
):
    """Preserve first-appearance target order with absent or partially supplied custom names."""
    from tudatpy.data_input.tracking_data.mpc import BatchMPC

    table = pd.concat([optical_table, optical_table, optical_table], ignore_index=True)
    table["number"] = ["433", "433", "1"]
    table["custom_name"] = custom_names
    batch = BatchMPC()
    batch._table = optical.create_augmented_optical_table(table)
    batch._refresh_metadata()
    # Refreshing metadata must neither sort identifiers nor lose resolved names.
    assert batch.MPC_objects == expected_names
