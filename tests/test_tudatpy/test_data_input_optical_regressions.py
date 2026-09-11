import importlib

import numpy as np
import pandas as pd
import pytest
from astropy.table import Table

from tudatpy.dynamics import environment
from tudatpy.estimation.observations import create_observation_collection_from_tracking_data
from tudatpy.data_input.tracking_data.optical_utilities import read_optical_data

optical = importlib.import_module(
    "tudatpy.data_input.tracking_data.optical_utilities.optical_utilities"
)


@pytest.fixture
def optical_table():
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


def test_catalog_bias_is_subtracted_in_final_collection(optical_table, monkeypatch):
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
        data, environment.SystemOfBodies(), apply_corrections=True
    )
    np.testing.assert_allclose(
        np.array(collection.concatenated_observations).reshape(-1) - raw,
        -arcsec * np.array([1.0, 2.0]),
        rtol=1.0e-10,
        atol=1.0e-16,
    )


def weighting_bodies():
    from tudatpy.dynamics import environment_setup

    bodies = environment.SystemOfBodies()
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
    optical_table["note2"] = technique
    data, _ = read_optical_data(optical_table, add_weights=True, add_ancillary_data=ancillary)
    metadata = data[0].get_ancillary_settings_string_vector()
    assert metadata["note2"] == [technique]
    assert metadata["catalog"] == ["U"]
    collection = create_observation_collection_from_tracking_data(data, weighting_bodies())
    expected = 1.0 / np.deg2rad(sigma / 3600.0) ** 2
    np.testing.assert_allclose(collection.concatenated_weights, expected, rtol=1.0e-13)


@pytest.mark.parametrize("source", ["pandas", "astropy", "mpc"])
def test_custom_target_name_reaches_optical_link(optical_table, source):
    if source == "mpc":
        from tudatpy.data_input.tracking_data.mpc import BatchMPC

        batch = BatchMPC()
        batch._table = optical.create_augmented_optical_table(optical_table, custom_name="Eros")
        batch._refresh_metadata()
        data, _ = batch.to_tracking_dataset()
    else:
        table = Table.from_pandas(optical_table) if source == "astropy" else optical_table
        data, _ = read_optical_data(table, custom_name="Eros")
    assert (("Eros", ""), "transmitter") in data[0].link_ends
    assert data[0].get_ancillary_settings_string_vector()["number"] == ["433"]


def test_conflicting_optical_target_names_are_rejected(optical_table):
    table = pd.concat([optical_table, optical_table], ignore_index=True)
    table["custom_name"] = ["Eros", "Other"]
    table["observatory"] = ["500", "501"]
    with pytest.raises(ValueError, match="Conflicting custom names"):
        read_optical_data(table)
    table["custom_name"] = ["Eros", None]
    data, _ = read_optical_data(table)
    assert all((("Eros", ""), "transmitter") in item.link_ends for item in data)


def test_batch_metadata_resolves_partial_custom_names(optical_table):
    from tudatpy.data_input.tracking_data.mpc import BatchMPC

    table = pd.concat([optical_table, optical_table, optical_table], ignore_index=True)
    table["number"] = ["433", "433", "1"]
    table["custom_name"] = ["Eros", None, None]
    batch = BatchMPC()
    batch._table = optical.create_augmented_optical_table(table)
    batch._refresh_metadata()
    assert set(batch.MPC_objects) == {"Eros", "1"}
