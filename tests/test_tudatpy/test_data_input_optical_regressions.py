"""Check names, catalogue corrections, and weights when optical observations are read."""

import importlib

import numpy as np
import pandas as pd
import pytest

from tudatpy.dynamics import environment_setup
from tudatpy.estimation.observations import create_observation_collection_from_tracking_data
from tudatpy.data_input.tracking_data.optical_utilities import read_optical_data

optical = importlib.import_module(
    "tudatpy.data_input.tracking_data.optical_utilities.optical_utilities"
)


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


@pytest.mark.parametrize("technique,sigma", [("C", 1.0), ("P", 2.5)])
def test_vfcc17_assigns_weight_without_requesting_extra_details(
    optical_table, technique, sigma
):
    """Without extra details, assign catalogue U's weight for CCD and photographic data."""
    optical_table["note2"] = technique
    data, _ = read_optical_data(optical_table, add_weights=True, add_ancillary_data=False)
    # The observing method and catalogue must remain available to choose the weight.
    metadata = data[0].get_ancillary_settings_string_vector()
    assert metadata["note2"] == [technique]
    assert metadata["catalog"] == ["U"]
    # The final weight must match the published value for the selected observing method.
    collection = create_observation_collection_from_tracking_data(data, weighting_bodies())
    expected = 1.0 / np.deg2rad(sigma / 3600.0) ** 2
    np.testing.assert_allclose(collection.concatenated_weights, expected, rtol=1.0e-13)


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
