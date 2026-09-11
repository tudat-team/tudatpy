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
