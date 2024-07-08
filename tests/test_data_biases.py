from tudatpy.data._biases import get_biases_EFCC18, load_bias_file, BIAS_LOWRES_FILE

import numpy as np
import pytest

from numpy.testing import assert_allclose


# 94% coverage of data/_biases.py
table_input = [
    (3988, "a", np.array([0.440, -0.322, -0.30, 0.53])),
    (0, "a", np.array([0.279, 0.326, -0.49, -1.04])),
]


bias_input = [
    # pixel: 3988    0.994838    +0.988506
    # RA DEC PMRA PMDEC: 0.440 -0.322  -0.30   0.53
    (
        1.0,
        1.0,
        0,
        "a",
        (0.440 / np.cos(1.0)) * ((2 * np.pi) / 1296000),
        (-0.322) * ((2 * np.pi) / 1296000),
    ),
    # pixel:       0    0.785398    +1.558038
    # RA DEC PMRA PMDEC: 0.279  0.326  -0.49  -1.04
    (
        0.785398,
        1.558038,
        0,
        "a",
        0.279 / np.cos(1.558038) * ((2 * np.pi) / 1296000),
        (0.326) * ((2 * np.pi) / 1296000),
    ),
    # 25 years after j2000
    # pixel:       0    0.785398    +1.558038
    # RA DEC PMRA PMDEC: 0.279  0.326  -0.49  -1.04
    (
        0.785398,
        1.558038,
        int(25 * 365.25 * 86400),
        "a",
        (0.279 + (25 * -0.49 / 1000)) / np.cos(1.558038) * ((2 * np.pi) / 1296000),
        (0.326 + (25 * -1.04 / 1000)) * ((2 * np.pi) / 1296000),
    ),
]


@pytest.mark.parametrize("pixel,catalog,expected", table_input)
def test_bias_table(pixel, catalog, expected):
    """test if the biases are retrieved for the correct pixels"""
    df, nside = load_bias_file(BIAS_LOWRES_FILE)

    assert nside == 64

    result = df.loc[(pixel, catalog), ["RA", "DEC", "PMRA", "PMDEC"]].to_numpy()

    assert_allclose(result, expected)


@pytest.mark.parametrize("RA,DEC,Time,Catalog,RA_cor,DEC_cor", bias_input)
def test_bias_values(RA, DEC, Time, Catalog, RA_cor, DEC_cor):
    """Check if values match hand calculation"""
    RA_cor_calc, DEC_cor_calc = get_biases_EFCC18(RA, DEC, Time, Catalog)

    assert (RA_cor_calc[0] - RA_cor) == pytest.approx(0.00)
    assert (DEC_cor_calc[0] - DEC_cor) == pytest.approx(0.00)
