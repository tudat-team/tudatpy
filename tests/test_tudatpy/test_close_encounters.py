"""B-plane coordinates checked against the GMAT unit test for B-plane parameters.

Reference: NASA GMAT, src/UnitTests/TestParam/TestBplaneParam.cpp, expected values credited
there to S. Hughes (2005).
"""

import numpy as np
import pytest

from tudatpy.astro import close_encounters

GM_EARTH = 398600.4415e9  # m^3/s^2, GMAT default

# Hyperbolic geocentric state, EarthMJ2000Eq [m, m/s]
RELATIVE_STATE = np.array(
    [
        233410.6846140172e3,
        83651.0868276347170e3,
        -168884.42195943173e3,
        -0.4038280708568842e3,
        2.0665425988121107e3,
        0.4654706868112324e3,
    ]
)

# GMAT expected values [m]
B_DOT_T = 365738.686341826e3
B_DOT_R = 276107.260600374e3


def test_b_plane_coordinates_match_gmat():
    b_plane = close_encounters.compute_b_plane_coordinates(RELATIVE_STATE, GM_EARTH)

    assert b_plane.b_t == pytest.approx(B_DOT_T, rel=1e-12)
    assert b_plane.b_r == pytest.approx(B_DOT_R, rel=1e-12)
