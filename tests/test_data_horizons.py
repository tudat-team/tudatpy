from tudatpy.data.horizons import HorizonsQuery

from tudatpy.interface import spice
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
import datetime


spice.load_standard_kernels()

# 81% test coverage, Horizons Batch not tested


def test_compare_horizons_spice(pos_tolerance=1e3, vel_tolerance=1e-3):
    """Compares the states from an array of times between HorizonsQuery and SPICE"""
    start = 20 * 365 * 86400
    end = 23 * 365 * 86400

    npoints = 100

    times = np.linspace(start, end, npoints)

    query = HorizonsQuery(
        query_id="299",
        location="500@399",
        epoch_list=list(times),
        extended_query=True,
    )
    ref = "ECLIPJ2000"

    horizons_states = query.cartesian(frame_orientation=ref)[:, 1:]
    spice_states = []

    for time in times:
        # from spice
        state_spice = spice.get_body_cartesian_state_at_epoch(
            "Venus", "Earth", ref, "NONE", time
        )
        spice_states.append(state_spice)

    diff = np.abs(horizons_states - spice_states)
    max_diff_pos = diff.max(axis=0)[0:3].max()
    max_diff_vel = diff.max(axis=0)[3:6].max()

    assert max_diff_pos < pos_tolerance
    assert max_diff_vel < vel_tolerance

