from datetime import datetime
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from tudatpy.astro import time_representation as time
from tudatpy.data_input.tracking_data.gaia import GaiaAstrometry, load_gaia_astrometry
from tudatpy.dynamics import environment_setup
from tudatpy.estimation import observations

ARCHIVE = Path(__file__).with_name("gaia_astro_archive_for_tests.parquet")


def test_local_loading_utc_dates_matches_existing_filters():
    raw = GaiaAstrometry.load_from_local_archive(ARCHIVE, 673)
    start = datetime(2015, 1, 1)
    end = datetime(2016, 1, 1)
    converter = time.default_time_scale_converter()
    start_tdb = converter.convert_time(
        time.utc_scale, time.tdb_scale, time.DateTime.from_python_datetime(start).epoch()
    )
    end_tdb = converter.convert_time(
        time.utc_scale, time.tdb_scale, time.DateTime.from_python_datetime(end).epoch()
    )
    raw.apply_filters(epoch_start=start_tdb, epoch_end=end_tdb)
    loaded = load_gaia_astrometry(
        673, archive_path=ARCHIVE, epoch_start=start, epoch_end=end, time_scale=time.utc_scale
    )
    assert loaded is not None
    assert loaded.table.equals(raw.table)
    expected_tracking, _ = raw.to_tracking_data()
    actual_tracking, supplementary = loaded.to_tracking_data()
    for expected, actual in zip(expected_tracking, actual_tracking):
        np.testing.assert_array_equal(expected.observations, actual.observations)
    assert len(actual_tracking) == len(expected_tracking)
    bodies = environment_setup.create_system_of_bodies(
        environment_setup.get_default_body_settings([], "SSB", "J2000")
    )
    expected_dataset = observations.create_observation_dataset_from_tracking_data(
        expected_tracking, bodies
    )
    actual_dataset = observations.create_observation_dataset_from_tracking_data(
        actual_tracking, bodies
    )
    np.testing.assert_array_equal(
        expected_dataset.get_weight_matrix().toarray(), actual_dataset.get_weight_matrix().toarray()
    )
    assert supplementary[0].translational_state_supplementary_data.frame_origin == "Earth"


def test_local_loading_accepts_tdb_datetime_bounds():
    loaded = load_gaia_astrometry(
        673,
        archive_path=ARCHIVE,
        epoch_start=time.DateTime(2014, 1, 1),
        epoch_end=time.DateTime(2018, 1, 1),
    )
    assert loaded is not None and len(loaded.table) > 0


def test_empty_local_selection_is_allowed():
    assert (
        load_gaia_astrometry(
            673,
            archive_path=ARCHIVE,
            epoch_start=datetime(1900, 1, 1),
            epoch_end=datetime(1901, 1, 1),
            time_scale=time.utc_scale,
        )
        is None
    )


def test_missing_target_is_allowed():
    assert load_gaia_astrometry(999999, archive_path=ARCHIVE) is None


@pytest.mark.parametrize("message", ["Archive query failed", "No observations found for [673]"])
def test_query_failure_is_distinct_from_empty_data(message):
    with patch.object(GaiaAstrometry, "load_from_astroquery", side_effect=RuntimeError(message)):
        if message.startswith("No observations found"):
            assert load_gaia_astrometry(673) is None
        else:
            with pytest.raises(RuntimeError, match=message):
                load_gaia_astrometry(673)


def test_unreadable_archive_is_not_treated_as_empty():
    with patch.object(
        GaiaAstrometry, "load_from_local_archive", side_effect=OSError("archive unreadable")
    ):
        with pytest.raises(OSError, match="archive unreadable"):
            load_gaia_astrometry(673, archive_path=ARCHIVE)


def test_reversed_interval_fails_before_query():
    with patch.object(GaiaAstrometry, "load_from_astroquery") as query:
        with pytest.raises(ValueError, match="start"):
            load_gaia_astrometry(673, epoch_start=20.0, epoch_end=10.0)
        query.assert_not_called()
