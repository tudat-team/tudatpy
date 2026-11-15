from datetime import datetime
from io import BytesIO
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd
import pytest
from astropy.table import Table
from astropy.io.votable import from_table, writeto
from astropy.io.votable.tree import Info

from tudatpy.astro import time_representation as time
from tudatpy.data_input.tracking_data.gaia import GaiaAstrometry, read_gaia_astrometry_data
from tudatpy.dynamics import environment_setup
from tudatpy.estimation import observations

ARCHIVE = Path(__file__).with_name("gaia_astro_archive_for_tests.parquet")


def aip_response(table, status="OK"):
    result = from_table(table)
    # AIP assigns a VOTable field ID that differs from the source_id column name.
    for field in result.get_first_table().fields:
        if field.name == "source_id":
            field.ID = "datalinkID"
    result.resources[0].infos.append(Info(name="QUERY_STATUS", value=status))
    buffer = BytesIO()
    writeto(result, buffer)
    return Mock(content=buffer.getvalue())


def transmitter_name(tracking_data):
    return dict((role, link) for link, role in tracking_data.link_ends)["transmitter"][0]


def assert_same_tracking_data(expected, actual):
    assert len(actual) == len(expected)
    for expected_data, actual_data in zip(expected, actual):
        np.testing.assert_array_equal(expected_data.observations, actual_data.observations)
        np.testing.assert_array_equal(expected_data.epochs, actual_data.epochs)
        assert transmitter_name(actual_data) == transmitter_name(expected_data)


def test_local_loading_utc_dates_matches_existing_filters():
    raw = GaiaAstrometry.from_parquet(673, ARCHIVE)
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
    actual_tracking, supplementary = read_gaia_astrometry_data(
        673, parquet_path=ARCHIVE, epoch_start=start, epoch_end=end, time_scale=time.utc_scale
    )
    expected_tracking, _ = raw.to_tracking_data(light_deflection_bodies=("Sun",))
    assert len(actual_tracking) > 0
    assert_same_tracking_data(expected_tracking, actual_tracking)
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
    tracking, supplementary = read_gaia_astrometry_data(
        673,
        parquet_path=ARCHIVE,
        epoch_start=time.DateTime(2014, 1, 1),
        epoch_end=time.DateTime(2018, 1, 1),
    )
    assert len(tracking) > 0 and len(supplementary) == 1


def test_default_and_custom_target_body_name():
    tracking, _ = read_gaia_astrometry_data(673, parquet_path=ARCHIVE)
    assert {transmitter_name(data) for data in tracking} == {"673"}
    tracking, _ = read_gaia_astrometry_data(673, target_body_name="Edda", parquet_path=ARCHIVE)
    assert {transmitter_name(data) for data in tracking} == {"Edda"}


def test_empty_local_selection_is_allowed():
    assert read_gaia_astrometry_data(
        673,
        parquet_path=ARCHIVE,
        epoch_start=datetime(1900, 1, 1),
        epoch_end=datetime(1901, 1, 1),
        time_scale=time.utc_scale,
    ) == ([], [])


def test_missing_target_is_allowed():
    assert read_gaia_astrometry_data(999999, parquet_path=ARCHIVE) == ([], [])


@pytest.mark.parametrize("message", ["Archive query failed", "No observations found for [673]"])
def test_query_failure_is_distinct_from_empty_data(message):
    with patch.object(GaiaAstrometry, "from_aip", side_effect=RuntimeError(message)):
        if message.startswith("No observations found"):
            assert read_gaia_astrometry_data(673) == ([], [])
        else:
            with pytest.raises(RuntimeError, match=message):
                read_gaia_astrometry_data(673)


def test_unreadable_archive_is_not_treated_as_empty():
    with patch.object(GaiaAstrometry, "from_parquet", side_effect=OSError("archive unreadable")):
        with pytest.raises(OSError, match="archive unreadable"):
            read_gaia_astrometry_data(673, parquet_path=ARCHIVE)


def test_reversed_interval_fails_before_query():
    with patch.object(GaiaAstrometry, "from_aip") as query:
        with pytest.raises(ValueError, match="start"):
            read_gaia_astrometry_data(673, epoch_start=20.0, epoch_end=10.0)
        query.assert_not_called()


@pytest.mark.parametrize("target", [673, 779])
def test_online_response_is_cached_in_original_units_before_filtering(target, tmp_path):
    raw = pd.read_parquet(ARCHIVE, filters=[("number_mp", "==", target)])
    archive_path = tmp_path / f"gaia_{target}_fpr.parquet"
    source_ids = Table({"source_id": raw.source_id.unique()})
    responses = [aip_response(source_ids), aip_response(Table.from_pandas(raw))]
    with patch("requests.get", side_effect=responses) as query:
        # An empty selected interval must not discard observations from the cache.
        assert read_gaia_astrometry_data(
            target, parquet_path=archive_path, epoch_start=0, epoch_end=1
        ) == ([], [])
        assert query.call_count == 2
        assert "WHERE source_id IN" in query.call_args.kwargs["params"]["QUERY"]
        pd.testing.assert_frame_equal(pd.read_parquet(archive_path), raw.reset_index(drop=True))
        query.reset_mock()
        loaded, _ = read_gaia_astrometry_data(target, parquet_path=archive_path)
        query.assert_not_called()
    expected = GaiaAstrometry.from_parquet(target, ARCHIVE)
    expected.apply_filters()
    expected_tracking, _ = expected.to_tracking_data(light_deflection_bodies=("Sun",))
    assert_same_tracking_data(expected_tracking, loaded)


def test_empty_online_response_is_cached_to_avoid_repeated_queries(tmp_path):
    raw = pd.read_parquet(ARCHIVE).iloc[:0]
    archive_path = tmp_path / "gaia_999999_fpr.parquet"
    source_ids = Table({"source_id": np.array([], dtype=np.int64)})
    responses = [aip_response(source_ids), aip_response(Table.from_pandas(raw))]
    with patch("requests.get", side_effect=responses) as query:
        assert read_gaia_astrometry_data(999999, parquet_path=archive_path) == ([], [])
        assert query.call_count == 2
        assert "WHERE 1 = 0" in query.call_args.kwargs["params"]["QUERY"]
        assert archive_path.is_file() and pd.read_parquet(archive_path).empty
        query.reset_mock()
        assert read_gaia_astrometry_data(999999, parquet_path=archive_path) == ([], [])
        query.assert_not_called()


def test_failed_query_does_not_create_archive(tmp_path):
    archive_path = tmp_path / "gaia_673_fpr.parquet"
    with patch("requests.get", side_effect=RuntimeError("Archive query failed")):
        with pytest.raises(RuntimeError, match="Archive query failed"):
            read_gaia_astrometry_data(673, parquet_path=archive_path)
    assert not archive_path.exists()


def test_aip_queries_by_number_and_denomination(tmp_path):
    raw = pd.read_parquet(ARCHIVE)
    sources = Table({"source_id": raw.source_id.unique()})
    responses = [aip_response(sources), aip_response(Table.from_pandas(raw))]
    with patch("requests.get", side_effect=responses) as query:
        loaded = GaiaAstrometry.from_aip([673, "Nina"])
    source_query = query.call_args_list[0].kwargs["params"]["QUERY"]
    assert "number_mp IN (673) OR denomination IN ('nina')" in source_query
    assert set(loaded.table["number_mp"]) == {673, 779}


@pytest.mark.parametrize("status", ["ERROR", "OVERFLOW"])
def test_failed_or_truncated_aip_response_is_not_saved(status, tmp_path):
    archive_path = tmp_path / "gaia_673_fpr.parquet"
    raw = pd.read_parquet(ARCHIVE, filters=[("number_mp", "==", 673)])
    sources = Table({"source_id": raw.source_id.unique()})
    responses = [aip_response(sources), aip_response(Table.from_pandas(raw), status)]
    with patch("requests.get", side_effect=responses):
        with pytest.raises(RuntimeError, match=status):
            read_gaia_astrometry_data(673, parquet_path=archive_path)
    assert not archive_path.exists()
