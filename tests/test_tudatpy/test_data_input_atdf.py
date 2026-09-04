import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tudatpy.astro import time_representation
from tudatpy.data_input import resource_paths
from tudatpy.data_input.environment_data import spice
from tudatpy.data_input.tracking_data.atdf import _processor as processor_module
from tudatpy.data_input.tracking_data.atdf import read_atdf_data
from tudatpy.data_input.tracking_data.atdf._converters import (
    AtdfNwayDopplerConverter,
    AtdfNwayRangeConverter,
    AtdfRampConverter,
)
from tudatpy.data_input.tracking_data.atdf._processor import AtdfTrackingDataProcessor
from tudatpy.dynamics.environment_setup import (
    create_system_of_bodies,
    get_default_body_settings,
    ground_station,
)
from tudatpy.estimation.observable_models_setup import links
from tudatpy.estimation.observable_models_setup.model_settings import ObservableType
from tudatpy.estimation.observations import create_observation_collection_from_tracking_data
from tudatpy.estimation.observations_setup import ancillary_settings

SPACECRAFT = "MGS"
EARTH = "Earth"

# ``atdf2ascii`` output for PDS file 9066071a.tdf (Mars Global Surveyor, March
# 1999), taken from ``atdf2ascii/Examples``
FIXTURE_DIR = Path(resource_paths.get_test_data_path()) / "atdf"
FIXTURE_STEM = "9066071a"

# Row counts of the fixture's ``.msr`` table per "Data Type", and of its``.ramp`` table per station
MSR_ROWS_PER_DATA_TYPE = {
    "2-Way-Doppler": 525,
    "3-Way-Doppler": 6,
    "2-Way-Range": 138,
    "1-Way-Range": 20,
}
MSR_ROW_COUNT = sum(MSR_ROWS_PER_DATA_TYPE.values())
NWAY_DOPPLER_ROW_COUNT = (
    MSR_ROWS_PER_DATA_TYPE["2-Way-Doppler"] + MSR_ROWS_PER_DATA_TYPE["3-Way-Doppler"]
)
NWAY_RANGE_ROW_COUNT = MSR_ROWS_PER_DATA_TYPE["2-Way-Range"]
RAMP_ROWS_PER_STATION = {"DSS-34": 746, "DSS-45": 629, "DSS-15": 179, "DSS-54": 10}
RAMP_ROW_COUNT = sum(RAMP_ROWS_PER_STATION.values())

FIXTURE_COUNT_TIME = 60.0
FIXTURE_EXCITER_BAND = "S"
FIXTURE_LOWEST_RANGING_COMPONENT = 21

FIXTURE_DOPPLER_GROUP_SIZES = [5, 6, 6, 26, 56, 67, 67, 99, 199]

FIRST_DOPPLER_RECORD = {
    "epoch": pd.Timestamp("1999-03-07 19:27:35"),
    "observable": -19094.1917333329,
    "reference_frequency": 2114118912.0,
}
FIRST_RANGE_RECORD = {
    "epoch": pd.Timestamp("1999-03-07 19:30:05"),
    "observable": 28529144.580163553,
}


# -----------------------------------------------------------------------------
# Row builders
# -----------------------------------------------------------------------------
def msr_row(**overrides) -> dict:
    """A single processed ``.msr`` row, defaulting to a 2-way X-band Doppler point.

    Only the columns the converters actually read are included; ``overrides``
    lets each test state the one field whose effect it is checking.
    """
    row = {
        "Data Type": "2-Way-Doppler",
        "time_tag (DT)": pd.Timestamp("1999-03-07 19:27:35"),
        "Xmtr": "DSS 34",
        "Rcvr": "DSS 34",
        "UL": "X",
        "DL": "X",
        "Ex": "S",
        "CT (sec)": 60.0,
        "Rng-LC": 21.0,
        "Observed": -19094.19,
        "Ref-Freq (Hz)": 2114118912.0,
        "XmtrDly (nsec)": 0.0,
        "RcvrDly (nsec)": 0.0,
        "ScDly (nsec)": 0.0,
    }
    row.update(overrides)
    return row


def msr_frame(*rows: dict) -> pd.DataFrame:
    """Build a processed ``.msr`` DataFrame from :func:`msr_row` overrides."""
    return pd.DataFrame([msr_row(**row) for row in rows] or [msr_row()])


def ramp_row(**overrides) -> dict:
    """A single processed ``.ramp`` row."""
    row = {
        "Start-Time (DT)": pd.Timestamp("1999-03-07 11:46:54"),
        "End-Time (DT)": pd.Timestamp("1999-03-07 11:47:00"),
        "Station": " DSS 34",
        "Frequency (Hz)": 7164234321.75112,
        "Rate (Hz/sec)": 0.0,
    }
    row.update(overrides)
    return row


def ramp_frame(*rows: dict) -> pd.DataFrame:
    return pd.DataFrame([ramp_row(**row) for row in rows] or [ramp_row()])


def expected_epoch(timestamp: pd.Timestamp) -> float:
    """The UTC seconds since J2000 that the converters should produce for a time tag."""
    return time_representation.DateTime.from_python_datetime(timestamp.to_pydatetime()).to_epoch()


# -----------------------------------------------------------------------------
# Fixtures
# -----------------------------------------------------------------------------
@pytest.fixture
def make_processor():
    """Factory for a processor reading the MGS fixture, with overridable flags."""

    def _make(**overrides) -> AtdfTrackingDataProcessor:
        return AtdfTrackingDataProcessor(
            atdf_file_path=[Path(FIXTURE_STEM)], spacecraft_name=SPACECRAFT, **overrides
        )

    return _make


@pytest.fixture(scope="module")
def msr_df() -> pd.DataFrame:
    """The fixture's ``.msr`` table, parsed once. Treated as read-only by tests."""
    processor = AtdfTrackingDataProcessor([Path(FIXTURE_STEM)], SPACECRAFT)
    processor.read_atdf_ascii_msr(FIXTURE_DIR)
    return processor.df_processed_msr


@pytest.fixture(scope="module")
def ramp_df() -> pd.DataFrame:
    """The fixture's ``.ramp`` table, parsed once. Treated as read-only by tests."""
    processor = AtdfTrackingDataProcessor([Path(FIXTURE_STEM)], SPACECRAFT)
    processor.read_atdf_ascii_ramp(FIXTURE_DIR)
    return processor.df_processed_rmp


@pytest.fixture
def stub_atdf_to_ascii(monkeypatch):
    """Replace ``atdf2ascii.atdf_to_ascii`` with a recorder that emits the fixture.

    Decoding a binary ATDF file is ``atdf2ascii``'s responsibility, so the tests
    here only need the call to be made correctly and the expected ASCII tables to
    appear in ``output_dir``. The returned list collects the keyword arguments of
    every call, in order.
    """
    calls: list[dict] = []

    def _fake_atdf_to_ascii(*, input_file, output_dir, **kwargs):
        calls.append({"input_file": input_file, "output_dir": output_dir, **kwargs})
        for suffix in ("msr", "ramp"):
            shutil.copy(
                FIXTURE_DIR / f"{FIXTURE_STEM}.{suffix}",
                Path(output_dir) / f"{Path(input_file).stem}.{suffix}",
            )

    monkeypatch.setattr(processor_module, "atdf_to_ascii", _fake_atdf_to_ascii)
    return calls


@pytest.fixture(scope="module")
def dsn_bodies():
    """Earth with the DSN ground stations, plus an empty spacecraft body."""
    spice.load_standard_kernels()
    body_settings = get_default_body_settings([EARTH], "SSB", "J2000")
    body_settings.get(EARTH).ground_station_settings = ground_station.dsn_stations()
    body_settings.add_empty_settings(SPACECRAFT)
    return create_system_of_bodies(body_settings)


# -----------------------------------------------------------------------------
# Test helper functions shared by the converters
# -----------------------------------------------------------------------------
class TestConverterHelpers:
    """``AtdfConverter`` methods translating single ATDF fields to tudat conventions."""

    @pytest.fixture
    def converter(self) -> AtdfNwayRangeConverter:
        # The helpers under test live on the shared AtdfConverter base class; the
        # range converter is used simply as a concrete instantiation of it.
        return AtdfNwayRangeConverter()

    @pytest.mark.parametrize(
        "atdf_name, tudat_name",
        [("DSS 34", "DSS-34"), (" DSS 65 ", "DSS-65"), ("DSS-45", "DSS-45")],
    )
    def test_station_names_are_stripped_and_hyphenated(self, converter, atdf_name, tudat_name):
        assert converter.atdf_station_to_tudat(atdf_name) == tudat_name

    def test_two_way_link_ends(self, converter):
        row = pd.Series(msr_row(Xmtr="DSS 34", Rcvr="DSS 34"))
        assert converter.get_link_ends(row, SPACECRAFT) == ("DSS-34", SPACECRAFT, "DSS-34")

    def test_three_way_link_ends(self, converter):
        row = pd.Series(msr_row(Xmtr="DSS 45", Rcvr="DSS 54"))
        assert converter.get_link_ends(row, SPACECRAFT) == ("DSS-45", SPACECRAFT, "DSS-54")

    def test_link_end_delays_are_converted_from_nanoseconds(self, converter):
        row = pd.Series(
            msr_row(**{"XmtrDly (nsec)": 100.0, "ScDly (nsec)": 200.0, "RcvrDly (nsec)": 50.0})
        )
        assert converter.get_link_end_delays(row) == pytest.approx((1e-7, 2e-7, 5e-8))

    def test_link_ends_structure_for_nway_obs(self, converter):
        assert converter.build_link_ends(("DSS-45", SPACECRAFT, "DSS-54"), EARTH) == [
            ((EARTH, "DSS-45"), "transmitter"),
            ((SPACECRAFT, ""), "reflector_1"),
            ((EARTH, "DSS-54"), "receiver"),
        ]

    def test_epoch_conversion_yields_utc_seconds_since_j2000(self, converter):
        timestamp = pd.Timestamp("1999-03-07 19:27:35")
        assert converter.to_utc_epoch(timestamp.to_pydatetime()).to_float() == pytest.approx(
            -25893145.0, rel=0, abs=1e-6
        )

    @pytest.mark.parametrize(
        "band_code, band_name", [("S", "S-band"), ("X", "X-band"), ("Ka", "Ka-band")]
    )
    def test_supported_frequency_bands_map_to_tudat_names(self, converter, band_code, band_name):
        assert converter.frequency_bands_mapping[band_code] == band_name


# -----------------------------------------------------------------------------
# ASCII table reading
# -----------------------------------------------------------------------------
class TestAsciiTableReading:
    """Parsing of the intermediate ``.msr``/``.ramp`` tables written by ``atdf2ascii``."""

    def test_msr_columns_come_from_the_table_header(self, msr_df):
        assert set(msr_df.columns) >= {
            "time_tag (UTC)",
            "Data Type",
            "Xmtr",
            "Rcvr",
            "UL",
            "DL",
            "Ex",
            "CT (sec)",
            "Rng-LC",
            "Ref-Freq (Hz)",
            "XmtrDly (nsec)",
            "RcvrDly (nsec)",
            "ScDly (nsec)",
        }

    def test_observable_column_is_renamed_to_a_unit_agnostic_name(self, msr_df):
        """Doppler is archived in Hz and range in RU, but both land in "Observed"."""
        assert "Observed" in msr_df.columns
        assert not {"Observed (Hz)", "Observed (RU)"} & set(msr_df.columns)

    def test_range_records_appended_below_the_doppler_table_are_read(self, msr_df):
        """``atdf2ascii`` restarts the header mid-file when it appends the range table.

        That second header line must be skipped rather than parsed as data, and
        the range rows below it must still be picked up.
        """
        assert len(msr_df) == MSR_ROW_COUNT
        assert msr_df["Data Type"].value_counts().to_dict() == MSR_ROWS_PER_DATA_TYPE

    def test_msr_time_tags_are_parsed_into_datetimes(self, msr_df):
        assert msr_df["time_tag (DT)"].iloc[0] == pd.Timestamp("1999-03-07 19:27:35")
        assert msr_df["time_tag (DT)"].iloc[-1] == pd.Timestamp("1999-03-11 22:39:53")

    def test_ramp_table_is_parsed_with_both_interval_bounds(self, ramp_df):
        assert len(ramp_df) == RAMP_ROW_COUNT
        assert ramp_df["Start-Time (DT)"].iloc[0] == pd.Timestamp("1999-03-07 11:46:54")
        assert ramp_df["End-Time (DT)"].iloc[-1] == pd.Timestamp("1999-03-12 11:45:00")
        assert (ramp_df["End-Time (DT)"] >= ramp_df["Start-Time (DT)"]).all()


# -----------------------------------------------------------------------------
# N-way range converter
# -----------------------------------------------------------------------------
class TestNwayRangeConverter:
    @pytest.fixture
    def converter(self) -> AtdfNwayRangeConverter:
        return AtdfNwayRangeConverter()

    def test_extracted_frame_carries_the_fields_needed_for_grouping(self, converter):
        extracted = converter.extract(msr_frame({"Data Type": "2-Way-Range"}), SPACECRAFT)
        assert set(extracted.columns) == {
            "epoch",
            "link_ends",
            "band",
            "link_end_delays",
            "lowest_ranging_component",
            "obs",
        }

    def test_each_observation_becomes_its_own_tracking_data_object(self, converter):
        """Range observations are not grouped, so that the range conversion factor can be stored per observation."""
        df = msr_frame(
            {"Data Type": "2-Way-Range", "Observed": 1.0},
            {"Data Type": "2-Way-Range", "Observed": 2.0},
        )
        tracking_data = converter.process(converter.extract(df, SPACECRAFT), EARTH)

        assert len(tracking_data) == 2
        assert [np.asarray(td.observations).item() for td in tracking_data] == [1.0, 2.0]

    def test_observation_carries_epoch_reference_link_end_and_time_scale(self, converter):
        df = msr_frame({"Data Type": "2-Way-Range"})
        (tracking_data,) = converter.process(converter.extract(df, SPACECRAFT), EARTH)

        assert tracking_data.observable_type == "DsnNWayRange"
        assert tracking_data.reference_link_end == "receiver"
        assert tracking_data.time_scale == "UTC"
        assert tracking_data.epochs[0].to_float() == pytest.approx(
            expected_epoch(pd.Timestamp("1999-03-07 19:27:35")), rel=0, abs=1e-6
        )

    def test_ancillary_settings_describe_bands_delays_and_ranging_component(self, converter):
        df = msr_frame(
            {
                "Data Type": "2-Way-Range",
                "UL": "X",
                "DL": "S",
                "Rng-LC": 21.0,
                "XmtrDly (nsec)": 100.0,
                "ScDly (nsec)": 200.0,
                "RcvrDly (nsec)": 50.0,
            }
        )
        (tracking_data,) = converter.process(converter.extract(df, SPACECRAFT), EARTH)

        assert tracking_data.get_ancillary_settings_string_vector()["frequency bands"] == [
            "X-band",
            "S-band",
        ]
        assert tracking_data.get_ancillary_settings_double()[
            "DSN sequential range lowest ranging component"
        ] == pytest.approx(21.0)
        assert tracking_data.get_ancillary_settings_double_vector()[
            "link ends time delays"
        ] == pytest.approx([1e-7, 2e-7, 5e-8])

    def test_fixture_yields_one_tracking_data_object_per_two_way_range_record(
        self, converter, msr_df
    ):
        tracking_data = converter.process(converter.extract(msr_df, SPACECRAFT), EARTH)
        assert len(tracking_data) == NWAY_RANGE_ROW_COUNT
        assert {td.observable_type for td in tracking_data} == {"DsnNWayRange"}

    def test_first_fixture_observation_matches_the_table(
        self, converter: AtdfNwayRangeConverter, msr_df
    ):
        """Test the observable value and epoch of the first record."""
        tracking_data = converter.process(converter.extract(msr_df, SPACECRAFT), EARTH)

        assert np.asarray(tracking_data[0].observations).item() == pytest.approx(
            FIRST_RANGE_RECORD["observable"]
        )
        assert tracking_data[0].epochs[0].to_float() == pytest.approx(
            expected_epoch(FIRST_RANGE_RECORD["epoch"]), rel=0, abs=1e-6
        )

    def test_fixture_shares_one_lowest_ranging_component(
        self, converter: AtdfNwayRangeConverter, msr_df
    ):
        tracking_data = converter.process(converter.extract(msr_df, SPACECRAFT), EARTH)
        components = {
            td.get_ancillary_settings_double()["DSN sequential range lowest ranging component"]
            for td in tracking_data
        }
        assert components == {float(FIXTURE_LOWEST_RANGING_COMPONENT)}


# -----------------------------------------------------------------------------
# N-way Doppler converter
# -----------------------------------------------------------------------------
class TestNwayDopplerConverter:
    @pytest.fixture
    def converter(self) -> AtdfNwayDopplerConverter:
        return AtdfNwayDopplerConverter()

    def test_two_and_three_way_doppler_records_are_extracted_together(self, converter, msr_df):
        assert len(converter.extract(msr_df, SPACECRAFT)) == NWAY_DOPPLER_ROW_COUNT

    def test_range_records_are_not_extracted(self, converter, msr_df):
        extracted = converter.extract(msr_df, SPACECRAFT)
        doppler_rows = msr_df[msr_df["Data Type"].str.contains("Doppler")]

        assert extracted["epoch"].to_list() == doppler_rows["time_tag (DT)"].to_list()

    def test_extracted_frame_carries_the_fields_needed_for_tracking_data(self, converter):
        extracted = converter.extract(msr_frame(), SPACECRAFT)
        assert set(extracted.columns) == {
            "epoch",
            "link_ends",
            "band",
            "link_end_delays",
            "count_time",
            "obs",
            "ref_freq",
            "ex",
        }

    def test_compatible_observations_are_merged_into_one_tracking_data_object(self, converter):
        df = msr_frame(
            {"time_tag (DT)": pd.Timestamp("1999-03-07 19:27:35"), "Observed": -1.0},
            {"time_tag (DT)": pd.Timestamp("1999-03-07 19:28:35"), "Observed": -2.0},
        )
        (tracking_data,) = converter.process(converter.extract(df, SPACECRAFT), EARTH)

        assert tracking_data.observable_type == "DsnNWayAveragedDoppler"
        assert np.asarray(tracking_data.observations).flatten() == pytest.approx([-1.0, -2.0])
        assert [epoch.to_float() for epoch in tracking_data.epochs] == pytest.approx(
            [
                expected_epoch(pd.Timestamp("1999-03-07 19:27:35")),
                expected_epoch(pd.Timestamp("1999-03-07 19:28:35")),
            ],
            rel=0,
            abs=1e-6,
        )

    def test_ancillary_settings_describe_bands_exciter_frequency_and_integration_time(
        self, converter
    ):
        df = msr_frame({"UL": "X", "DL": "X", "Ex": "S", "CT (sec)": 60.0})
        (tracking_data,) = converter.process(converter.extract(df, SPACECRAFT), EARTH)
        double_settings = tracking_data.get_ancillary_settings_double()

        assert tracking_data.get_ancillary_settings_string_vector()["frequency bands"] == [
            "X-band",
            "X-band",
        ]
        assert double_settings["DSN Doppler reference frequency"] == pytest.approx(2114118912.0)
        assert double_settings["Doppler observable integration time"] == pytest.approx(60.0)

    def test_three_way_observation_keeps_transmitting_and_receiving_station_apart(
        self, converter, msr_df
    ):
        """The fixture contains a 3-way pass, transmitted by DSS 45 and received by DSS 54."""
        tracking_data = converter.process(converter.extract(msr_df, SPACECRAFT), EARTH)
        three_way = [td for td in tracking_data if td.link_ends[0][0][1] != td.link_ends[2][0][1]]

        assert len(three_way) == 1
        assert three_way[0].link_ends == [
            ((EARTH, "DSS-45"), "transmitter"),
            ((SPACECRAFT, ""), "reflector_1"),
            ((EARTH, "DSS-54"), "receiver"),
        ]
        assert np.asarray(three_way[0].observations).size == MSR_ROWS_PER_DATA_TYPE["3-Way-Doppler"]

    def test_grouping_preserves_every_fixture_observation(self, converter, msr_df):
        tracking_data = converter.process(converter.extract(msr_df, SPACECRAFT), EARTH)
        total = sum(np.asarray(td.observations).size for td in tracking_data)

        assert total == NWAY_DOPPLER_ROW_COUNT

    def test_first_fixture_observation_matches_the_table(self, converter, msr_df):
        tracking_data = converter.process(converter.extract(msr_df, SPACECRAFT), EARTH)

        assert np.asarray(tracking_data[0].observations).flatten()[0] == pytest.approx(
            FIRST_DOPPLER_RECORD["observable"]
        )
        assert tracking_data[0].epochs[0].to_float() == pytest.approx(
            expected_epoch(FIRST_DOPPLER_RECORD["epoch"]), rel=0, abs=1e-6
        )
        assert tracking_data[0].get_ancillary_settings_double()[
            "DSN Doppler reference frequency"
        ] == pytest.approx(FIRST_DOPPLER_RECORD["reference_frequency"])


# -----------------------------------------------------------------------------
# Ramp converter
# -----------------------------------------------------------------------------
class TestRampConverter:
    @pytest.fixture
    def converter(self) -> AtdfRampConverter:
        return AtdfRampConverter()

    def test_extract_normalises_station_names(self, converter):
        extracted = converter.extract(ramp_frame({"Station": " DSS 34"}))
        assert extracted["station"].to_list() == ["DSS-34"]

    def test_extract_keeps_interval_bounds_frequency_and_rate(self, converter):
        extracted = converter.extract(
            ramp_frame({"Frequency (Hz)": 7.16e9, "Rate (Hz/sec)": -57.565952})
        )
        assert extracted["freq"].iloc[0] == pytest.approx(7.16e9)
        assert extracted["rate"].iloc[0] == pytest.approx(-57.565952)
        assert extracted["start_time"].iloc[0] == pd.Timestamp("1999-03-07 11:46:54")
        assert extracted["end_time"].iloc[0] == pd.Timestamp("1999-03-07 11:47:00")

    def test_one_supplementary_data_object_is_created_per_station(self, converter):
        df = ramp_frame({"Station": " DSS 34"}, {"Station": " DSS 45"}, {"Station": " DSS 34"})
        supplementary_data = converter.process(converter.extract(df), EARTH)

        assert {sd.reference_point_name for sd in supplementary_data} == {"DSS-34", "DSS-45"}
        assert {sd.body_name for sd in supplementary_data} == {EARTH}

    def test_ramps_are_stored_with_epoch_seconds_frequency_and_rate(self, converter):
        df = ramp_frame(
            {
                "Start-Time (DT)": pd.Timestamp("1999-03-07 11:46:54"),
                "End-Time (DT)": pd.Timestamp("1999-03-07 11:47:00"),
                "Frequency (Hz)": 7164234321.75112,
                "Rate (Hz/sec)": -57.565952,
            }
        )
        (supplementary_data,) = converter.process(converter.extract(df), EARTH)
        (ramp,) = supplementary_data.frequency_supplementary_data[0].frequency_ramps

        assert ramp.start_time == pytest.approx(
            expected_epoch(pd.Timestamp("1999-03-07 11:46:54")), rel=0, abs=1e-6
        )
        assert ramp.end_time - ramp.start_time == pytest.approx(6.0)
        assert ramp.start_frequency == pytest.approx(7164234321.75112, rel=0, abs=1e-6)
        assert ramp.frequency_rate == pytest.approx(-57.565952, rel=0, abs=1e-6)

    def test_every_fixture_ramp_is_assigned_to_its_station(self, converter, ramp_df):
        supplementary_data = converter.process(converter.extract(ramp_df), EARTH)
        ramps_per_station = {
            sd.reference_point_name: len(sd.frequency_supplementary_data[0].frequency_ramps)
            for sd in supplementary_data
        }
        assert ramps_per_station == RAMP_ROWS_PER_STATION


# -----------------------------------------------------------------------------
# Public entry point
# -----------------------------------------------------------------------------
class TestReadAtdfData:
    def test_decoding_and_conversion_are_chained(self, stub_atdf_to_ascii, tmp_path):
        tracking_data, supplementary_data = read_atdf_data(
            [Path(f"{FIXTURE_STEM}.tdf")], SPACECRAFT, output_dir=tmp_path
        )

        assert len(stub_atdf_to_ascii) == 1
        assert [td.observable_type for td in tracking_data].count(
            "DsnNWayRange"
        ) == NWAY_RANGE_ROW_COUNT
        assert (
            sum(
                np.asarray(td.observations).size
                for td in tracking_data
                if td.observable_type == "DsnNWayAveragedDoppler"
            )
            == NWAY_DOPPLER_ROW_COUNT
        )
        assert {sd.reference_point_name for sd in supplementary_data} == set(RAMP_ROWS_PER_STATION)

    def test_spacecraft_name_is_used_in_the_link_definitions(self, stub_atdf_to_ascii, tmp_path):
        tracking_data, _ = read_atdf_data(
            [Path(f"{FIXTURE_STEM}.tdf")], SPACECRAFT, output_dir=tmp_path
        )
        assert all(td.link_ends[1] == ((SPACECRAFT, ""), "reflector_1") for td in tracking_data)

    def test_reader_options_reach_atdf2ascii(self, stub_atdf_to_ascii, tmp_path):
        read_atdf_data(
            [Path(f"{FIXTURE_STEM}.tdf")],
            SPACECRAFT,
            output_dir=tmp_path,
            count_time=[60.0],
            proc_count=2,
            doppler_three_way=False,
        )

        (call,) = stub_atdf_to_ascii
        assert call["count_time"] == [60.0]
        assert call["proc_count"] == 2
        assert call["doppler_three_way"] is False


# -----------------------------------------------------------------------------
# Integration with the estimation interface
# -----------------------------------------------------------------------------
@pytest.fixture(scope="module")
def observation_collection(dsn_bodies):
    """The full fixture converted and handed to tudat's estimation interface."""
    tracking_data, _ = AtdfTrackingDataProcessor([Path(FIXTURE_STEM)], SPACECRAFT).process(
        FIXTURE_DIR
    )
    return create_observation_collection_from_tracking_data(tracking_data, dsn_bodies)


@pytest.fixture(scope="module")
def observation_sets(observation_collection):
    return observation_collection.get_single_observation_sets()


class TestObservationCollectionIntegration:
    """Conversion from ATDF TrackingData sets to ObservationCollection"""

    def test_observation_collection_observation_count(self, observation_collection):
        assert (
            observation_collection.concatenated_observations.size
            == NWAY_DOPPLER_ROW_COUNT + NWAY_RANGE_ROW_COUNT
        )

    def test_both_dsn_observable_types_are_present(self, observation_sets):
        types = {obs_set.observable_type for obs_set in observation_sets}
        assert types == {
            ObservableType.dsn_n_way_range_type,
            ObservableType.dsn_n_way_averaged_doppler_type,
        }

    def test_range_ancillary_settings_conversion(self, observation_sets):
        variables = ancillary_settings.ObservationAncillarySimulationVariable
        range_sets = [
            set
            for set in observation_sets
            if set.observable_type == ObservableType.dsn_n_way_range_type
        ]
        assert len(range_sets) == NWAY_RANGE_ROW_COUNT

        # Every range record in the fixture shares one band pair, delay triplet and
        # lowest ranging component
        for obs_set in range_sets:
            settings = obs_set.ancillary_settings
            assert settings.get_float_settings(
                variables.sequential_range_lowest_ranging_component
            ) == float(FIXTURE_LOWEST_RANGING_COMPONENT)
            assert settings.get_float_list_settings(variables.frequency_bands) == [
                ancillary_settings.FrequencyBands.x_band,
                ancillary_settings.FrequencyBands.x_band,
            ]
            assert settings.get_float_list_settings(variables.link_ends_delays) == [0.0, 0.0, 0.0]

    def test_doppler_ancillary_settings_conversion(self, observation_sets):
        variables = ancillary_settings.ObservationAncillarySimulationVariable
        doppler_sets = [
            set
            for set in observation_sets
            if set.observable_type == ObservableType.dsn_n_way_averaged_doppler_type
        ]
        assert len(doppler_sets) == len(FIXTURE_DOPPLER_GROUP_SIZES)

        for obs_set in doppler_sets:
            settings = obs_set.ancillary_settings
            assert (
                settings.get_float_settings(variables.doppler_integration_time)
                == FIXTURE_COUNT_TIME
            )
            # The fixture's exciter band is S, while the up- and downlink are X band.
            assert (
                settings.get_float_settings(variables.reception_reference_frequency_band)
                == ancillary_settings.FrequencyBands.s_band
            )
            assert settings.get_float_list_settings(variables.frequency_bands) == [
                ancillary_settings.FrequencyBands.x_band,
                ancillary_settings.FrequencyBands.x_band,
            ]
            assert settings.get_float_list_settings(variables.link_ends_delays) == [0.0, 0.0, 0.0]

        # Reference frequency differs between passes
        reference_frequencies = {
            s.ancillary_settings.get_float_settings(variables.doppler_reference_frequency)
            for s in doppler_sets
        }
        assert len(reference_frequencies) == len(doppler_sets)
        assert FIRST_DOPPLER_RECORD["reference_frequency"] in reference_frequencies

    def test_link_definitions_resolve_to_dsn_ground_stations(self, observation_sets):
        three_way = [
            obs_set
            for obs_set in observation_sets
            if obs_set.link_definition.link_end_id(links.LinkEndType.transmitter).reference_point
            != obs_set.link_definition.link_end_id(links.LinkEndType.receiver).reference_point
        ]
        assert len(three_way) == 1

        link_definition = three_way[0].link_definition
        assert link_definition.link_end_id(links.LinkEndType.transmitter).body_name == EARTH
        assert (
            link_definition.link_end_id(links.LinkEndType.transmitter).reference_point == "DSS-45"
        )
        assert link_definition.link_end_id(links.LinkEndType.reflector1).body_name == SPACECRAFT
        assert link_definition.link_end_id(links.LinkEndType.receiver).reference_point == "DSS-54"


if __name__ == "__main__":
    pytest.main([__file__])
