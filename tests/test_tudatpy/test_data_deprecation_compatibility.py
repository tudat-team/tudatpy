import importlib
import inspect
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

COMPATIBILITY_MODULES = (
    "tudatpy.data",
    "tudatpy.data.ancillary",
    "tudatpy.data.ancillary.ancillary_downloader",
    "tudatpy.data.coma_model",
    "tudatpy.data.discos",
    "tudatpy.data.discos.discos",
    "tudatpy.data.horizons",
    "tudatpy.data.horizons.horizons",
    "tudatpy.data.mission_data_downloader",
    "tudatpy.data.mission_data_downloader.mission_data_downloader",
    "tudatpy.data.mpc",
    "tudatpy.data.mpc.mpc",
    "tudatpy.data.mpc.parser_80col",
    "tudatpy.data.mpc.parser_80col.parsers",
    "tudatpy.data.mpc.parser_80col.unpackers",
    "tudatpy.data.processTrk234",
    "tudatpy.data.processTrk234.processor",
    "tudatpy.data.processTrk234.converters",
    "tudatpy.data.processTrk234.converters.converter",
    "tudatpy.data.processTrk234.converters.derivedDoppler",
    "tudatpy.data.processTrk234.converters.derivedSraRange",
    "tudatpy.data.processTrk234.converters.radioBase",
    "tudatpy.data.processTrk234.converters.ramp",
    "tudatpy.data.sbdb",
    "tudatpy.data.sbdb.sbdb",
    "tudatpy.data.spacetrack",
    "tudatpy.data.spacetrack.spacetrack",
)

DEVELOP_DATA_TOP_LEVEL_NAMES = (
    "CrdFullRateRecord",
    "CrdMeteoRecord",
    "CrdNormalPointRecord",
    "CrdPass",
    "CrdPassConfigurationData",
    "CrdPassData",
    "IlrsStationRegistryEntry",
    "OdfCommonDataBlock",
    "OdfDataBlock",
    "OdfDataSpecificBlock",
    "OdfDataType",
    "OdfDopplerDataBlock",
    "OdfRampBlock",
    "OdfRawFileContents",
    "SinexStationEccentricity",
    "SinexStationState",
    "SolarActivityContainer",
    "SolarActivityData",
    "TrackingDataType",
    "TrackingTxtFileContents",
    "TrackingTxtFileReadFilterType",
    "coma_model",
    "convert_crd_two_way_time_of_flight_to_slr_range",
    "convert_sinex_datetime_to_seconds_since_epoch",
    "extract_full_rate_measurements",
    "extract_full_rate_measurements_from_passes",
    "extract_normal_point_measurements",
    "extract_normal_point_measurements_from_passes",
    "get_atmosphere_tables_path",
    "get_earth_orientation_path",
    "get_ephemeris_path",
    "get_gravity_models_path",
    "get_quadrature_path",
    "get_resource_path",
    "get_space_weather_path",
    "get_spice_kernel_path",
    "get_station_wavelengths",
    "grail_antenna_file_reader",
    "grail_mass_level_0_file_reader",
    "grail_mass_level_1_file_reader",
    "group_crd_data_per_station",
    "group_crd_data_per_target",
    "read_crd_file",
    "read_crd_files",
    "read_domes_id_numbers",
    "read_fdets_file",
    "read_ground_station_names",
    "read_ifms_file",
    "read_ilrs_station_registry_from_sinex_site_id",
    "read_matrix_history_from_file",
    "read_monument_numbers",
    "read_odf_file",
    "read_sinex_station_data",
    "read_sinex_station_eccentricities",
    "read_solar_activity_data",
    "read_tracking_txt_file",
    "read_vector_history_from_file",
    "set_dsn_weather_data_in_ground_stations",
    "set_estrack_weather_data_in_ground_stations",
)


def _resolve_dotted_name(dotted_name):
    module_name, object_name = dotted_name.rsplit(".", 1)
    return getattr(importlib.import_module(module_name), object_name)


def _access_deprecated_data_alias_for_warning_location():
    import tudatpy.data as data

    expected_lineno = inspect.currentframe().f_lineno + 1
    deprecated_object = data.get_resource_path
    return expected_lineno, deprecated_object


def _call_deprecated_non_one_to_one_function_for_warning_location():
    import tudatpy.data as data
    import tudatpy.kernel.data as kernel_data

    with patch.object(kernel_data, "read_ifms_file", lambda *args, **kwargs: "parsed ifms"):
        expected_lineno = inspect.currentframe().f_lineno + 1
        data.read_ifms_file("ifms.tab")
    return expected_lineno


def _compatibility_cases():
    cases = []
    for module_name in COMPATIBILITY_MODULES:
        module = importlib.import_module(module_name)
        for old_name, new_name in module._ALIASES.items():
            cases.append((module_name, old_name, new_name))
    return cases


def test_kernel_data_submodule_is_exposed():
    import tudatpy.kernel as kernel

    assert hasattr(kernel, "data")


def test_develop_data_top_level_names_are_covered_by_compatibility_layer():
    import tudatpy.data as data

    covered_names = set(data._ALIASES) | data._CUSTOM_DEPRECATION_NAMES | {"coma_model"}

    assert set(DEVELOP_DATA_TOP_LEVEL_NAMES) <= covered_names


@pytest.mark.parametrize(
    ("module_name", "old_name", "new_name"),
    _compatibility_cases(),
)
def test_deprecated_data_alias_resolves_to_new_object(module_name, old_name, new_name):
    module = importlib.import_module(module_name)
    expected_object = _resolve_dotted_name(new_name)

    with pytest.warns(DeprecationWarning, match="is deprecated"):
        deprecated_object = getattr(module, old_name)

    assert deprecated_object is expected_object


@pytest.mark.parametrize(
    ("module_name", "old_name", "new_name"),
    (
        (
            "tudatpy.data",
            "get_resource_path",
            "tudatpy.data_input.resource_paths.get_resource_path",
        ),
        (
            "tudatpy.data.horizons",
            "HorizonsQuery",
            "tudatpy.data_input.environment_data.horizons.HorizonsQuery",
        ),
        ("tudatpy.data.mpc", "BatchMPC", "tudatpy.data_input.tracking_data.mpc.BatchMPC"),
        (
            "tudatpy.data.processTrk234",
            "Trk234Processor",
            "tudatpy.data_input.tracking_data.tnf.TnfTrackingDataProcessor",
        ),
    ),
)
def test_from_import_deprecated_data_alias_warns(module_name, old_name, new_name):
    module = importlib.import_module(module_name)
    expected_object = _resolve_dotted_name(new_name)

    with pytest.warns(DeprecationWarning, match=f"{module_name}.{old_name} is deprecated"):
        imported_object = getattr(module, old_name)

    assert imported_object is expected_object


def test_deprecated_data_alias_warning_points_to_user_code():
    with pytest.warns(DeprecationWarning, match="get_resource_path is deprecated") as record:
        expected_lineno, deprecated_object = _access_deprecated_data_alias_for_warning_location()

    assert deprecated_object is _resolve_dotted_name(
        "tudatpy.data_input.resource_paths.get_resource_path"
    )
    assert Path(record[0].filename).resolve() == Path(__file__).resolve()
    assert record[0].lineno == expected_lineno


def test_deprecated_non_one_to_one_warning_points_to_user_code():
    with pytest.warns(DeprecationWarning, match="read_ifms_data") as record:
        expected_lineno = _call_deprecated_non_one_to_one_function_for_warning_location()

    assert Path(record[0].filename).resolve() == Path(__file__).resolve()
    assert record[0].lineno == expected_lineno


def test_deprecated_observations_wrapper_alias_warns():
    module = importlib.import_module("tudatpy.estimation.observations_setup.observations_wrapper")
    observations = importlib.import_module("tudatpy.estimation.observations")

    with pytest.warns(
        DeprecationWarning,
        match=(
            "tudatpy.estimation.observations_setup.observations_wrapper."
            "create_observation_collection_from_tracking_data is deprecated"
        ),
    ):
        imported_object = getattr(module, "create_observation_collection_from_tracking_data")

    assert imported_object is observations.create_observation_collection_from_tracking_data


def test_deprecated_crd_single_file_reader_warns_and_delegates():
    import tudatpy.data as data
    import tudatpy.kernel.data as kernel_data

    def read_crd_file(file_name):
        return ("parsed slr", file_name)

    with patch.object(kernel_data, "read_crd_file", read_crd_file):
        with pytest.warns(DeprecationWarning, match="read_slr_data"):
            parsed_data = data.read_crd_file("example.crd")

    assert parsed_data == ("parsed slr", "example.crd")


def test_deprecated_tracking_txt_single_file_reader_warns_and_delegates():
    import tudatpy.data as data
    import tudatpy.kernel.data as kernel_data

    calls = []

    def read_tracking_txt_file(*args):
        calls.append(args)
        return "parsed tracking text"

    with patch.object(kernel_data, "read_tracking_txt_file", read_tracking_txt_file):
        with pytest.warns(DeprecationWarning, match="read_generic_text_data"):
            parsed_data = data.read_tracking_txt_file(
                "tracking.txt",
                ["year", "month"],
                data_filter_method="filter",
            )

    assert parsed_data == "parsed tracking text"
    assert calls == [
        (
            "tracking.txt",
            ["year", "month"],
            "#",
            ",:\t ",
            False,
            "filter",
        )
    ]


def test_deprecated_ifms_raw_reader_warns_and_delegates():
    import tudatpy.data as data
    import tudatpy.kernel.data as kernel_data

    calls = []

    def read_ifms_file(*args):
        calls.append(args)
        return "parsed ifms"

    with patch.object(kernel_data, "read_ifms_file", read_ifms_file):
        with pytest.warns(DeprecationWarning, match="read_ifms_data"):
            parsed_data = data.read_ifms_file("ifms.tab", False, True)

    assert parsed_data == "parsed ifms"
    assert calls == [("ifms.tab", False, True)]


def test_deprecated_fdets_raw_reader_warns_and_delegates():
    import tudatpy.data as data
    import tudatpy.kernel.data as kernel_data

    calls = []

    def read_fdets_file(*args):
        calls.append(args)
        return "parsed fdets"

    with patch.object(kernel_data, "read_fdets_file", read_fdets_file):
        with pytest.warns(DeprecationWarning, match="read_fdets_data"):
            parsed_data = data.read_fdets_file("fdets.txt", ["utc_datetime_string"])

    assert parsed_data == "parsed fdets"
    assert calls == [("fdets.txt", ["utc_datetime_string"])]


def test_deprecated_dsn_weather_setter_warns_and_delegates():
    import tudatpy.data as data
    import tudatpy.kernel.data as kernel_data

    calls = []

    def set_dsn_weather_data_in_ground_stations(*args, **kwargs):
        calls.append((args, kwargs))
        return "dsn weather set"

    with patch.object(
        kernel_data,
        "set_dsn_weather_data_in_ground_stations",
        set_dsn_weather_data_in_ground_stations,
    ):
        with pytest.warns(
            DeprecationWarning, match="set_dsn_weather_data_in_ground_station_settings"
        ):
            result = data.set_dsn_weather_data_in_ground_stations(
                "bodies", ["weather.tab"], body_with_ground_stations_name="Mars"
            )

    assert result == "dsn weather set"
    assert calls == [
        (
            ("bodies", ["weather.tab"]),
            {"body_with_ground_stations_name": "Mars"},
        )
    ]


def test_deprecated_estrack_weather_setter_warns_and_delegates():
    import tudatpy.data as data
    import tudatpy.kernel.data as kernel_data

    calls = []

    def set_estrack_weather_data_in_ground_stations(*args, **kwargs):
        calls.append((args, kwargs))
        return "estrack weather set"

    with patch.object(
        kernel_data,
        "set_estrack_weather_data_in_ground_stations",
        set_estrack_weather_data_in_ground_stations,
    ):
        with pytest.warns(
            DeprecationWarning, match="set_estrack_weather_data_in_ground_station_settings"
        ):
            result = data.set_estrack_weather_data_in_ground_stations(
                "bodies", ["weather.tab"], "NWNORCIA", body_with_ground_stations_name="Mars"
            )

    assert result == "estrack weather set"
    assert calls == [
        (
            ("bodies", ["weather.tab"], "NWNORCIA"),
            {"body_with_ground_stations_name": "Mars"},
        )
    ]


@pytest.mark.parametrize("module_name", ("tudatpy.data.mpc", "tudatpy.data.mpc.mpc"))
def test_deprecated_mpc_vfcc17_weight_helper_warns_and_works(module_name):
    module = importlib.import_module(module_name)

    class ObservatoryCodes:
        def to_pandas(self):
            return pd.DataFrame({"Code": ["704"], "Longitude": [0.0]})

    from astroquery.mpc import MPC

    with patch.object(MPC, "get_observatory_codes", lambda: ObservatoryCodes()):
        with pytest.warns(DeprecationWarning, match="no one-to-one equivalent"):
            weights = module.get_weights_VFCC17(
                MPC_codes=["433"],
                epoch=[2459000.0],
                observation_type=["C"],
                observatory=["704"],
                star_catalog=["U"],
            )

    expected_weight = 1.0 / np.square(np.deg2rad(1.0 / 3600.0))
    assert weights.shape == (1,)
    assert np.allclose(weights[0], expected_weight)
