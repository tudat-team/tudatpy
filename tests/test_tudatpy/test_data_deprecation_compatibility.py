import importlib

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


def _resolve_dotted_name(dotted_name):
    module_name, object_name = dotted_name.rsplit(".", 1)
    return getattr(importlib.import_module(module_name), object_name)


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


def test_deprecated_crd_single_file_reader_warns_and_delegates(monkeypatch):
    import tudatpy.data as data
    import tudatpy.data_input.tracking_data.slr as slr

    def read_slr_data(file_names):
        return ("parsed slr", file_names)

    monkeypatch.setattr(slr, "read_slr_data", read_slr_data)

    with pytest.warns(DeprecationWarning, match="read_slr_data"):
        parsed_data = data.read_crd_file("example.crd")

    assert parsed_data == ("parsed slr", ["example.crd"])


def test_deprecated_tracking_txt_single_file_reader_warns_and_delegates(monkeypatch):
    import tudatpy.data as data
    import tudatpy.data_input.tracking_data.generic_text_file as generic_text_file

    calls = []

    def read_generic_text_data(*args):
        calls.append(args)
        return ["parsed tracking text"]

    monkeypatch.setattr(generic_text_file, "read_generic_text_data", read_generic_text_data)

    with pytest.warns(DeprecationWarning, match="read_generic_text_data"):
        parsed_data = data.read_tracking_txt_file(
            "tracking.txt",
            ["year", "month"],
            data_filter_method="filter",
        )

    assert parsed_data == "parsed tracking text"
    assert calls == [
        (
            ["tracking.txt"],
            ["year", "month"],
            "#",
            ",:\t ",
            False,
            "filter",
        )
    ]


@pytest.mark.parametrize(
    ("old_name", "replacement"),
    (
        ("read_ifms_file", "read_ifms_data"),
        ("read_fdets_file", "read_fdets_data"),
        (
            "set_dsn_weather_data_in_ground_stations",
            "set_dsn_weather_data_in_ground_station_settings",
        ),
        (
            "set_estrack_weather_data_in_ground_stations",
            "set_estrack_weather_data_in_ground_station_settings",
        ),
    ),
)
def test_non_one_to_one_deprecated_data_functions_warn_and_stop(old_name, replacement):
    import tudatpy.data as data

    old_function = getattr(data, old_name)

    with pytest.warns(DeprecationWarning, match=replacement):
        with pytest.raises(NotImplementedError):
            old_function()


def test_data_mpc_weight_helper_is_not_reintroduced():
    with pytest.raises(AttributeError):
        getattr(importlib.import_module("tudatpy.data.mpc"), "get_weights_VFCC17")
