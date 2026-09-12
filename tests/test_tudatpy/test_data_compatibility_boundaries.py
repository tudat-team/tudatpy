"""Check the Python compatibility boundary without loading the Tudat extension."""

import ast
import importlib
import inspect
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
import warnings

import pytest

SOURCE = Path(__file__).resolve().parents[2] / "src" / "tudatpy"


@pytest.fixture
def compatibility(monkeypatch):
    """Load source compatibility modules independently of the installed extension."""
    # Use the actual source modules, even when the installed extension predates
    # this refactor. Restore the original module cache after each test.
    for name in list(sys.modules):
        if name == "tudatpy" or name.startswith("tudatpy."):
            monkeypatch.delitem(sys.modules, name)
    package = ModuleType("tudatpy")
    package.__path__ = [str(SOURCE)]
    monkeypatch.setitem(sys.modules, "tudatpy", package)
    try:
        yield importlib.import_module("tudatpy.data._compat")
    finally:
        for name in list(sys.modules):
            if name.startswith("tudatpy."):
                del sys.modules[name]


@pytest.mark.parametrize(
    "module_name,symbol,replacement",
    [
        (
            "tudatpy.data.discos.discos",
            "DiscosQuery",
            "tudatpy.data_input.environment_data.discos.DiscosQuery",
        ),
        (
            "tudatpy.data.sbdb.sbdb",
            "SBDBquery",
            "tudatpy.data_input.environment_data.sbdb.SBDBquery",
        ),
        (
            "tudatpy.data.ancillary.ancillary_downloader",
            "download_ionex",
            "tudatpy.data_input.data_retrieval.media_corrections.download_ionex",
        ),
        (
            "tudatpy.data.mission_data_downloader.mission_data_downloader",
            "LoadPDS",
            "tudatpy.data_input.data_retrieval.missions.LoadPDS",
        ),
        (
            "tudatpy.data",
            "get_resource_path",
            "tudatpy.data_input.resource_paths.get_resource_path",
        ),
        ("tudatpy.data.mpc", "BatchMPC", "tudatpy.data_input.tracking_data.mpc.BatchMPC"),
        ("tudatpy.data.mpc.mpc", "BatchMPC", "tudatpy.data_input.tracking_data.mpc.BatchMPC"),
        (
            "tudatpy.data.mpc.mpc",
            "get_biases_EFCC18",
            "tudatpy.data_input.tracking_data.optical_utilities.get_biases_EFCC18",
        ),
        (
            "tudatpy.data.horizons.horizons",
            "HorizonsQuery",
            "tudatpy.data_input.environment_data.horizons.HorizonsQuery",
        ),
        (
            "tudatpy.data.spacetrack.spacetrack",
            "OMMUtils",
            "tudatpy.data_input.environment_data.spacetrack.OMMUtils",
        ),
        (
            "tudatpy.data.processTrk234",
            "Trk234Processor",
            "tudatpy.data_input.tracking_data.tnf.read_tnf_data",
        ),
        (
            "tudatpy.data.processTrk234.processor",
            "Trk234Processor",
            "tudatpy.data_input.tracking_data.tnf.read_tnf_data",
        ),
        (
            "tudatpy.data.processTrk234.converters.derivedDoppler",
            "DerivedDopplerConverter",
            "tudatpy.data_input.tracking_data.tnf.read_tnf_data",
        ),
        (
            "tudatpy.data.processTrk234.converters.derivedSraRange",
            "DerivedSraRangeConverter",
            "tudatpy.data_input.tracking_data.tnf.read_tnf_data",
        ),
    ],
)
@pytest.mark.parametrize("access", ["attribute", "from_import"])
def test_legacy_access_warns_at_user_line_with_public_replacement(
    compatibility, monkeypatch, module_name, symbol, replacement, access
):
    """Warn at the real user access/import line and name a public migration target."""
    module = importlib.import_module(module_name)
    sentinel = object()
    target_name = module._ALIASES[symbol].rsplit(".", 1)[1]
    # Only the destination needs a stub; alias lookup and Python's actual
    # from-import machinery both run unchanged.
    monkeypatch.setattr(
        compatibility,
        "importlib",
        SimpleNamespace(import_module=lambda name: SimpleNamespace(**{target_name: sentinel})),
    )
    statement = (
        f"value = module.{symbol}"
        if access == "attribute"
        else f"from {module_name} import {symbol} as value"
    )
    namespace = {"module": module, "__name__": "user_script"}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", DeprecationWarning)
        exec(compile("\n" * 6 + statement, "user_script.py", "exec"), namespace)
    # Import machinery may look up an alias twice; every warning must name the user line.
    assert namespace["value"] is sentinel
    assert caught
    for warning in caught:
        assert warning.category is DeprecationWarning
        assert (
            str(warning.message)
            == f"{module_name}.{symbol} is deprecated. Use {replacement} instead."
        )
        assert warning.filename == "user_script.py"
        assert warning.lineno == 7


def test_missing_names_and_introspection_do_not_warn(compatibility):
    """Keep discovery and unknown-attribute checks silent while exposing legacy names."""
    module = importlib.import_module("tudatpy.data.mpc.mpc")
    # Listing exports and probing absent names must not look up deprecated targets.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        assert set(module.__all__) <= set(dir(module))
        with pytest.raises(AttributeError, match="no attribute 'missing'"):
            getattr(module, "missing")
    assert not caught


def test_nested_legacy_warning_points_to_user_call(compatibility):
    """Attribute warnings from nested legacy helpers to the external call site."""
    legacy = ModuleType("tudatpy.data.mpc._legacy")
    legacy.warn = compatibility.warn_custom_deprecation
    exec(
        "def to_tudat():\n"
        "    warn('tudatpy.data.mpc', 'get_weights_VFCC17', 'Request add_weights=True.')\n",
        legacy.__dict__,
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", DeprecationWarning)
        expected_line = inspect.currentframe().f_lineno + 1
        legacy.to_tudat()
    # Skip the nested compatibility frame while retaining the warning and migration hint.
    assert len(caught) == 1
    assert caught[0].filename == __file__
    assert caught[0].lineno == expected_line
    assert "get_weights_VFCC17 is deprecated" in str(caught[0].message)


def test_non_one_to_one_reader_warns_at_call_and_preserves_arguments(compatibility, monkeypatch):
    """Warn on legacy reader calls without changing their arguments or result."""
    data = importlib.import_module("tudatpy.data")
    kernel = ModuleType("tudatpy.kernel")
    calls = []
    kernel.data = SimpleNamespace(read_ifms_file=lambda *args: calls.append(args) or "parsed")
    monkeypatch.setitem(sys.modules, "tudatpy.kernel", kernel)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", DeprecationWarning)
        expected_line = inspect.currentframe().f_lineno + 1
        result = data.read_ifms_file("input.tab", False, True)
    # Forward the old signature unchanged and attribute the warning to this call.
    assert result == "parsed"
    assert calls == [("input.tab", False, True)]
    assert len(caught) == 1
    assert caught[0].filename == __file__
    assert caught[0].lineno == expected_line
    assert "read_ifms_data" in str(caught[0].message)


@pytest.mark.parametrize("action,count", [("default", 1), ("ignore", 0), ("always", 2)])
def test_warning_filters_and_caller_registry_are_respected(compatibility, action, count):
    """Respect user warning filters and deduplicate repeated calls at one source line."""
    namespace = {"warn": compatibility.warn_custom_deprecation, "__name__": "user_script"}
    code = compile(
        "warn('tudatpy.data', 'old_reader', 'Use the current reader.')", "user.py", "exec"
    )
    # Execute the same call twice without resetting the caller's warning registry.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter(action, DeprecationWarning)
        exec(code, namespace)
        exec(code, namespace)
    assert len(caught) == count


def test_warning_filter_can_promote_deprecation_to_error(compatibility):
    """Allow users to reject deprecated calls through the normal warnings-as-errors filter."""
    # The compatibility helper must not suppress or replace the requested exception.
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        with pytest.raises(DeprecationWarning, match="old_reader is deprecated"):
            compatibility.warn_custom_deprecation(
                "tudatpy.data", "old_reader", "Use its replacement."
            )


@pytest.mark.filterwarnings("ignore:invalid escape sequence:DeprecationWarning")
@pytest.mark.filterwarnings("ignore:invalid escape sequence:SyntaxWarning")
def test_current_python_api_does_not_import_removable_data_package():
    """Keep current Python imports independent of the deprecated package and kernel bridge."""
    dependencies = []
    # Inspect source imports, resolving relative paths without importing extension modules.
    for path in SOURCE.rglob("*.py"):
        if path.is_relative_to(SOURCE / "data"):
            continue
        relative = path.relative_to(SOURCE).with_suffix("")
        parts = list(relative.parts)
        if parts[-1] == "__init__":
            parts.pop()
        else:
            parts = parts[:-1]
        package = ".".join(["tudatpy", *parts])
        for node in ast.walk(ast.parse(path.read_text(), filename=str(path))):
            names = []
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                module = node.module or ""
                if node.level:
                    module = importlib.util.resolve_name("." * node.level + module, package)
                names = [module, *(f"{module}.{alias.name}" for alias in node.names)]
            for name in names:
                if (
                    name == "tudatpy.data"
                    or name.startswith("tudatpy.data.")
                    or name == "tudatpy.kernel.data"
                ):
                    dependencies.append(f"{path.relative_to(SOURCE)}:{node.lineno}: {name}")
    # New backwards dependencies would prevent deleting the compatibility package as a unit.
    assert (
        not dependencies
    ), "Current APIs must remain independent of deprecated data adapters:\n" + "\n".join(
        dependencies
    )
