"""Optional-feature selection must follow the loaded binary, not missing APIs."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "serialization_policy_conftest", Path(__file__).with_name("conftest.py")
)
conftest = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(conftest)


@pytest.mark.parametrize("enabled", [False, True])
def test_explicit_boolean_build_flag(enabled):
    kernel = SimpleNamespace(_build_with_serialization=enabled)
    assert conftest._serialization_enabled(kernel) is enabled


@pytest.mark.parametrize("value", [None, 0, 1, "OFF", "ON"])
def test_unknown_or_non_boolean_build_flag_fails_loudly(value):
    with pytest.raises(pytest.UsageError, match="Rebuild the kernel"):
        conftest._serialization_enabled(SimpleNamespace(_build_with_serialization=value))


def test_missing_build_flag_is_not_treated_as_disabled():
    with pytest.raises(pytest.UsageError, match="import path"):
        conftest._serialization_enabled(SimpleNamespace())


class _Item:
    def __init__(self, needs_serialization):
        self.needs_serialization = needs_serialization
        self.added = []

    def get_closest_marker(self, name):
        if name == "requires_serialization" and self.needs_serialization:
            return object()
        return None

    def add_marker(self, marker):
        self.added.append(marker)


@pytest.mark.parametrize("enabled", [False, True])
def test_collection_only_skips_optional_tests_when_build_option_is_off(monkeypatch, enabled):
    from tudatpy import kernel

    monkeypatch.setattr(kernel, "_build_with_serialization", enabled)
    marked, ordinary = _Item(True), _Item(False)
    config = SimpleNamespace(
        getoption=lambda name: False,
        pluginmanager=SimpleNamespace(get_plugin=lambda name: None),
    )
    conftest.pytest_collection_modifyitems(config, [marked, ordinary])
    assert not ordinary.added
    if enabled:
        # ON must run the tests even if a serialization method is broken/missing.
        assert not marked.added
    else:
        assert len(marked.added) == 1
        assert marked.added[0].name == "skip"
        assert "TUDAT_BUILD_WITH_SERIALIZATION=OFF" in marked.added[0].kwargs["reason"]


def test_loaded_build_flag_and_basic_binding_agree():
    from tudatpy import kernel
    from tudatpy.math import root_finders

    enabled = conftest._serialization_enabled(kernel)
    # This assertion fails for an ON build with a broken registration; it never
    # changes the feature flag or converts a broken binding into a skip.
    assert hasattr(root_finders.RootFinderSettings, "save_to_json") == enabled
    assert hasattr(root_finders.RootFinderSettings, "load_from_json") == enabled
