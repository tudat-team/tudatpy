import importlib.util
import json
from pathlib import Path

import pytest

MODULE_PATH = Path(__file__).parents[2] / "src" / "tudatpy" / "resource_catalog.py"
SPEC = importlib.util.spec_from_file_location("resource_catalog_under_test", MODULE_PATH)
resource_catalog = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(resource_catalog)


class Response:
    def __init__(self, text):
        self.text = text

    def raise_for_status(self):
        return None


def manifest(record_id, archive, path):
    return (
        f"10.5281/zenodo.{record_id}\n"
        f"{archive},,\n"
        f"{path},2026-08-03T00:00:00+00:00,https://example.test/resource/{path}\n"
    )


def test_missing_default_catalog_explains_how_to_create_it(monkeypatch, tmp_path):
    """Direct users to manifest mode when no local catalog exists."""
    missing_catalog = tmp_path / "resource_catalog.csv"
    monkeypatch.setattr(resource_catalog, "USER_RESOURCE_CATALOG", missing_catalog)

    with pytest.raises(FileNotFoundError, match="--mode manifest") as error:
        resource_catalog.load_resource_catalog()

    assert str(missing_catalog) in str(error.value)


def test_missing_explicit_catalog_includes_catalog_command(tmp_path):
    """Include the selected path in instructions for creating a custom catalog."""
    missing_catalog = tmp_path / "custom.csv"

    with pytest.raises(FileNotFoundError, match="--catalog") as error:
        resource_catalog.load_resource_catalog(missing_catalog)

    assert str(missing_catalog) in str(error.value)


def test_update_resource_catalog_writes_catalog_and_hashes(monkeypatch, tmp_path):
    responses = iter(
        [
            Response(manifest("21261530", "rolling.tar.gz", "rolling.dat")),
            Response("a" * 64 + "  rolling.tar.gz\n"),
            Response(manifest("21277280", "static.tar.gz", "static.dat")),
            Response("b" * 64 + " *static.tar.gz\n"),
        ]
    )
    monkeypatch.setattr(resource_catalog.requests, "get", lambda *args, **kwargs: next(responses))
    catalog_path = tmp_path / "resource_catalog.csv"
    hash_path = tmp_path / "resource_hashes.json"

    count = resource_catalog.update_resource_catalog(catalog_path, hash_path)

    assert count == 2
    catalog = resource_catalog.load_resource_catalog(catalog_path)
    assert set(catalog) == {"rolling.dat", "static.dat"}
    assert json.loads(hash_path.read_text()) == {
        "rolling.tar.gz": "a" * 64,
        "static.tar.gz": "b" * 64,
    }


def test_update_resource_catalog_preserves_files_on_invalid_hashes(monkeypatch, tmp_path):
    catalog_path = tmp_path / "resource_catalog.csv"
    hash_path = tmp_path / "resource_hashes.json"
    catalog_path.write_text("old catalog\n")
    hash_path.write_text("old hashes\n")
    responses = iter(
        [
            Response(manifest("21261530", "rolling.tar.gz", "rolling.dat")),
            Response("not a SHA256SUM file\n"),
        ]
    )
    monkeypatch.setattr(resource_catalog.requests, "get", lambda *args, **kwargs: next(responses))

    with pytest.raises(ValueError, match="Invalid SHA256SUM"):
        resource_catalog.update_resource_catalog(catalog_path, hash_path)

    assert catalog_path.read_text() == "old catalog\n"
    assert hash_path.read_text() == "old hashes\n"
