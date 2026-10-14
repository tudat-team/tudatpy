"""Pinned provider excerpts stay unchanged; provenance lists exact source ranges."""

import hashlib
import json
from pathlib import Path


def test_ocean_pole_gravity_provider_extracts():
    fixtures = Path(__file__).parents[1] / "test_tudat/reference_files/ocean_pole_gravity"
    manifest = json.loads((fixtures / "manifest.json").read_text())
    for name, record in manifest.items():
        assert (
            hashlib.sha256((fixtures / name).read_bytes()).hexdigest() == record["extract_sha256"]
        )
        assert len(record["source_sha256"]) == 64
        assert record["source_lines"]
        assert record["edits"] == []
    lines = (fixtures / "desaiscopole2004.txt").read_text().splitlines()
    assert sum(len(line.split()) == 9 and line.split()[0].isdigit() for line in lines) == 366
