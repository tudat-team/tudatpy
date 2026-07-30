from pathlib import Path

import tudatpy


def test_tudatpy_version_matches_version_file():
    repository_root = Path(__file__).resolve().parents[2]
    expected_version = (
        (repository_root / "version").read_text(encoding="utf-8").splitlines()[0].strip()
    )

    assert tudatpy.__version__ == expected_version
