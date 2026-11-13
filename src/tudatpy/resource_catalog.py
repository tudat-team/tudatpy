"""Resource catalog loading and manifest refresh support."""

import csv
import json
import os
import re
import tempfile
from pathlib import Path
from typing import Dict, Iterable, Optional

import requests

RESOURCE_METADATA_URLS = (
    (
        "https://zenodo.org/api/records/21261530/files/manifest-rolling.txt/content",
        "https://zenodo.org/api/records/21261530/files/SHA256SUM-rolling/content",
    ),
    (
        "https://zenodo.org/api/records/22128149/files/manifest-static.txt/content",
        "https://zenodo.org/api/records/22128149/files/SHA256SUM-static/content",
    ),
)
USER_RESOURCE_CATALOG = Path(
    os.environ.get("TUDATPY_RESOURCE_CATALOG", "~/.cache/tudatpy_resources/resource_catalog.csv")
).expanduser()
USER_RESOURCE_HASHES = Path(
    os.environ.get("TUDATPY_RESOURCE_HASHES", "~/.cache/tudatpy_resources/resource_hashes.json")
).expanduser()


def resource_catalog_path() -> Path:
    """Return the user-refreshed catalog or explain how to create it."""
    if not USER_RESOURCE_CATALOG.is_file():
        raise FileNotFoundError(
            f"Resource catalog not found: {USER_RESOURCE_CATALOG}. Create it with "
            "'python -m tudatpy.resource_installer --mode manifest'."
        )
    return USER_RESOURCE_CATALOG


def load_resource_catalog(catalog_path: Optional[Path] = None) -> Dict[str, str]:
    """Load the local resource catalog without making network requests."""
    catalog_path = catalog_path or resource_catalog_path()
    if not catalog_path.is_file():
        raise FileNotFoundError(
            f"Resource catalog not found: {catalog_path}. Create it with "
            "'python -m tudatpy.resource_installer --mode manifest "
            f"--catalog {catalog_path}'."
        )
    with catalog_path.open(newline="", encoding="utf-8") as catalog_file:
        rows = csv.DictReader(catalog_file)
        if rows.fieldnames is None or not {"path", "modified", "url"}.issubset(rows.fieldnames):
            raise ValueError(
                f"Resource catalog {catalog_path} must contain path, modified and url columns."
            )
        catalog = {}
        for row in rows:
            path = row["path"].strip()
            url = row["url"].strip()
            if not path or not url:
                raise ValueError(f"Invalid resource catalog row in {catalog_path}: {row}")
            if path in catalog:
                raise ValueError(f"Duplicate resource path '{path}' in {catalog_path}")
            catalog[path] = url
        return catalog


def _parse_manifest(manifest: str, source: str) -> Iterable[Dict[str, str]]:
    """Yield installable file rows from a Zenodo ``manifest.txt`` file."""
    rows = csv.reader(manifest.splitlines())
    try:
        doi = next(rows)[0]
    except StopIteration:
        raise ValueError(f"Manifest from {source} is empty.")

    record_id = doi.rsplit(".", 1)[-1]
    if not record_id.isdecimal():
        raise ValueError(f"Invalid record DOI in manifest from {source}: {doi}")
    archive = None

    for line_number, row in enumerate(rows, start=2):
        if len(row) != 3:
            raise ValueError(f"Invalid manifest row {line_number} from {source}: {row}")
        path, modified, source_url = (value.strip() for value in row)
        if not source_url:
            # The static manifest currently has a malformed star-catalog
            # archive descriptor. The record exposes the canonical name.
            archive = "star_catalog_biases.tar.gz" if path == ".tar.gz_catalog_biases.tar" else path
            continue
        if not path or not modified:
            raise ValueError(f"Invalid manifest row {line_number} from {source}: {row}")
        if not archive:
            raise ValueError(f"File row {line_number} from {source} has no preceding archive.")

        resource_marker = "/resource/"
        if (
            "raw.githubusercontent.com/tudat-team/tudat-resources/" in source_url
            and resource_marker in source_url
        ):
            path = source_url.split(resource_marker, 1)[1].lstrip("/")
        yield {
            "path": path,
            "modified": modified,
            "source_url": source_url,
            "url": f"https://zenodo.org/api/records/{record_id}/files/{archive}/content",
        }


def _parse_sha256sums(checksums: str, source: str) -> Dict[str, str]:
    """Parse a file in the format emitted by ``sha256sum``."""
    hashes: Dict[str, str] = {}
    for line_number, line in enumerate(checksums.splitlines(), start=1):
        if not line.strip():
            continue
        match = re.fullmatch(r"([0-9a-fA-F]{64})\s+[ *](.+)", line)
        if match is None:
            raise ValueError(f"Invalid SHA256SUM line {line_number} from {source}: {line!r}")
        digest, name = match.groups()
        name = name.strip()
        if not name:
            raise ValueError(f"Empty filename on SHA256SUM line {line_number} from {source}.")
        digest = digest.lower()
        if name in hashes and hashes[name] != digest:
            raise ValueError(f"Conflicting SHA256 hashes for '{name}' in {source}.")
        hashes[name] = digest
    if not hashes:
        raise ValueError(f"SHA256SUM file from {source} is empty.")
    return hashes


def _write_catalog(path: Path, catalog_rows: Dict[str, Dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as catalog_file:
        writer = csv.DictWriter(
            catalog_file,
            fieldnames=["path", "modified", "source_url", "url"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(catalog_rows[item] for item in sorted(catalog_rows))


def update_resource_catalog(
    catalog_path: Path = USER_RESOURCE_CATALOG,
    hash_path: Optional[Path] = None,
) -> int:
    """Fetch the authoritative manifests and checksums and update the local metadata."""
    catalog_path = Path(catalog_path).expanduser()
    hash_path = (
        Path(hash_path).expanduser()
        if hash_path is not None
        else catalog_path.with_name("resource_hashes.json")
    )
    catalog_rows: Dict[str, Dict[str, str]] = {}
    hashes: Dict[str, str] = {}
    for manifest_url, hashes_url in RESOURCE_METADATA_URLS:
        manifest_response = requests.get(manifest_url, timeout=30)
        manifest_response.raise_for_status()
        hashes_response = requests.get(hashes_url, timeout=30)
        hashes_response.raise_for_status()

        for row in _parse_manifest(manifest_response.text, manifest_url):
            path = row["path"]
            if path in catalog_rows:
                raise ValueError(f"Duplicate resource path '{path}' in Zenodo manifests.")
            catalog_rows[path] = row
        for name, digest in _parse_sha256sums(hashes_response.text, hashes_url).items():
            if name in hashes and hashes[name] != digest:
                raise ValueError(f"Conflicting SHA256 hashes for '{name}' in Zenodo records.")
            hashes[name] = digest

    # Write both complete results to temporary files before replacing either
    # local metadata file. Temporary files live beside their targets so that
    # os.replace remains atomic on each filesystem.
    catalog_path.parent.mkdir(parents=True, exist_ok=True)
    hash_path.parent.mkdir(parents=True, exist_ok=True)
    catalog_temp = None
    hash_temp = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=catalog_path.parent, prefix=f".{catalog_path.name}.", delete=False
        ) as temporary:
            catalog_temp = Path(temporary.name)
        with tempfile.NamedTemporaryFile(
            dir=hash_path.parent, prefix=f".{hash_path.name}.", delete=False
        ) as temporary:
            hash_temp = Path(temporary.name)
        _write_catalog(catalog_temp, catalog_rows)
        with hash_temp.open("w", encoding="utf-8") as hash_file:
            json.dump(dict(sorted(hashes.items())), hash_file, indent=2)
            hash_file.write("\n")
        os.replace(catalog_temp, catalog_path)
        catalog_temp = None
        os.replace(hash_temp, hash_path)
        hash_temp = None
    finally:
        if catalog_temp is not None:
            catalog_temp.unlink(missing_ok=True)
        if hash_temp is not None:
            hash_temp.unlink(missing_ok=True)
    return len(catalog_rows)
