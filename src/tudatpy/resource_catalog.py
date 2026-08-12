"""Resource catalog loading and manifest refresh support."""

import csv
import os
from importlib import resources
from importlib.resources.abc import Traversable
from pathlib import Path
from typing import Dict, Iterable, Optional, Union

import requests

MANIFEST_URLS = (
    "https://zenodo.org/api/records/21277280/files/manifest.txt/content",
    "https://zenodo.org/api/records/21261530/files/manifest.txt/content",
)
USER_RESOURCE_CATALOG = Path(
    os.environ.get("TUDATPY_RESOURCE_CATALOG", "~/.cache/tudatpy_resources/resource_catalog.csv")
).expanduser()

CatalogPath = Union[Path, Traversable]


def packaged_resource_catalog() -> Traversable:
    """Return the catalog bundled with the installed ``tudatpy`` package."""
    return resources.files("tudatpy").joinpath("resource_catalog.csv")


def resource_catalog_path() -> CatalogPath:
    """Return a user-refreshed catalog when available, else the bundled one."""
    return USER_RESOURCE_CATALOG if USER_RESOURCE_CATALOG.exists() else packaged_resource_catalog()


def load_resource_catalog(catalog_path: Optional[CatalogPath] = None) -> Dict[str, str]:
    """Load the local resource catalog without making network requests."""
    catalog_path = catalog_path or resource_catalog_path()
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


def update_resource_catalog(catalog_path: Path = USER_RESOURCE_CATALOG) -> int:
    """Fetch both authoritative manifests and write a local catalog override."""
    catalog_rows: Dict[str, Dict[str, str]] = {}
    for manifest_url in MANIFEST_URLS:
        response = requests.get(manifest_url, timeout=30)
        response.raise_for_status()
        for row in _parse_manifest(response.text, manifest_url):
            path = row["path"]
            if path in catalog_rows:
                raise ValueError(f"Duplicate resource path '{path}' in Zenodo manifests.")
            catalog_rows[path] = row

    catalog_path.parent.mkdir(parents=True, exist_ok=True)
    with catalog_path.open("w", newline="", encoding="utf-8") as catalog_file:
        writer = csv.DictWriter(
            catalog_file,
            fieldnames=["path", "modified", "source_url", "url"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(catalog_rows[path] for path in sorted(catalog_rows))
    return len(catalog_rows)
