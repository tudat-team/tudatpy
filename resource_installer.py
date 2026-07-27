"""Resource installer module for downloading and extracting files."""

import argparse
import csv
import os
import tarfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import hashlib
import json
from urllib.parse import urlparse

import requests

try:
    import pooch
except ImportError:  # pragma: no cover
    pooch = None

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover
    tqdm = None

DEFAULT_DEST = Path(os.environ.get("TUDATPY_RESOURCE_DIR", "~/.tudat_resources")).expanduser()
DEFAULT_CACHE = Path(
    os.environ.get("TUDATPY_RESOURCE_CACHE", "~/.cache/tudatpy_resources")
).expanduser()
RESOURCE_CATALOG = Path(__file__).with_name("resource_catalog.csv")
MANIFEST_URLS = (
    "https://zenodo.org/api/records/21277280/files/manifest.txt/content",
    "https://zenodo.org/api/records/21261530/files/manifest.txt/content",
)


def load_resource_catalog(catalog_path: Path = RESOURCE_CATALOG) -> Dict[str, str]:
    """Load the local, offline resource catalog.

    The catalog contains the manifest path, timestamp and source URL, plus the
    Zenodo tarball URL used for installation. The latter avoids relying on
    upstream source URLs that may no longer be available.
    """
    try:
        with catalog_path.open(newline="", encoding="utf-8") as catalog_file:
            rows = csv.DictReader(catalog_file)
            if rows.fieldnames is None or not {"path", "modified", "url"}.issubset(rows.fieldnames):
                raise ValueError(
                    f"Resource catalog {catalog_path} must contain path, modified and url columns."
                )
            registry = {}
            for row in rows:
                path = row["path"].strip()
                url = row["url"].strip()
                if not path or not url:
                    raise ValueError(f"Invalid resource catalog row in {catalog_path}: {row}")
                if path in registry:
                    raise ValueError(f"Duplicate resource path '{path}' in {catalog_path}")
                registry[path] = url
            return registry
    except FileNotFoundError as error:
        raise FileNotFoundError(f"Resource catalog not found: {catalog_path}") from error


def _parse_manifest(manifest: str, source: str) -> Iterable[Dict[str, str]]:
    """Yield file rows from a Zenodo ``manifest.txt`` file.

    The first line is the record DOI. Tarball descriptor rows determine the
    download URL for the file rows that follow them.
    """
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
        path, modified, url = (value.strip() for value in row)
        if not url:
            # The static manifest currently has a malformed star-catalog
            # archive descriptor. The record itself exposes the canonical name.
            archive = "star_catalog_biases.tar.gz" if path == ".tar.gz_catalog_biases.tar" else path
            continue
        if not path or not modified:
            raise ValueError(f"Invalid manifest row {line_number} from {source}: {row}")
        if not archive:
            raise ValueError(f"File row {line_number} from {source} has no preceding archive.")
        # The resource repository URL is authoritative when the manifest path
        # is malformed (currently the star-catalog-biases entries have this
        # issue). This avoids preserving a bad install path in the catalog.
        resource_marker = "/resource/"
        if (
            "raw.githubusercontent.com/tudat-team/tudat-resources/" in url
            and resource_marker in url
        ):
            path = url.split(resource_marker, 1)[1].lstrip("/")
        yield {
            "path": path,
            "modified": modified,
            "source_url": url,
            "url": f"https://zenodo.org/api/records/{record_id}/files/{archive}/content",
        }


def update_resource_catalog(catalog_path: Path = RESOURCE_CATALOG) -> int:
    """Fetch both authoritative Zenodo manifests and replace the local catalog."""
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


# Keep module import usable before the catalog has been created for the first
# time. All normal CLI modes load it explicitly in ``main``.
DATA_REGISTRY = load_resource_catalog() if RESOURCE_CATALOG.exists() else {}


def _has_progressbar() -> bool:
    return tqdm is not None


def find_in_registry(search_string: str) -> Dict[str, str]:
    """Return registry entries whose keys contain ``search_string``.

    Parameters
    ----------
    search_string:
        Substring to search for in registry keys.

    Returns
    -------
    Dict[str, str]
        Mapping of matching registry keys to their URLs.
    """
    return {key: url for key, url in DATA_REGISTRY.items() if search_string in key}


def resolve_keys(keys: Optional[List[str]] = None, search: Optional[str] = None) -> Dict[str, str]:
    """Resolve registry keys using search string or key matching.

    Parameters
    ----------
    keys:
        List of registry keys or partial key strings to match.
    search:
        Substring to search for in registry keys.

    Returns
    -------
    Dict[str, str]
        Mapping of matched registry keys to their URLs.
    """
    resolved: Dict[str, str] = {}
    if search:
        resolved.update(find_in_registry(search))

    if keys:
        for key in keys:
            if key in DATA_REGISTRY:
                resolved[key] = DATA_REGISTRY[key]
                continue

            matches = {
                registry_key: url
                for registry_key, url in DATA_REGISTRY.items()
                if registry_key.startswith(key) or key in registry_key
            }
            if not matches:
                raise ValueError(f"No registry entries match '{key}'")
            resolved.update(matches)

    return resolved


def _tarball_cache_name(url: str) -> str:
    parsed = urlparse(url)
    name = Path(parsed.path).name
    if name == "content":
        parent = Path(parsed.path).parent
        if parent.name.endswith((".tar.gz", ".tgz", ".tar")):
            return parent.name
    return name


def _is_tarball_url(url: str) -> bool:
    parsed = urlparse(url)
    name = Path(parsed.path).name
    if name.endswith((".tar.gz", ".tgz", ".tar")):
        return True
    if name == "content" and Path(parsed.path).parent.name.endswith((".tar.gz", ".tgz", ".tar")):
        return True
    return False


def _split_tarball_files(files: Dict[str, str]) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    regular: Dict[str, str] = {}
    tarball_groups: Dict[str, List[str]] = {}
    for path, url in files.items():
        if _is_tarball_url(url):
            tarball_groups.setdefault(url, []).append(path)
        else:
            regular[path] = url
    return regular, tarball_groups


def _download_with_requests(url: str, dest_path: Path) -> None:
    print(f"Downloading {dest_path.name}")
    response = requests.get(url, stream=True, timeout=30)
    response.raise_for_status()
    total = int(response.headers.get("Content-Length", 0))
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    progress = None
    if _has_progressbar():
        progress = tqdm(
            total=total or None,
            unit="B",
            unit_scale=True,
            desc=f"Downloading {dest_path.name}",
        )
    try:
        with dest_path.open("wb") as fd:
            for chunk in response.iter_content(chunk_size=32_768):
                if chunk:
                    fd.write(chunk)
                    if progress is not None:
                        progress.update(len(chunk))
    finally:
        if progress is not None:
            progress.close()


def _sha256_file(path: Path) -> str:
    """Return the SHA256 hex digest of a file."""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _verify_sha256(path: Path, expected_hash: str) -> None:
    """Verify a file against an expected SHA256 digest."""
    got = _sha256_file(path)
    if got.lower() != expected_hash.lower():
        raise RuntimeError(f"SHA256 mismatch for {path}: expected {expected_hash}, got {got}")


def _download_file(
    url: str, dest_path: Path, force: bool = False, expected_hash: Optional[str] = None
) -> bool:
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    if dest_path.exists():
        if not force:
            return False
        dest_path.unlink()

    if pooch is not None:
        pooch.retrieve(url, path=str(dest_path), progressbar=True, retry_if_failed=3)
    else:
        _download_with_requests(url, dest_path)

    if expected_hash:
        try:
            _verify_sha256(dest_path, expected_hash)
        except RuntimeError:
            dest_path.unlink(missing_ok=True)
            raise
        print(f"Verified SHA256 for {dest_path}")
    return True


def _download_tarball(
    url: str, cache_dir: Path, force: bool = False, expected_hash: Optional[str] = None
) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    tar_path = cache_dir / _tarball_cache_name(url)
    if tar_path.exists() and not force:
        print(f"Using cached archive {tar_path.name}")
        return tar_path
    if tar_path.exists():
        tar_path.unlink()

    if pooch is not None:
        pooch.retrieve(url, path=str(tar_path), progressbar=True, retry_if_failed=3)
    else:
        _download_with_requests(url, tar_path)

    if expected_hash:
        try:
            _verify_sha256(tar_path, expected_hash)
        except RuntimeError:
            tar_path.unlink(missing_ok=True)
            raise
        print(f"Verified SHA256 for archive {tar_path.name}")
    return tar_path


def _find_tar_member(tar: tarfile.TarFile, target: str) -> tarfile.TarInfo:
    normalized_target = Path(target).as_posix()
    fallback = None
    for member in tar.getmembers():
        if member.isdir():
            continue
        member_path = Path(member.name).as_posix()
        if member_path == normalized_target:
            return member
        if Path(member.name).name == Path(target).name and fallback is None:
            fallback = member
    if fallback is not None:
        return fallback
    raise KeyError(f"File '{target}' not found in tarball")


def _extract_tarball_members(
    tar_path: Path,
    targets: List[str],
    dest_root: Path,
    force: bool = False,
    hashes: Optional[Dict[str, str]] = None,
) -> int:
    installed = 0
    verified = 0
    print(f"Extracting {len(targets)} resources from {tar_path.name}")
    with tarfile.open(tar_path, mode="r:*") as tar:
        for target in targets:
            dest_path = dest_root / target
            if dest_path.exists() and not force:
                if hashes and (expected_hash := hashes.get(target)):
                    _verify_sha256(dest_path, expected_hash)
                    verified += 1
                continue
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            member = _find_tar_member(tar, target)
            extracted = tar.extractfile(member)
            if extracted is None:
                raise RuntimeError(f"Cannot extract '{member.name}' from {tar_path}")
            dest_path.write_bytes(extracted.read())
            if hashes and (expected_hash := hashes.get(target)):
                try:
                    _verify_sha256(dest_path, expected_hash)
                except RuntimeError:
                    dest_path.unlink(missing_ok=True)
                    raise
                verified += 1
            installed += 1
    if verified:
        print(f"Verified SHA256 for {verified} resources from {tar_path.name}")
    return installed


def install_files(
    files: Dict[str, str],
    dest_path: Path,
    cache_dir: Path,
    force: bool = False,
    hashes: Optional[Dict[str, str]] = None,
) -> int:
    """Install files from URLs and tarballs to destination directory."""
    dest_path = dest_path.expanduser()
    cache_dir = cache_dir.expanduser()
    dest_path.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    regular_files, tarball_groups = _split_tarball_files(files)
    installed = 0

    for relative_path, url in regular_files.items():
        target_path = dest_path / relative_path
        expected = None
        if hashes:
            expected = hashes.get(relative_path) or hashes.get(url)
        if _download_file(url, target_path, force=force, expected_hash=expected):
            installed += 1

    for tarball_url, targets in tarball_groups.items():
        expected = None
        if hashes:
            expected = hashes.get(tarball_url) or hashes.get(_tarball_cache_name(tarball_url))
        tar_path = _download_tarball(tarball_url, cache_dir, force=force, expected_hash=expected)
        installed += _extract_tarball_members(
            tar_path, targets, dest_path, force=force, hashes=hashes
        )

    return installed


def download_extra_file(
    url: str,
    dest_path: Path,
    cache_dir: Path,
    force: bool = False,
    hashes: Optional[Dict[str, str]] = None,
) -> Path:
    """Download an extra file to the destination path.

    If dest_path is a directory, the filename is derived from the URL.
    The destination directory is created if it does not exist.
    """
    dest_path = dest_path.expanduser()
    if dest_path.is_dir():
        dest_path = dest_path / Path(url).name
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    expected = None
    if hashes:
        expected = hashes.get(url) or hashes.get(Path(url).name)
    if _download_file(url, dest_path, force=force, expected_hash=expected):
        return dest_path
    return dest_path


def list_registry(search: Optional[str] = None) -> None:
    """Print registry entries, optionally filtered by search string."""
    entries = DATA_REGISTRY
    if search:
        entries = find_in_registry(search)
    for key in sorted(entries):
        print(key)


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments for the resource installer."""
    parser = argparse.ArgumentParser(
        description="Install or update TudatPy common resource files with cache and progress support."
    )
    parser.add_argument(
        "--mode",
        choices=["scratch", "missing", "update", "extra", "list", "manifest"],
        default="missing",
        help=(
            "scratch=overwrite all; missing=download only missing; update=selected subset; "
            "extra=download one extra URL; list=show the local catalog; "
            "manifest=refresh the local catalog from Zenodo manifests."
        ),
    )
    parser.add_argument(
        "--dest", default=str(DEFAULT_DEST), help="Destination root for resource files."
    )
    parser.add_argument(
        "--cache-dir",
        default=str(DEFAULT_CACHE),
        help="Cache directory for downloads and tarballs.",
    )
    parser.add_argument(
        "--files", nargs="+", help="Exact registry keys or prefixes for update mode."
    )
    parser.add_argument(
        "--search", help="Substring search to select registry files for update or list."
    )
    parser.add_argument("--extra-url", help="URL of an extra file to download in extra mode.")
    parser.add_argument(
        "--extra-dest", help="Relative or absolute destination path for extra mode."
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing files.")
    parser.add_argument(
        "--list-search", help="Search string to filter registry output in list mode."
    )
    parser.add_argument(
        "--hash-file",
        help="Path to a JSON file mapping registry keys or URLs to SHA256 hex digests.",
    )
    return parser.parse_args()


def main() -> None:
    """Run the resource installer command-line interface."""
    global DATA_REGISTRY
    args = parse_arguments()

    if args.mode == "manifest":
        updated = update_resource_catalog()
        print(f"Updated {RESOURCE_CATALOG} with {updated} resources from Zenodo manifests")
        return

    DATA_REGISTRY = load_resource_catalog()
    dest_path = Path(args.dest).expanduser()
    cache_dir = Path(args.cache_dir).expanduser()
    hash_map: Optional[Dict[str, str]] = None
    if getattr(args, "hash_file", None):
        p = Path(args.hash_file).expanduser()
        if not p.exists():
            raise FileNotFoundError(f"Hash file not found: {p}")
        hash_map = json.loads(p.read_text())
        print(f"Loaded {len(hash_map)} SHA256 hashes from {p}")

    if args.mode == "list":
        list_registry(args.list_search or args.search)
        return

    if args.mode == "scratch":
        installed = install_files(DATA_REGISTRY, dest_path, cache_dir, force=True, hashes=hash_map)
        print(f"Installed {installed} resources to {dest_path}")
        return

    if args.mode == "missing":
        installed = install_files(DATA_REGISTRY, dest_path, cache_dir, force=False, hashes=hash_map)
        print(f"Installed {installed} missing resources to {dest_path}")
        return

    if args.mode == "update":
        if not args.files and not args.search:
            raise ValueError("Update mode requires --files or --search.")
        files = resolve_keys(args.files, args.search)
        installed = install_files(files, dest_path, cache_dir, force=True, hashes=hash_map)
        print(f"Updated {installed} selected resources to {dest_path}")
        return

    if args.mode == "extra":
        if not args.extra_url:
            raise ValueError("Extra mode requires --extra-url.")
        extra_dest = Path(args.extra_dest) if args.extra_dest else dest_path
        if not extra_dest.is_absolute():
            extra_dest = dest_path / extra_dest
        downloaded = download_extra_file(
            args.extra_url, extra_dest, cache_dir, force=args.force, hashes=hash_map
        )
        print(f"Downloaded extra resource to {downloaded}")
        return

    raise RuntimeError(f"Unknown mode: {args.mode}")


if __name__ == "__main__":
    main()
