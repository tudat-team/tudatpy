"""Resource installer module for downloading and extracting files."""

import argparse
import os
import tarfile
import warnings
from importlib import resources
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import hashlib
import json
from urllib.parse import urlparse

import requests

from tudatpy.resource_catalog import (
    USER_RESOURCE_CATALOG,
    USER_RESOURCE_HASHES,
    load_resource_catalog,
    update_resource_catalog,
)

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover
    tqdm = None

DEFAULT_DEST = Path(os.environ.get("TUDATPY_RESOURCE_DIR", "~/.tudat/resource")).expanduser()
DEFAULT_CACHE = Path(
    os.environ.get("TUDATPY_RESOURCE_CACHE", "~/.cache/tudatpy_resources")
).expanduser()
# Loading is deferred until after manifest mode has had the opportunity to
# create the local catalog on a fresh installation.
RESOURCE_CATALOG: Dict[str, str] = {}


class ChecksumMismatchError(RuntimeError):
    """Raised when a file does not match its expected SHA-256 digest."""


class ResourceInstallationError(RuntimeError):
    """Raised after all resources have been attempted and some have failed."""

    def __init__(self, failures: List[str]):
        self.failures = failures
        details = "\n".join(f"- {failure}" for failure in failures)
        super().__init__(f"Failed to install {len(failures)} resource(s):\n{details}")


def _has_progressbar() -> bool:
    return tqdm is not None


def find_in_catalog(search_string: str) -> Dict[str, str]:
    """Return catalog entries whose keys contain ``search_string``.

    Parameters
    ----------
    search_string:
        Substring to search for in catalog keys.

    Returns
    -------
    Dict[str, str]
        Mapping of matching catalog keys to their URLs.
    """
    return {key: url for key, url in RESOURCE_CATALOG.items() if search_string in key}


def resolve_catalog_keys(keys: List[str]) -> Dict[str, str]:
    """Resolve catalog keys by substring matching.

    Parameters
    ----------
    keys:
        Substrings to search for in catalog keys.

    Returns
    -------
    Dict[str, str]
        Mapping of matched catalog keys to their URLs.
    """
    resolved: Dict[str, str] = {}
    for key in keys:
        resolved.update(find_in_catalog(key))

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


def _download_with_requests(url: str, dest_path: Path, retries: int = 3) -> None:
    """Download ``url`` atomically, retrying failed transfers."""
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = dest_path.with_name(f".{dest_path.name}.part")

    for attempt in range(1, retries + 1):
        print(f"Downloading {dest_path.name} (attempt {attempt}/{retries})")
        progress = None
        try:
            with requests.get(url, stream=True, timeout=30) as response:
                response.raise_for_status()
                total = int(response.headers.get("Content-Length", 0))
                if _has_progressbar():
                    progress = tqdm(
                        total=total or None,
                        unit="B",
                        unit_scale=True,
                        desc=f"Downloading {dest_path.name}",
                    )
                with temporary_path.open("wb") as fd:
                    for chunk in response.iter_content(chunk_size=32_768):
                        if chunk:
                            fd.write(chunk)
                            if progress is not None:
                                progress.update(len(chunk))
            temporary_path.replace(dest_path)
            return
        except (requests.RequestException, OSError):
            temporary_path.unlink(missing_ok=True)
            if attempt == retries:
                raise
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


def _verify_sha256(path: Path, expected_hash: str) -> bool:
    """Return whether a file matches its expected SHA-256 digest."""
    return _sha256_file(path).lower() == expected_hash.lower()


def _checksum_mismatch_message(path: Path, expected_hash: str) -> str:
    """Return a diagnostic for a file that failed checksum verification."""
    return f"SHA256 mismatch for {path}: expected {expected_hash}"


def _download_file(
    url: str, dest_path: Path, force: bool = False, expected_hash: Optional[str] = None
) -> bool:
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    if dest_path.exists():
        if not force:
            return False
        dest_path.unlink()

    _download_with_requests(url, dest_path)

    if expected_hash and not _verify_sha256(dest_path, expected_hash):
        dest_path.unlink(missing_ok=True)
        raise ChecksumMismatchError(_checksum_mismatch_message(dest_path, expected_hash))
    if expected_hash:
        print(f"Verified SHA256 for {dest_path}")
    return True


def _download_tarball(
    url: str, cache_dir: Path, force: bool = False, expected_hash: Optional[str] = None
) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    tar_path = cache_dir / _tarball_cache_name(url)
    if tar_path.exists() and not force:
        if expected_hash:
            if not _verify_sha256(tar_path, expected_hash):
                print(f"Discarding outdated cached archive {tar_path.name}")
                tar_path.unlink()
            else:
                print(f"Using verified cached archive {tar_path.name}")
                return tar_path
        else:
            # A cached archive cannot be considered current without an
            # authoritative checksum against which to verify it.
            tar_path.unlink()
    if tar_path.exists():
        tar_path.unlink()

    _download_with_requests(url, tar_path)

    if expected_hash and not _verify_sha256(tar_path, expected_hash):
        tar_path.unlink(missing_ok=True)
        raise ChecksumMismatchError(_checksum_mismatch_message(tar_path, expected_hash))
    if expected_hash:
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
) -> Tuple[int, List[str]]:
    installed = 0
    verified = 0
    failures: List[str] = []
    print(f"Extracting {len(targets)} resources from {tar_path.name}")
    with tarfile.open(tar_path, mode="r:*") as tar:
        for target in targets:
            dest_path = dest_root / target
            if dest_path.exists() and not force:
                if hashes and (expected_hash := hashes.get(target)):
                    if not _verify_sha256(dest_path, expected_hash):
                        warning = _checksum_mismatch_message(dest_path, expected_hash)
                        warnings.warn(warning, RuntimeWarning, stacklevel=2)
                        dest_path.unlink(missing_ok=True)
                    else:
                        verified += 1
                        continue
                else:
                    continue
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            try:
                member = _find_tar_member(tar, target)
                extracted = tar.extractfile(member)
                if extracted is None:
                    raise RuntimeError(f"Cannot extract '{member.name}' from {tar_path}")
                dest_path.write_bytes(extracted.read())
            except (KeyError, OSError, RuntimeError, tarfile.TarError) as error:
                dest_path.unlink(missing_ok=True)
                failure = f"{target}: {error}"
                warnings.warn(failure, RuntimeWarning, stacklevel=2)
                failures.append(failure)
                continue
            if hashes and (expected_hash := hashes.get(target)):
                if not _verify_sha256(dest_path, expected_hash):
                    dest_path.unlink(missing_ok=True)
                    failure = f"{target}: {_checksum_mismatch_message(dest_path, expected_hash)}"
                    warnings.warn(failure, RuntimeWarning, stacklevel=2)
                    failures.append(failure)
                    continue
                verified += 1
            installed += 1
    if verified:
        print(f"Verified SHA256 for {verified} resources from {tar_path.name}")
    return installed, failures


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
    failures: List[str] = []

    for relative_path, url in regular_files.items():
        target_path = dest_path / relative_path
        expected = None
        if hashes:
            expected = hashes.get(relative_path) or hashes.get(url)
        try:
            if _download_file(url, target_path, force=force, expected_hash=expected):
                installed += 1
        except ChecksumMismatchError as error:
            failure = f"{relative_path}: {error}"
            warnings.warn(failure, RuntimeWarning, stacklevel=2)
            failures.append(failure)

    for tarball_url, targets in tarball_groups.items():
        expected = None
        if hashes:
            expected = hashes.get(tarball_url) or hashes.get(_tarball_cache_name(tarball_url))
        try:
            tar_path = _download_tarball(
                tarball_url, cache_dir, force=force, expected_hash=expected
            )
            extracted, extraction_failures = _extract_tarball_members(
                tar_path, targets, dest_path, force=force, hashes=hashes
            )
        except (ChecksumMismatchError, OSError, tarfile.TarError) as error:
            failure = f"{_tarball_cache_name(tarball_url)}: {error}"
            warnings.warn(failure, RuntimeWarning, stacklevel=2)
            failures.append(failure)
            continue
        installed += extracted
        failures.extend(extraction_failures)

    if failures:
        raise ResourceInstallationError(failures)

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


def list_catalog(keys: Optional[List[str]] = None) -> None:
    """Print catalog entries, optionally filtered by key substrings."""
    entries = resolve_catalog_keys(keys) if keys else RESOURCE_CATALOG
    for key in sorted(entries):
        print(key)


def _load_hash_file(hash_file: str) -> Tuple[Dict[str, str], str]:
    """Load a user-supplied hash file or the hash file bundled with TudatPy."""
    path = Path(hash_file).expanduser()
    if path.exists():
        return json.loads(path.read_text()), str(path)

    if path.name == hash_file:
        packaged_hashes = resources.files("tudatpy").joinpath(hash_file)
        if packaged_hashes.is_file():
            with packaged_hashes.open(encoding="utf-8") as hash_stream:
                return json.load(hash_stream), f"package resource {hash_file}"

    raise FileNotFoundError(f"Hash file not found: {path}")


def _automatic_hash_file(catalog_path: Optional[Path]) -> Optional[Path]:
    """Return the checksum file associated with the selected local catalog."""
    if catalog_path is not None:
        candidate = catalog_path.with_name(USER_RESOURCE_HASHES.name)
    elif USER_RESOURCE_CATALOG.exists():
        candidate = USER_RESOURCE_HASHES
    else:
        return None
    return candidate if candidate.exists() else None


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
        "--files", nargs="+", help="Substrings selecting catalog files for update mode."
    )
    parser.add_argument("--extra-url", help="URL of an extra file to download in extra mode.")
    parser.add_argument(
        "--extra-dest", help="Relative or absolute destination path for extra mode."
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing files.")
    parser.add_argument(
        "--list-search",
        nargs="+",
        help="Substrings used to filter catalog output in list mode.",
    )
    parser.add_argument(
        "--hash-file",
        help="Path to a JSON file mapping catalog keys or URLs to SHA256 hex digests.",
    )
    parser.add_argument(
        "--catalog",
        help=(
            "Path to a resource catalog. With --mode manifest, this is the output path; "
            "otherwise it overrides the user or packaged catalog."
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Run the resource installer command-line interface."""
    global RESOURCE_CATALOG
    args = parse_arguments()
    catalog_path = Path(args.catalog).expanduser() if args.catalog else None

    if args.mode == "manifest":
        output_path = catalog_path or USER_RESOURCE_CATALOG
        hash_path = (
            output_path.with_name(USER_RESOURCE_HASHES.name)
            if catalog_path is not None
            else USER_RESOURCE_HASHES
        )
        updated = update_resource_catalog(output_path, hash_path)
        print(
            f"Updated {output_path} with {updated} resources and wrote SHA256 hashes "
            f"to {hash_path}"
        )
        return

    RESOURCE_CATALOG = load_resource_catalog(catalog_path)
    dest_path = Path(args.dest).expanduser()
    cache_dir = Path(args.cache_dir).expanduser()
    hash_map: Optional[Dict[str, str]] = None
    if getattr(args, "hash_file", None):
        hash_map, hash_source = _load_hash_file(args.hash_file)
        print(f"Loaded {len(hash_map)} SHA256 hashes from {hash_source}")
    elif automatic_hash_file := _automatic_hash_file(catalog_path):
        hash_map, hash_source = _load_hash_file(str(automatic_hash_file))
        print(f"Loaded {len(hash_map)} SHA256 hashes from {hash_source}")

    if args.mode == "list":
        list_catalog(args.list_search)
        return

    if args.mode == "scratch":
        installed = install_files(
            RESOURCE_CATALOG, dest_path, cache_dir, force=True, hashes=hash_map
        )
        print(f"Installed {installed} resources to {dest_path}")
        return

    if args.mode == "missing":
        installed = install_files(
            RESOURCE_CATALOG, dest_path, cache_dir, force=False, hashes=hash_map
        )
        print(f"Installed {installed} missing resources to {dest_path}")
        return

    if args.mode == "update":
        if not args.files:
            raise ValueError("Update mode requires --files.")
        files = resolve_catalog_keys(args.files)
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
