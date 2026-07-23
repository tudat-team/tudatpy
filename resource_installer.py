"""Resource installer module for downloading and extracting files."""

import argparse
import os
import tarfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple
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

from resource_data_registry import DATA_REGISTRY

DEFAULT_DEST = Path(os.environ.get("TUDATPY_RESOURCE_DIR", "~/.tudat_resources")).expanduser()
DEFAULT_CACHE = Path(
    os.environ.get("TUDATPY_RESOURCE_CACHE", "~/.cache/tudatpy_resources")
).expanduser()


def _has_progressbar() -> bool:
    return tqdm is not None


def _progress_bar(iterable, total: Optional[int], description: str):
    if _has_progressbar():
        return tqdm(iterable, total=total, unit="B", unit_scale=True, desc=description)
    return iterable


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
    response = requests.get(url, stream=True, timeout=30)
    response.raise_for_status()
    total = int(response.headers.get("Content-Length", 0))
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    with dest_path.open("wb") as fd:
        for chunk in _progress_bar(
            response.iter_content(chunk_size=32_768),
            total=total or None,
            description=dest_path.name,
        ):
            if chunk:
                fd.write(chunk)


def _sha256_file(path: Path) -> str:
    """Return the SHA256 hex digest of a file."""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


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
        got = _sha256_file(dest_path)
        if got.lower() != expected_hash.lower():
            dest_path.unlink(missing_ok=True)
            raise RuntimeError(
                (f"SHA256 mismatch for {dest_path}: " f" expected {expected_hash}, got {got}")
            )
    return True


def _download_tarball(
    url: str, cache_dir: Path, force: bool = False, expected_hash: Optional[str] = None
) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    tar_path = cache_dir / _tarball_cache_name(url)
    if tar_path.exists() and not force:
        return tar_path
    if tar_path.exists():
        tar_path.unlink()

    if pooch is not None:
        pooch.retrieve(url, path=str(tar_path), progressbar=True, retry_if_failed=3)
    else:
        _download_with_requests(url, tar_path)

    if expected_hash:
        got = _sha256_file(tar_path)
        if got.lower() != expected_hash.lower():
            tar_path.unlink(missing_ok=True)
            raise RuntimeError(
                (
                    f"SHA256 mismatch for tarball {tar_path}: "
                    f" expected {expected_hash}, got {got}"
                )
            )
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
    tar_path: Path, targets: List[str], dest_root: Path, force: bool = False
) -> int:
    installed = 0
    with tarfile.open(tar_path, mode="r:*") as tar:
        for target in targets:
            dest_path = dest_root / target
            if dest_path.exists() and not force:
                continue
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            member = _find_tar_member(tar, target)
            extracted = tar.extractfile(member)
            if extracted is None:
                raise RuntimeError(f"Cannot extract '{member.name}' from {tar_path}")
            dest_path.write_bytes(extracted.read())
            installed += 1
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
        installed += _extract_tarball_members(tar_path, targets, dest_path, force=force)

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
        choices=["scratch", "missing", "update", "extra", "list"],
        default="missing",
        help="scratch=overwrite all; missing=download only missing; update=selected subset; extra=download one extra URL; list=show registry entries.",
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
    args = parse_arguments()
    dest_path = Path(args.dest).expanduser()
    cache_dir = Path(args.cache_dir).expanduser()
    hash_map: Optional[Dict[str, str]] = None
    if getattr(args, "hash_file", None):
        p = Path(args.hash_file).expanduser()
        if not p.exists():
            raise FileNotFoundError(f"Hash file not found: {p}")
        hash_map = json.loads(p.read_text())

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
