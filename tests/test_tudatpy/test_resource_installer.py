"""Tests for downloading and caching TudatPy resources."""

import argparse
import hashlib
import io
import json
import tarfile

import pytest
import requests

from tudatpy import resource_installer


class _Response:
    """Minimal streaming response used by the download tests."""

    headers = {"Content-Length": "8"}

    def __enter__(self):
        """Enter the mocked response context."""
        return self

    def __exit__(self, *args):
        """Exit the mocked response context."""
        pass

    def raise_for_status(self):
        """Represent a successful HTTP response."""
        pass

    def iter_content(self, chunk_size):
        """Yield the mocked response body."""
        yield b"resource"


def test_download_retries_and_replaces_atomically(monkeypatch, tmp_path):
    """Retry transient failures and expose only the completed download."""
    attempts = [requests.ConnectionError(), requests.Timeout(), _Response()]

    def get(*args, **kwargs):
        """Return the next failure or response in the attempt sequence."""
        result = attempts.pop(0)
        if isinstance(result, Exception):
            raise result
        return result

    destination = tmp_path / "resource.dat"
    monkeypatch.setattr(resource_installer.requests, "get", get)
    monkeypatch.setattr(resource_installer, "tqdm", None)

    resource_installer._download_with_requests("https://example.invalid/resource.dat", destination)

    assert destination.read_bytes() == b"resource"
    assert not destination.with_name(".resource.dat.part").exists()
    assert not attempts


def test_download_failure_removes_partial_file(monkeypatch, tmp_path):
    """Remove the temporary file when all download attempts fail."""

    class InterruptedResponse(_Response):
        """Mock a response whose stream is interrupted."""

        def iter_content(self, chunk_size):
            """Yield partial content and then interrupt the transfer."""
            yield b"partial"
            raise requests.ConnectionError()

    destination = tmp_path / "resource.dat"
    monkeypatch.setattr(
        resource_installer.requests, "get", lambda *args, **kwargs: InterruptedResponse()
    )
    monkeypatch.setattr(resource_installer, "tqdm", None)

    with pytest.raises(requests.ConnectionError):
        resource_installer._download_with_requests(
            "https://example.invalid/resource.dat", destination
        )

    assert not destination.exists()
    assert not destination.with_name(".resource.dat.part").exists()


def test_verify_sha256_returns_boolean(tmp_path):
    """Report checksum matches and mismatches without raising an exception."""
    resource = tmp_path / "resource.dat"
    resource.write_bytes(b"resource")
    expected_hash = hashlib.sha256(b"resource").hexdigest()

    assert resource_installer._verify_sha256(resource, expected_hash)
    assert not resource_installer._verify_sha256(resource, "0" * 64)


def test_list_catalog_accepts_multiple_key_substrings(monkeypatch, capsys):
    """List the union of resources matching multiple key substrings."""
    monkeypatch.setattr(
        resource_installer,
        "RESOURCE_CATALOG",
        {
            "ephemerides/de430.dat": "https://example.invalid/de430.dat",
            "earth/gravity.dat": "https://example.invalid/gravity.dat",
            "earth/rotation.dat": "https://example.invalid/rotation.dat",
        },
    )

    resource_installer.list_catalog(["de430", "earth"])

    assert capsys.readouterr().out.splitlines() == [
        "earth/gravity.dat",
        "earth/rotation.dat",
        "ephemerides/de430.dat",
    ]


def test_main_uses_files_to_filter_list_mode(monkeypatch, tmp_path):
    """Use the shared files option to select resources in list mode."""
    selected = ["de430", "earth"]
    captured = {}
    monkeypatch.setattr(
        resource_installer,
        "parse_arguments",
        lambda: argparse.Namespace(
            mode="list",
            catalog=None,
            dest=str(tmp_path / "resources"),
            cache_dir=str(tmp_path / "cache"),
            hash_file=None,
            files=selected,
        ),
    )
    monkeypatch.setattr(resource_installer, "load_resource_catalog", lambda path: {})
    monkeypatch.setattr(resource_installer, "_automatic_hash_file", lambda path: None)
    monkeypatch.setattr(resource_installer, "list_catalog", lambda keys: captured.update(keys=keys))

    resource_installer.main()

    assert captured["keys"] == selected


def test_download_tarball_reuses_verified_cache(monkeypatch, tmp_path):
    """Reuse a cached archive only after its checksum has been verified."""
    archive = tmp_path / "resources.tar.gz"
    archive.write_bytes(b"current archive")
    expected_hash = hashlib.sha256(archive.read_bytes()).hexdigest()

    def unexpected_download(*args, **kwargs):
        pytest.fail("A verified archive should not be downloaded again")

    monkeypatch.setattr(resource_installer, "_download_with_requests", unexpected_download)

    result = resource_installer._download_tarball(
        "https://example.invalid/resources.tar.gz",
        tmp_path,
        expected_hash=expected_hash,
    )

    assert result == archive
    assert result.read_bytes() == b"current archive"


def test_download_tarball_replaces_invalid_cache(monkeypatch, tmp_path):
    """Replace a cached archive when it does not match the current checksum."""
    archive = tmp_path / "resources.tar.gz"
    archive.write_bytes(b"outdated archive")
    current_content = b"current archive"
    expected_hash = hashlib.sha256(current_content).hexdigest()
    downloads = []

    def download(url, destination):
        downloads.append((url, destination))
        destination.write_bytes(current_content)

    monkeypatch.setattr(resource_installer, "_download_with_requests", download)

    result = resource_installer._download_tarball(
        "https://example.invalid/resources.tar.gz",
        tmp_path,
        expected_hash=expected_hash,
    )

    assert result.read_bytes() == current_content
    assert downloads == [("https://example.invalid/resources.tar.gz", archive)]


def test_download_tarball_does_not_trust_cache_without_hash(monkeypatch, tmp_path):
    """Redownload a cached archive when no authoritative checksum is available."""
    archive = tmp_path / "resources.tar.gz"
    archive.write_bytes(b"unknown archive")

    def download(url, destination):
        destination.write_bytes(b"fresh archive")

    monkeypatch.setattr(resource_installer, "_download_with_requests", download)

    result = resource_installer._download_tarball(
        "https://example.invalid/resources.tar.gz", tmp_path
    )

    assert result.read_bytes() == b"fresh archive"


def test_main_automatically_loads_manifest_hashes(monkeypatch, tmp_path):
    """Use checksums written beside an explicitly selected catalog."""
    catalog_path = tmp_path / "resource_catalog.csv"
    catalog_path.write_text(
        "path,modified,url\n"
        "resource.dat,2026-08-07T00:00:00+00:00,https://example.invalid/resources.tar.gz\n"
    )
    hash_path = tmp_path / resource_installer.USER_RESOURCE_HASHES.name
    hashes = {"resources.tar.gz": "a" * 64}
    hash_path.write_text(json.dumps(hashes))
    captured = {}

    monkeypatch.setattr(
        resource_installer,
        "parse_arguments",
        lambda: argparse.Namespace(
            mode="missing",
            catalog=str(catalog_path),
            dest=str(tmp_path / "resources"),
            cache_dir=str(tmp_path / "cache"),
            hash_file=None,
            files=None,
            extra_url=None,
            extra_dest=None,
            force=False,
        ),
    )

    def install(files, dest_path, cache_dir, force=False, hashes=None):
        captured["hashes"] = hashes
        return 0

    monkeypatch.setattr(resource_installer, "install_files", install)

    resource_installer.main()

    assert captured["hashes"] == hashes


def test_extract_continues_after_checksum_mismatch(tmp_path):
    """Try every archive member and report failed checksum verification."""
    archive = tmp_path / "resources.tar.gz"
    contents = {"bad.dat": b"corrupt", "good.dat": b"valid"}
    with tarfile.open(archive, "w:gz") as tar:
        for name, content in contents.items():
            member = tarfile.TarInfo(name)
            member.size = len(content)
            tar.addfile(member, io.BytesIO(content))

    hashes = {
        "bad.dat": hashlib.sha256(b"expected").hexdigest(),
        "good.dat": hashlib.sha256(contents["good.dat"]).hexdigest(),
    }
    destination = tmp_path / "destination"

    with pytest.warns(RuntimeWarning, match="bad.dat"):
        installed, failures = resource_installer._extract_tarball_members(
            archive, ["bad.dat", "good.dat"], destination, hashes=hashes
        )

    assert installed == 1
    assert len(failures) == 1
    assert "bad.dat" in failures[0]
    assert not (destination / "bad.dat").exists()
    assert (destination / "good.dat").read_bytes() == b"valid"


def test_extract_replaces_existing_file_with_invalid_checksum(tmp_path):
    """Replace an invalid existing target from the already verified archive."""
    archive = tmp_path / "resources.tar.gz"
    content = b"valid"
    with tarfile.open(archive, "w:gz") as tar:
        member = tarfile.TarInfo("resource.dat")
        member.size = len(content)
        tar.addfile(member, io.BytesIO(content))

    destination = tmp_path / "destination"
    destination.mkdir()
    target = destination / "resource.dat"
    target.write_bytes(b"corrupt")
    hashes = {"resource.dat": hashlib.sha256(content).hexdigest()}

    with pytest.warns(RuntimeWarning, match="resource.dat"):
        installed, failures = resource_installer._extract_tarball_members(
            archive, ["resource.dat"], destination, hashes=hashes
        )

    assert installed == 1
    assert not failures
    assert target.read_bytes() == content


def test_install_attempts_all_tarballs_before_raising(monkeypatch, tmp_path):
    """Aggregate failures only after all tarball groups have been processed."""
    attempted = []

    def download(url, cache_dir, force=False, expected_hash=None):
        return cache_dir / resource_installer._tarball_cache_name(url)

    def extract(tar_path, targets, dest_root, force=False, hashes=None):
        attempted.append(tar_path.name)
        if tar_path.name == "first.tar.gz":
            return 0, ["first.dat: checksum mismatch"]
        return 1, []

    monkeypatch.setattr(resource_installer, "_download_tarball", download)
    monkeypatch.setattr(resource_installer, "_extract_tarball_members", extract)

    files = {
        "first.dat": "https://example.invalid/first.tar.gz",
        "second.dat": "https://example.invalid/second.tar.gz",
    }
    with pytest.raises(resource_installer.ResourceInstallationError) as error:
        resource_installer.install_files(files, tmp_path / "dest", tmp_path / "cache")

    assert attempted == ["first.tar.gz", "second.tar.gz"]
    assert error.value.failures == ["first.dat: checksum mismatch"]
