"""Tests for downloading and caching TudatPy resources."""

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
