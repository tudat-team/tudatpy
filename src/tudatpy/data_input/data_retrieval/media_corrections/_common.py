"""
Shared download infrastructure for media-correction ancillary data.

This module holds the functionality common to the IONEX and VMF downloaders:
the HTTP/curl download layer, the generic daily-file synchronisation
orchestrator, the shared return type, and the common exceptions and helpers.
IONEX-specific code lives in :mod:`.ionex`, VMF-specific code in :mod:`.vmf`.
"""

from __future__ import annotations

import gzip
import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from datetime import date, datetime, timedelta
from pathlib import Path
from typing import Callable, Sequence

try:
    import requests as _requests
except ImportError:
    _requests = None  # type: ignore[assignment]

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class AncillaryDownloadError(Exception):
    """Base exception for download failures."""


class AuthenticationError(AncillaryDownloadError):
    """Raised when authentication credentials are missing or invalid.

    IONEX downloads from NASA CDDIS require Earthdata credentials stored in
    a ``~/.netrc`` file. See the ``download_ionex`` docstring for setup
    instructions.
    """


# ---------------------------------------------------------------------------
# Return type
# ---------------------------------------------------------------------------


@dataclass
class DownloadResult:
    """Result of a download operation.

    Attributes
    ----------
    downloaded : list[Path]
        Files newly fetched in this call.
    existing : list[Path]
        Files that already existed locally and were not re-downloaded.
    failed : list[str]
        Human-readable failure descriptions for days where no file could
        be obtained.
    """

    downloaded: list[Path] = field(default_factory=list)
    existing: list[Path] = field(default_factory=list)
    failed: list[str] = field(default_factory=list)

    @property
    def all_files(self) -> list[Path]:
        """**read-only**

        All available files (downloaded + existing), sorted by name.

        :type: list[Path]
        """
        return sorted(set(self.downloaded + self.existing))

    @property
    def success(self) -> bool:
        """**read-only**

        True when at least one file is available and nothing failed.

        :type: bool
        """
        return len(self.all_files) > 0 and len(self.failed) == 0


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_HTTP_CHUNK_SIZE = 1024 * 1024
_HTTP_MAX_ATTEMPTS = 3


@dataclass(frozen=True)
class _RemoteCandidate:
    local_name: str
    url: str
    download_name: str | None = None
    gunzip: bool = False


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _safe_unlink(path: Path) -> None:
    try:
        path.unlink(missing_ok=True)
    except Exception:
        pass


def _format_doy(doy: int) -> str:
    return f"{doy:03d}"


def _date_range_to_year_doys(
    start: datetime,
    end: datetime,
    padding: int = 0,
) -> list[tuple[int, str]]:
    """Convert a UTC datetime range to a list of ``(year, doy_str)`` tuples."""
    start_date = start.date() - timedelta(days=padding)
    end_date = end.date() + timedelta(days=padding)
    result: list[tuple[int, str]] = []
    day = start_date
    while day <= end_date:
        doy = (day - date(day.year, 1, 1)).days + 1
        result.append((day.year, _format_doy(doy)))
        day += timedelta(days=1)
    return result


def _resolve_path(p: Path | str) -> Path:
    return Path(p).expanduser().resolve()


# ---------------------------------------------------------------------------
# HTTP download layer
# ---------------------------------------------------------------------------


def _http_download(
    url: str,
    dest: Path,
    *,
    netrc_path: Path | None = None,
    max_attempts: int = _HTTP_MAX_ATTEMPTS,
) -> tuple[bool, int | None, str | None]:
    """Download *url* to *dest*.  Returns ``(ok, status_code, error_string)``."""
    # For CDDIS, prefer curl (netrc + cookie handling).
    if "cddis.nasa.gov" in url:
        ok, status, err = _curl_download(
            url, dest, netrc_path=netrc_path, max_attempts=max_attempts
        )
        if ok or err != "curl-not-available":
            return ok, status, err

    return _requests_download(url, dest, max_attempts=max_attempts)


def _requests_download(
    url: str,
    dest: Path,
    *,
    max_attempts: int = _HTTP_MAX_ATTEMPTS,
) -> tuple[bool, int | None, str | None]:
    if _requests is None:
        log.warning("'requests' library not available; cannot fetch: %s", url)
        return False, None, "requests-not-available"

    tmp = dest.with_suffix(dest.suffix + ".part")
    last_error: str | None = None
    last_status: int | None = None

    for _ in range(max(1, max_attempts)):
        _safe_unlink(tmp)
        try:
            with _requests.get(url, stream=True, timeout=(10, 120)) as resp:
                status = resp.status_code
                last_status = status
                if status == 404:
                    return False, 404, None
                if status != 200:
                    if 300 <= status < 400 or status >= 500 or status == 429:
                        continue
                    return False, status, None

                content_length = resp.headers.get("content-length")
                expected = (
                    int(content_length) if content_length and content_length.isdigit() else None
                )
                written = 0
                with open(tmp, "wb") as fh:
                    for chunk in resp.iter_content(chunk_size=_HTTP_CHUNK_SIZE):
                        if chunk:
                            fh.write(chunk)
                            written += len(chunk)

                if expected is not None and written != expected:
                    raise IOError(f"Incomplete download: expected {expected} bytes, got {written}.")
                tmp.replace(dest)
                return True, 200, None
        except Exception as exc:
            last_error = str(exc)
            _safe_unlink(tmp)

    return False, last_status, last_error


def _curl_download(
    url: str,
    dest: Path,
    *,
    netrc_path: Path | None = None,
    max_attempts: int = _HTTP_MAX_ATTEMPTS,
) -> tuple[bool, int | None, str | None]:
    curl_bin = shutil.which("curl")
    if not curl_bin:
        return False, None, "curl-not-available"

    if netrc_path is None:
        netrc_path = Path.home() / ".netrc"
    else:
        netrc_path = _resolve_path(netrc_path)

    tmp = dest.with_suffix(dest.suffix + ".part")
    last_error: str | None = None
    last_status: int | None = None

    for _ in range(max(1, max_attempts)):
        _safe_unlink(tmp)
        cookie_path: Path | None = None
        try:
            with tempfile.NamedTemporaryFile(
                prefix="ancillary_cookie_",
                suffix=".txt",
                dir=tmp.parent if tmp.parent.exists() else None,
                delete=False,
            ) as cookie_fh:
                cookie_path = Path(cookie_fh.name)

            cmd: list[str] = [
                curl_bin,
                "-sS",
                "-L",
                "--connect-timeout",
                "10",
                "--max-time",
                "180",
                "--retry",
                "3",
                "--retry-delay",
                "1",
                "--retry-all-errors",
                "-b",
                str(cookie_path),
                "-c",
                str(cookie_path),
                "-o",
                str(tmp),
                "-w",
                "\\n%{http_code}\\n",
            ]
            if netrc_path.exists():
                cmd.extend(["--netrc-file", str(netrc_path)])
            cmd.append(url)

            result = subprocess.run(cmd, check=False, capture_output=True, text=True)
            stdout_lines = [ln.strip() for ln in (result.stdout or "").splitlines() if ln.strip()]
            status: int | None = None
            if stdout_lines and stdout_lines[-1].isdigit():
                status = int(stdout_lines[-1])
            last_status = status

            if result.returncode == 0 and status == 200 and tmp.exists():
                tmp.replace(dest)
                return True, 200, None

            stderr = (result.stderr or "").strip()
            if status == 404:
                return False, 404, stderr or None
            last_error = (
                stderr
                or (stdout_lines[-2] if len(stdout_lines) >= 2 else None)
                or f"curl-exit-{result.returncode}"
            )
            _safe_unlink(tmp)
        except Exception as exc:
            last_error = str(exc)
            _safe_unlink(tmp)
        finally:
            if cookie_path is not None:
                _safe_unlink(cookie_path)

    return False, last_status, last_error


def _gunzip(src_gz: Path, out_path: Path) -> bool:
    try:
        with gzip.open(src_gz, "rb") as src, open(out_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        return True
    except Exception as exc:
        _safe_unlink(out_path)
        log.warning("Failed to gunzip %s: %s", src_gz, exc)
        return False


# ---------------------------------------------------------------------------
# Generic daily-file sync orchestrator
# ---------------------------------------------------------------------------


def _download_candidate(
    candidate: _RemoteCandidate,
    target_dir: Path,
    *,
    netrc_path: Path | None = None,
) -> tuple[Path | None, str | None]:
    """Attempt to download a single remote candidate.  Returns (path, error)."""
    out_path = target_dir / candidate.local_name
    download_name = candidate.download_name or candidate.local_name
    download_path = target_dir / download_name

    ok, status, err = _http_download(candidate.url, download_path, netrc_path=netrc_path)
    if not ok:
        detail = (
            f"{candidate.url} (status: {status})"
            if status is not None
            else f"{candidate.url} ({err})"
        )
        return None, detail

    if not candidate.gunzip:
        return out_path, None

    if _gunzip(download_path, out_path):
        _safe_unlink(download_path)
        return out_path, None

    _safe_unlink(download_path)
    return None, f"{candidate.url} (gunzip failed)"


def _sync_daily_files(
    year_doys: Sequence[tuple[int, str]],
    *,
    product_label: str,
    target_dir: Path,
    list_local: Callable[[Path, int, str], list[Path]],
    select_preferred: Callable[[Sequence[Path]], Path],
    should_download: Callable[[date, Sequence[Path], date], bool],
    build_candidates: Callable[[int, str], list[_RemoteCandidate]],
    cleanup: Callable[[Path, int, str], None] | None = None,
    netrc_path: Path | None = None,
) -> DownloadResult:
    """Synchronise daily files for a single product type."""
    _ensure_dir(target_dir)
    reference_day = date.today()
    result = DownloadResult()

    for year, doy_str in year_doys:
        day = date(year, 1, 1) + timedelta(days=int(doy_str) - 1)
        local_matches = sorted(list_local(target_dir, year, doy_str))

        if not should_download(day, local_matches, reference_day):
            if local_matches:
                result.existing.append(select_preferred(local_matches))
            continue

        current_best = select_preferred(local_matches) if local_matches else None
        failures: list[str] = []

        for candidate in build_candidates(year, doy_str):
            if current_best is not None and candidate.local_name == current_best.name:
                result.existing.append(current_best)
                break
            candidate_path = target_dir / candidate.local_name
            if candidate_path.exists():
                result.existing.append(candidate_path)
                break

            downloaded_path, detail = _download_candidate(
                candidate, target_dir, netrc_path=netrc_path
            )
            if downloaded_path is not None:
                if cleanup is not None:
                    cleanup(target_dir, year, doy_str)
                result.downloaded.append(downloaded_path)
                log.info(
                    "\u2b07 %s %s%s downloaded: %s",
                    product_label,
                    year,
                    doy_str,
                    downloaded_path.name,
                )
                break
            if detail:
                failures.append(detail)
        else:
            if failures:
                msg = f"{product_label} {year}{doy_str}: all candidates failed:\n" + "\n".join(
                    f"  - {f}" for f in failures
                )
                log.warning(msg)
                result.failed.append(msg)

    return result
