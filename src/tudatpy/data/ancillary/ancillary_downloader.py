"""
Ancillary data downloader for space geodesy applications.

Download IONEX ionospheric maps and VMF tropospheric mapping functions
for any UTC date range.

Supported data products
-----------------------
- **IONEX** : Global ionospheric maps from NASA CDDIS (requires .netrc credentials).
  Products from JPL, IGS, CODE, and UPC analysis centres in final and rapid solutions.
- **VMF**  : Vienna Mapping Functions (VMF3) for GNSS, VLBI, DORIS, and SLR (optical)
  from the Vienna University of Technology.

References
----------
- S. Schaer, W. Gurtner, and J. Feltens (1998), IONEX: The IONosphere Map EXchange Format
  Version 1, Proceedings of the IGS AC Workshop, Darmstadt, Germany.
- D. Landskron and J. Bohm (2018), VMF3/GPT3: refined discrete and empirical troposphere
  mapping functions, J. Geodesy, 92(4), doi:10.1007/s00190-017-1066-2

Example
-------
>>> from datetime import datetime
>>> from tudatpy.data.ancillary import download_ionex, download_vmf, VmfTechnique
>>>
>>> ionex = download_ionex(datetime(2025, 3, 1), datetime(2025, 3, 7))
>>> vmf   = download_vmf(datetime(2025, 3, 1), datetime(2025, 3, 7),
...                       technique=VmfTechnique.VLBI)
>>> print(ionex.all_files)
>>> print(vmf.all_files)
"""

from __future__ import annotations

__all__ = [
    "IonexProduct",
    "IonexResolution",
    "VmfTechnique",
    "VmfProcessing",
    "DownloadResult",
    "AncillaryDownloadError",
    "AuthenticationError",
    "download_ionex",
    "download_vmf",
    "download_ancillary",
]

import gzip
import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from datetime import date, datetime, timedelta
from enum import Enum
from pathlib import Path
from typing import Callable, Sequence

try:
    import requests as _requests
except ImportError:
    _requests = None  # type: ignore[assignment]

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Public enums
# ---------------------------------------------------------------------------


class IonexProduct(Enum):
    """IONEX analysis centre and solution type, ordered by preference.

    Each member corresponds to a specific analysis centre (JPL, IGS, CODE, UPC)
    and solution type (final or rapid). The enum ordering reflects the default
    priority: final solutions from more precise centres are preferred.

    Members
    -------
    JPL_FINAL
        JPL final GIM product (``JPL0OPSFIN``). Bi-cubic spline model, Kalman
        filter, typically highest accuracy. 2-hour resolution.
    JPL_RAPID
        JPL rapid GIM product (``JPL0OPSRAP``). Same method as final, lower
        latency. 2-hour resolution.
    IGS_FINAL
        IGS combined final product (``IGS0OPSFIN``). Weighted combination of
        multiple analysis centres. 2-hour resolution.
    COD_FINAL
        CODE (Bern) final product (``COD0OPSFIN``). Spherical harmonic
        expansion (degree/order 15). 1-hour resolution.
    UPC_FINAL
        UPC-IonSAT final product (``UPC0OPSFIN``). Tomographic model. 2-hour
        resolution.
    IGS_RAPID
        IGS combined rapid product (``IGS0OPSRAP``). 2-hour resolution.
    COD_RAPID
        CODE rapid product (``COD0OPSRAP``). 1-hour resolution.
    """

    JPL_FINAL = "JPL0OPSFIN"
    JPL_RAPID = "JPL0OPSRAP"
    IGS_FINAL = "IGS0OPSFIN"
    COD_FINAL = "COD0OPSFIN"
    UPC_FINAL = "UPC0OPSFIN"
    IGS_RAPID = "IGS0OPSRAP"
    COD_RAPID = "COD0OPSRAP"


class IonexResolution(Enum):
    """Temporal resolution of IONEX maps.

    Members
    -------
    ONE_HOUR
        1-hour interval between TEC maps (25 maps per day). Available for
        CODE products.
    TWO_HOUR
        2-hour interval between TEC maps (13 maps per day). Available for
        JPL, IGS, UPC products.
    """

    ONE_HOUR = "01H"
    TWO_HOUR = "02H"


# Default resolution per product (CDDIS naming convention).
_IONEX_DEFAULT_RESOLUTION: dict[str, str] = {
    "JPL0OPSFIN": "02H",
    "JPL0OPSRAP": "02H",
    "IGS0OPSFIN": "02H",
    "COD0OPSFIN": "01H",
    "UPC0OPSFIN": "02H",
    "IGS0OPSRAP": "02H",
    "COD0OPSRAP": "01H",
}


class VmfTechnique(Enum):
    """Geodetic technique for Vienna Mapping Function products.

    Different techniques use different site networks and produce coefficient
    files with different column structures and station grids.

    Members
    -------
    GNSS
        GNSS station network coefficients (``.v3gr_g``).
    VLBI
        VLBI station network coefficients (``.v3gr_r``).
    DORIS
        DORIS station network coefficients (``.v3gr_d``).
    OPTICAL
        Optical/SLR station network coefficients (``.vmf3o``). Uses the VMF3o
        model with wavelength-dependent scaling.
    """

    GNSS = "gnss"
    VLBI = "vlbi"
    DORIS = "doris"
    OPTICAL = "optical"


class VmfProcessing(Enum):
    """VMF processing type.

    Members
    -------
    OPERATIONAL
        Near-real-time operational processing (``OP``). Lowest latency,
        uses operational NWP analysis fields.
    ERA_INTERIM
        ERA-Interim reanalysis processing (``EI``). Higher consistency
        for long-term studies, available only for historical periods.
    FORECAST
        Forecast processing (``FC``). Based on NWP forecast fields.
    """

    OPERATIONAL = "OP"
    ERA_INTERIM = "EI"
    FORECAST = "FC"


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
        """All available files (downloaded + existing), sorted by name."""
        return sorted(set(self.downloaded + self.existing))

    @property
    def success(self) -> bool:
        """True when at least one file is available and nothing failed."""
        return len(self.all_files) > 0 and len(self.failed) == 0


# ---------------------------------------------------------------------------
# VMF product registry
# ---------------------------------------------------------------------------

_VMF_BASE_URL = "https://vmf.geo.tuwien.ac.at/trop_products"


@dataclass(frozen=True)
class _VmfProductSpec:
    technique_dir: str
    product_dir: str
    suffix: str

    def url(self, processing: VmfProcessing, year: int, doy_str: str) -> str:
        name = f"{year}{doy_str}.{self.suffix}"
        return (
            f"{_VMF_BASE_URL}/{self.technique_dir}/{self.product_dir}/"
            f"{self.product_dir}_{processing.value}/daily/{year}/{name}"
        )


_VMF_REGISTRY: dict[VmfTechnique, _VmfProductSpec] = {
    VmfTechnique.GNSS: _VmfProductSpec("GNSS", "V3GR", "v3gr_g"),
    VmfTechnique.VLBI: _VmfProductSpec("VLBI", "V3GR", "v3gr_r"),
    VmfTechnique.DORIS: _VmfProductSpec("DORIS", "V3GR", "v3gr_d"),
    VmfTechnique.OPTICAL: _VmfProductSpec("SLR", "VMF3o", "vmf3o"),
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_HTTP_CHUNK_SIZE = 1024 * 1024
_HTTP_MAX_ATTEMPTS = 3
_IONEX_RECENT_REFRESH_DAYS = 7


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
                    int(content_length)
                    if content_length and content_length.isdigit()
                    else None
                )
                written = 0
                with open(tmp, "wb") as fh:
                    for chunk in resp.iter_content(chunk_size=_HTTP_CHUNK_SIZE):
                        if chunk:
                            fh.write(chunk)
                            written += len(chunk)

                if expected is not None and written != expected:
                    raise IOError(
                        f"Incomplete download: expected {expected} bytes, got {written}."
                    )
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
                "-sS", "-L",
                "--connect-timeout", "10",
                "--max-time", "180",
                "--retry", "3",
                "--retry-delay", "1",
                "--retry-all-errors",
                "-b", str(cookie_path),
                "-c", str(cookie_path),
                "-o", str(tmp),
                "-w", "\\n%{http_code}\\n",
            ]
            if netrc_path.exists():
                cmd.extend(["--netrc-file", str(netrc_path)])
            cmd.append(url)

            result = subprocess.run(cmd, check=False, capture_output=True, text=True)
            stdout_lines = [
                ln.strip() for ln in (result.stdout or "").splitlines() if ln.strip()
            ]
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

    ok, status, err = _http_download(
        candidate.url, download_path, netrc_path=netrc_path
    )
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
                    product_label, year, doy_str, downloaded_path.name,
                )
                break
            if detail:
                failures.append(detail)
        else:
            if failures:
                msg = (
                    f"{product_label} {year}{doy_str}: all candidates failed:\n"
                    + "\n".join(f"  - {f}" for f in failures)
                )
                log.warning(msg)
                result.failed.append(msg)

    return result


# ---------------------------------------------------------------------------
# IONEX
# ---------------------------------------------------------------------------

_IONEX_BASE_URL = "https://cddis.nasa.gov/archive/gnss/products/ionex"


def _ionex_local_pattern(year: int, doy_str: str) -> str:
    return f"*_{year}{doy_str}0000_*_GIM.INX"


def _build_ionex_candidates(
    year: int,
    doy_str: str,
    products: Sequence[IonexProduct],
    resolution: IonexResolution | None,
) -> list[_RemoteCandidate]:
    candidates: list[_RemoteCandidate] = []
    for product in products:
        res = (
            resolution.value
            if resolution is not None
            else _IONEX_DEFAULT_RESOLUTION.get(product.value, "02H")
        )
        local_name = f"{product.value}_{year}{doy_str}0000_01D_{res}_GIM.INX"
        name_gz = f"{local_name}.gz"
        candidates.append(
            _RemoteCandidate(
                local_name=local_name,
                url=f"{_IONEX_BASE_URL}/{year}/{doy_str}/{name_gz}",
                download_name=name_gz,
                gunzip=True,
            )
        )
    return candidates


def _ionex_preference_key(
    path: Path,
    products: Sequence[IonexProduct],
) -> tuple[int, str]:
    name = path.name
    for idx, product in enumerate(products):
        if name.startswith(f"{product.value}_"):
            return idx, name
    return len(products), name


def _should_refresh_ionex(
    day: date,
    matches: Sequence[Path],
    reference_day: date,
    *,
    products: Sequence[IonexProduct],
    refresh_within_days: int,
) -> bool:
    if not matches:
        return True
    if (reference_day - day).days > refresh_within_days:
        return False
    preferred = min(matches, key=lambda p: _ionex_preference_key(p, products))
    return not preferred.name.startswith(f"{products[0].value}_")


def _cleanup_ionex_partials(directory: Path, year: int, doy_str: str) -> None:
    for part in directory.glob(f"*_{year}{doy_str}0000_*_GIM.INX.gz.part"):
        _safe_unlink(part)


def download_ionex(
    start: datetime,
    end: datetime,
    *,
    directory: Path | str = "~/.tudat/ancillary/ionex",
    products: Sequence[IonexProduct] | None = None,
    resolution: IonexResolution | None = IonexResolution.TWO_HOUR,
    netrc_path: Path | str | None = None,
    refresh_within_days: int = _IONEX_RECENT_REFRESH_DAYS,
) -> DownloadResult:
    """Download IONEX ionospheric map files from NASA CDDIS.

    Downloads global ionospheric TEC maps in IONEX format for each day in the
    given date range. Files are fetched from the NASA CDDIS archive, which
    requires Earthdata credentials in a ``.netrc`` file.

    The downloader uses a **product priority fallback**: for each day it tries
    the preferred product first (e.g., JPL final), and if that is not available
    it falls back to the next product in the list (e.g., IGS final, then CODE,
    etc.). For recent days within ``refresh_within_days``, the downloader will
    attempt to upgrade from rapid to final products on subsequent calls.

    Parameters
    ----------
    start : datetime
        Start of the UTC date range (inclusive).
    end : datetime
        End of the UTC date range (inclusive).
    directory : Path or str, default = "~/.tudat/ancillary/ionex"
        Local storage directory (created if needed).
    products : sequence of IonexProduct, optional
        Ordered preference list of IONEX products to try. Defaults to all
        products in priority order (JPL final first, CODE rapid last).
    resolution : IonexResolution or None, default = IonexResolution.TWO_HOUR
        Temporal resolution override. ``None`` uses each product's native
        resolution (e.g., 1 h for CODE, 2 h for JPL).
    netrc_path : Path or str, optional
        Path to ``.netrc`` file for CDDIS authentication. Defaults to
        ``~/.netrc``.
    refresh_within_days : int, default = 7
        Re-attempt download for days within this window of today to upgrade
        from rapid to final products.

    Returns
    -------
    DownloadResult
        Result object with ``downloaded``, ``existing``, ``failed``, and
        ``all_files`` attributes.

    Raises
    ------
    AuthenticationError
        If no ``.netrc`` file is found (required for CDDIS access).

    Examples
    --------
    >>> from datetime import datetime
    >>> from tudatpy.data.ancillary import download_ionex, IonexProduct
    >>>
    >>> # Download with default settings (JPL preferred, 2h resolution)
    >>> result = download_ionex(datetime(2025, 3, 1), datetime(2025, 3, 7))
    >>> print(result.all_files)
    >>>
    >>> # Download only JPL products
    >>> result = download_ionex(
    ...     datetime(2025, 3, 1), datetime(2025, 3, 7),
    ...     products=[IonexProduct.JPL_FINAL, IonexProduct.JPL_RAPID]
    ... )
    """
    if products is None:
        products = list(IonexProduct)
    target_dir = _resolve_path(directory)
    resolved_netrc = _resolve_path(netrc_path) if netrc_path is not None else None

    if resolved_netrc is None:
        default_netrc = Path.home() / ".netrc"
        if not default_netrc.exists():
            raise AuthenticationError(
                "IONEX downloads from CDDIS require a .netrc file for "
                "authentication.\n\n"
                "Create ~/.netrc with the following content:\n\n"
                "  machine urs.earthdata.nasa.gov\n"
                "      login <your_earthdata_username>\n"
                "      password <your_earthdata_password>\n\n"
                "Register at https://urs.earthdata.nasa.gov if needed.\n"
                "Then: chmod 600 ~/.netrc"
            )

    year_doys = _date_range_to_year_doys(start, end)

    return _sync_daily_files(
        year_doys,
        product_label="IONEX",
        target_dir=target_dir,
        list_local=lambda d, y, doy: sorted(
            d.glob(_ionex_local_pattern(y, doy))
        ),
        select_preferred=lambda matches: min(
            matches, key=lambda p: _ionex_preference_key(p, products)
        ),
        should_download=lambda day, matches, ref: _should_refresh_ionex(
            day, matches, ref,
            products=products,
            refresh_within_days=refresh_within_days,
        ),
        build_candidates=lambda y, doy: _build_ionex_candidates(
            y, doy, products, resolution
        ),
        cleanup=_cleanup_ionex_partials,
        netrc_path=resolved_netrc,
    )


# ---------------------------------------------------------------------------
# VMF (Vienna Mapping Functions)
# ---------------------------------------------------------------------------


def _build_vmf_candidates(
    year: int,
    doy_str: str,
    spec: _VmfProductSpec,
    processing: VmfProcessing,
) -> list[_RemoteCandidate]:
    name = f"{year}{doy_str}.{spec.suffix}"
    return [
        _RemoteCandidate(
            local_name=name,
            url=spec.url(processing, year, doy_str),
        )
    ]


def _vmf_local_matches(
    directory: Path,
    year: int,
    doy_str: str,
    suffix: str,
) -> list[Path]:
    path = directory / f"{year}{doy_str}.{suffix}"
    return [path] if path.is_file() else []


def _should_download_missing_only(
    _day: date,
    matches: Sequence[Path],
    _reference_day: date,
) -> bool:
    return len(matches) == 0


def download_vmf(
    start: datetime,
    end: datetime,
    *,
    technique: VmfTechnique = VmfTechnique.GNSS,
    processing: VmfProcessing = VmfProcessing.OPERATIONAL,
    directory: Path | str = "~/.tudat/ancillary/vmf",
    day_padding: int = 1,
) -> DownloadResult:
    """Download VMF troposphere mapping function files from Vienna.

    Downloads daily VMF3 tropospheric mapping function coefficient files from
    the Vienna University of Technology server for each day in the given date
    range (plus optional padding days for interpolation continuity).

    Parameters
    ----------
    start : datetime
        Start of the UTC date range (inclusive).
    end : datetime
        End of the UTC date range (inclusive).
    technique : VmfTechnique, default = VmfTechnique.GNSS
        Geodetic technique: ``GNSS``, ``VLBI``, ``DORIS``, or ``OPTICAL``
        (SLR). Determines the station network and coefficient format.
    processing : VmfProcessing, default = VmfProcessing.OPERATIONAL
        Processing type: ``OPERATIONAL`` (near-real-time), ``ERA_INTERIM``
        (reanalysis, for historical studies), or ``FORECAST``.
    directory : Path or str, default = "~/.tudat/ancillary/vmf"
        Local storage directory (created if needed).
    day_padding : int, default = 1
        Extra days to download before and after the range to ensure
        interpolation continuity at arc boundaries.

    Returns
    -------
    DownloadResult
        Result object with ``downloaded``, ``existing``, ``failed``, and
        ``all_files`` attributes.

    Examples
    --------
    >>> from datetime import datetime
    >>> from tudatpy.data.ancillary import download_vmf, VmfTechnique
    >>>
    >>> # Download GNSS troposphere data (default)
    >>> result = download_vmf(datetime(2025, 3, 1), datetime(2025, 3, 7))
    >>>
    >>> # Download VLBI troposphere with ERA-Interim reanalysis
    >>> result = download_vmf(
    ...     datetime(2025, 3, 1), datetime(2025, 3, 7),
    ...     technique=VmfTechnique.VLBI,
    ...     processing=VmfProcessing.ERA_INTERIM
    ... )
    """
    spec = _VMF_REGISTRY[technique]
    target_dir = _resolve_path(directory)
    year_doys = _date_range_to_year_doys(start, end, padding=day_padding)

    return _sync_daily_files(
        year_doys,
        product_label=f"VMF ({technique.name})",
        target_dir=target_dir,
        list_local=lambda d, y, doy: _vmf_local_matches(d, y, doy, spec.suffix),
        select_preferred=lambda matches: matches[0],
        should_download=_should_download_missing_only,
        build_candidates=lambda y, doy: _build_vmf_candidates(
            y, doy, spec, processing
        ),
    )


# ---------------------------------------------------------------------------
# Convenience wrapper
# ---------------------------------------------------------------------------


def download_ancillary(
    start: datetime,
    end: datetime,
    *,
    ionex: bool = True,
    vmf: bool = True,
    vmf_technique: VmfTechnique = VmfTechnique.GNSS,
    vmf_processing: VmfProcessing = VmfProcessing.OPERATIONAL,
    ionex_products: Sequence[IonexProduct] | None = None,
    ionex_resolution: IonexResolution | None = IonexResolution.TWO_HOUR,
    base_directory: Path | str = "~/.tudat/ancillary",
    netrc_path: Path | str | None = None,
) -> dict[str, DownloadResult]:
    """Download all ancillary atmospheric and ionospheric data for a date range.

    Convenience wrapper around :func:`download_ionex` and :func:`download_vmf`
    that downloads all requested data products into organized sub-directories.

    Parameters
    ----------
    start : datetime
        Start of the UTC date range (inclusive).
    end : datetime
        End of the UTC date range (inclusive).
    ionex : bool, default = True
        Download IONEX ionospheric map files.
    vmf : bool, default = True
        Download VMF troposphere mapping function files.
    vmf_technique : VmfTechnique, default = VmfTechnique.GNSS
        Geodetic technique for VMF products.
    vmf_processing : VmfProcessing, default = VmfProcessing.OPERATIONAL
        VMF processing type.
    ionex_products : sequence of IonexProduct, optional
        Ordered preference list for IONEX. Defaults to all products.
    ionex_resolution : IonexResolution or None, default = IonexResolution.TWO_HOUR
        IONEX temporal resolution override.
    base_directory : Path or str, default = "~/.tudat/ancillary"
        Root directory; sub-folders ``ionex/`` and ``vmf/`` are created
        underneath.
    netrc_path : Path or str, optional
        Path to ``.netrc`` for CDDIS authentication.

    Returns
    -------
    dict[str, DownloadResult]
        Keys ``"ionex"`` and/or ``"vmf"`` (only present if requested), each
        mapping to a :class:`DownloadResult`.

    Examples
    --------
    >>> from datetime import datetime
    >>> from tudatpy.data.ancillary import download_ancillary, VmfTechnique
    >>>
    >>> # Download everything (IONEX + VMF GNSS)
    >>> results = download_ancillary(datetime(2025, 3, 1), datetime(2025, 3, 7))
    >>> print(results["ionex"].all_files)
    >>> print(results["vmf"].all_files)
    >>>
    >>> # Download only VMF with optical/SLR technique
    >>> results = download_ancillary(
    ...     datetime(2025, 3, 1), datetime(2025, 3, 7),
    ...     ionex=False, vmf_technique=VmfTechnique.OPTICAL
    ... )
    """
    base = _resolve_path(base_directory)
    results: dict[str, DownloadResult] = {}

    if ionex:
        results["ionex"] = download_ionex(
            start, end,
            directory=base / "ionex",
            products=ionex_products,
            resolution=ionex_resolution,
            netrc_path=netrc_path,
        )

    if vmf:
        results["vmf"] = download_vmf(
            start, end,
            technique=vmf_technique,
            processing=vmf_processing,
            directory=base / "vmf",
        )

    return results
