"""
VMF (Vienna Mapping Functions) tropospheric mapping downloader.

Downloads daily VMF3 tropospheric mapping function coefficient files from the
Vienna University of Technology server, for GNSS, VLBI, DORIS, and SLR (optical)
techniques.

References
----------
- D. Landskron and J. Bohm (2018), VMF3/GPT3: refined discrete and empirical
  troposphere mapping functions, J. Geodesy, 92(4),
  doi:10.1007/s00190-017-1066-2
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import date, datetime
from enum import Enum
from pathlib import Path
from typing import Sequence

from .._common import (
    DownloadResult,
    _RemoteCandidate,
    _date_range_to_year_doys,
    _resolve_path,
    _sync_daily_files,
)

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public enums
# ---------------------------------------------------------------------------


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
# VMF product registry
# ---------------------------------------------------------------------------

_VMF_BASE_URL = "https://vmf.geo.tuwien.ac.at/trop_products"


@dataclass(frozen=True)
class _VmfProductSpec:
    """Internal VMF product URL specification."""

    technique_dir: str
    product_dir: str
    suffix: str

    def url(self, processing: VmfProcessing, year: int, doy_str: str) -> str:
        """Build the remote VMF product URL.

        Parameters
        ----------
        processing : VmfProcessing
            VMF processing mode.
        year : int
            Calendar year.
        doy_str : str
            Day-of-year string.

        Returns
        -------
        str
            Remote VMF product URL.
        """
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


def _filter_vmf_stations(file_path: Path, stations: Sequence[str]) -> None:
    """Filter a VMF file to keep only the specified stations (case-insensitive)."""
    stations_upper = {s.upper() for s in stations}
    lines = file_path.read_text().splitlines()
    filtered = [ln for ln in lines if ln.split()[0].upper() in stations_upper]
    file_path.write_text("\n".join(filtered) + "\n")


def download_vmf(
    start: datetime,
    end: datetime,
    *,
    technique: VmfTechnique = VmfTechnique.GNSS,
    processing: VmfProcessing = VmfProcessing.OPERATIONAL,
    directory: Path | str = "~/.tudat/ancillary/vmf",
    day_padding: int = 1,
    stations_to_keep: Sequence[str] | None = None,
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
    stations_to_keep : sequence of str, optional
        If provided, filter each downloaded VMF file to keep only lines for
        these station codes (case-insensitive, e.g. ``["OPMT", "WTZR"]``).
        This significantly reduces memory usage when loading VMF data into
        Tudat via ``set_vmf_troposphere_data``, especially for global GNSS
        networks with 300+ stations.

    Returns
    -------
    DownloadResult
        Result object with ``downloaded``, ``existing``, ``failed``, and
        ``all_files`` attributes.

    Examples
    --------
    >>> from datetime import datetime
    >>> from tudatpy.data.data_retrieval.media_corrections.vmf import (
    ...     download_vmf, VmfTechnique)
    >>>
    >>> # Download GNSS troposphere data (default)
    >>> result = download_vmf(datetime(2025, 3, 1), datetime(2025, 3, 7))
    >>>
    >>> # Download only for specific stations
    >>> result = download_vmf(
    ...     datetime(2025, 3, 1), datetime(2025, 3, 7),
    ...     stations_to_keep=["OPMT", "WTZR", "HERS"]
    ... )
    """
    spec = _VMF_REGISTRY[technique]
    target_dir = _resolve_path(directory)
    year_doys = _date_range_to_year_doys(start, end, padding=day_padding)

    result = _sync_daily_files(
        year_doys,
        product_label=f"VMF ({technique.name})",
        target_dir=target_dir,
        list_local=lambda d, y, doy: _vmf_local_matches(d, y, doy, spec.suffix),
        select_preferred=lambda matches: matches[0],
        should_download=_should_download_missing_only,
        build_candidates=lambda y, doy: _build_vmf_candidates(y, doy, spec, processing),
    )

    # Filter stations after download
    if stations_to_keep is not None:
        for f in result.all_files:
            try:
                _filter_vmf_stations(f, stations_to_keep)
            except Exception as exc:
                log.warning("Failed to filter stations in %s: %s", f, exc)

    return result
