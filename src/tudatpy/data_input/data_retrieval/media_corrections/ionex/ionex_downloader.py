"""
IONEX ionospheric map downloader.

Downloads global ionospheric TEC maps in IONEX format from the NASA CDDIS
archive (requires ``.netrc`` credentials). Products from JPL, IGS, CODE, and
UPC analysis centres in final and rapid solutions are supported.

References
----------
- S. Schaer, W. Gurtner, and J. Feltens (1998), IONEX: The IONosphere Map
  EXchange Format Version 1, Proceedings of the IGS AC Workshop, Darmstadt,
  Germany.
"""

from __future__ import annotations

from datetime import date, datetime
from enum import Enum
from pathlib import Path
from typing import Sequence

from .._common import (
    AuthenticationError,
    DownloadResult,
    _RemoteCandidate,
    _date_range_to_year_doys,
    _resolve_path,
    _safe_unlink,
    _sync_daily_files,
)

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


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_IONEX_BASE_URL = "https://cddis.nasa.gov/archive/gnss/products/ionex"
_IONEX_RECENT_REFRESH_DAYS = 7


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
    >>> from tudatpy.data_input.data_retrieval.media_corrections.ionex import (
    ...     download_ionex, IonexProduct)
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
        list_local=lambda d, y, doy: sorted(d.glob(_ionex_local_pattern(y, doy))),
        select_preferred=lambda matches: min(
            matches, key=lambda p: _ionex_preference_key(p, products)
        ),
        should_download=lambda day, matches, ref: _should_refresh_ionex(
            day,
            matches,
            ref,
            products=products,
            refresh_within_days=refresh_within_days,
        ),
        build_candidates=lambda y, doy: _build_ionex_candidates(y, doy, products, resolution),
        cleanup=_cleanup_ionex_partials,
        netrc_path=resolved_netrc,
    )
