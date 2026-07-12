"""
Media-correction ancillary data downloaders.

Download IONEX ionospheric maps and VMF tropospheric mapping functions for any
UTC date range. IONEX-specific functionality lives in :mod:`.ionex`,
VMF-specific functionality in :mod:`.vmf`, and the shared download
infrastructure (HTTP layer, orchestrator, return type, exceptions) directly
under this package.

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
>>> from tudatpy.data.data_retrieval.media_corrections import (
...     download_ionex, download_vmf, VmfTechnique)
>>>
>>> ionex = download_ionex(datetime(2025, 3, 1), datetime(2025, 3, 7))
>>> vmf   = download_vmf(datetime(2025, 3, 1), datetime(2025, 3, 7),
...                       technique=VmfTechnique.VLBI)
>>> print(ionex.all_files)
>>> print(vmf.all_files)
"""

from __future__ import annotations

from tudatpy.kernel.data.data_retrieval.media_corrections import *

from pathlib import Path
from datetime import datetime
from typing import Sequence

from . import ionex
from . import vmf
from ._common import (
    DownloadResult,
    AncillaryDownloadError,
    AuthenticationError,
    _resolve_path,
)
from .ionex import IonexProduct, IonexResolution, download_ionex
from .vmf import VmfTechnique, VmfProcessing, download_vmf

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
    vmf_stations_to_keep: Sequence[str] | None = None,
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
    >>> from tudatpy.data.data_retrieval.media_corrections import (
    ...     download_ancillary, VmfTechnique)
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
            start,
            end,
            directory=base / "ionex",
            products=ionex_products,
            resolution=ionex_resolution,
            netrc_path=netrc_path,
        )

    if vmf:
        results["vmf"] = download_vmf(
            start,
            end,
            stations_to_keep=vmf_stations_to_keep,
            technique=vmf_technique,
            processing=vmf_processing,
            directory=base / "vmf",
        )

    return results
