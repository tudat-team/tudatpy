from __future__ import annotations

import shutil
from datetime import datetime
from pathlib import Path

from ._common import DownloadResult
from .ionex import IonexProduct, IonexResolution, download_ionex
from .vmf import VmfTechnique, download_vmf


class DownloadAtmosphericData:
    """Downloader for atmospheric ancillary data products."""

    def download_ionex_vmf3_files(
        self,
        start_utc: str,
        end_utc: str,
        dac: str = "JPL",
        typ: str = "FIN",
        smp: str = "02H",
        vmf_technique: str = "GNSS",
        ionex_repo: str = "Data/ionex_temp",
        vmf3_repo: str = "Data/vmf3_temp",
        ionex: bool = True,
        vmf3: bool = True,
        clear_repository: bool = False,
    ) -> dict[str, DownloadResult]:
        """Download IONEX and VMF3 files over a UTC date range.

        Parameters
        ----------
        start_utc : str
            Start date/time in ISO format.
        end_utc : str
            End date/time in ISO format.
        dac : str, default="JPL"
            IONEX data analysis center.
        typ : str, default="FIN"
            IONEX product type.
        smp : str, default="02H"
            IONEX sample interval.
        vmf_technique : str, default="GNSS"
            VMF technique name.
        ionex_repo : str, default="Data/ionex_temp"
            Local IONEX output directory.
        vmf3_repo : str, default="Data/vmf3_temp"
            Local VMF3 output directory.
        ionex : bool, default=True
            Whether IONEX files are downloaded.
        vmf3 : bool, default=True
            Whether VMF3 files are downloaded.
        clear_repository : bool, default=False
            Whether enabled output directories are removed before download.

        Returns
        -------
        dict[str, DownloadResult]
            Download results keyed by product type.
        """
        if clear_repository:
            for enabled, directory in ((ionex, ionex_repo), (vmf3, vmf3_repo)):
                if enabled:
                    shutil.rmtree(directory, ignore_errors=True)

        start = datetime.fromisoformat(start_utc)
        end = datetime.fromisoformat(end_utc)

        results: dict[str, DownloadResult] = {}
        if ionex:
            dac_name = {"CODE": "COD"}.get(dac.upper(), dac.upper())
            product_name = f"{dac_name}_{'FINAL' if typ.upper() == 'FIN' else 'RAPID'}"
            resolution_name = "ONE_HOUR" if smp.upper() == "01H" else "TWO_HOUR"
            results["ionex"] = download_ionex(
                start,
                end,
                directory=Path(ionex_repo),
                products=[IonexProduct[product_name]],
                resolution=IonexResolution[resolution_name],
            )
        if vmf3:
            results["vmf3"] = download_vmf(
                start,
                end,
                technique=VmfTechnique[vmf_technique.upper()],
                directory=Path(vmf3_repo),
                day_padding=0,
            )

        return results
