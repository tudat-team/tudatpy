from __future__ import annotations

import shutil
from datetime import datetime
from pathlib import Path

from tudatpy.data_access.downloading.media_corrections import (
    DownloadResult,
    IonexProduct,
    IonexResolution,
    VmfTechnique,
    download_ionex,
    download_vmf,
)


class DownloadAtmosphericData:
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
