"""Freeze the second consecutive LAGEOS-1 week before any residual evaluation."""

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

root = Path("/home/dominic/Tudat/laser-ranging-work/data/slr-week-20240407-13")
selected = []
for day in range(7, 14):
    name = f"lageos1_202404{day:02d}.np2"
    selected.append((name, f"https://edc.dgfi.tum.de/pub/slr/data/npt_crd_v2/lageos1/2024/{name}"))
for suffix in ("sp3", "sum"):
    name = f"ilrsa.orb.lageos1.240413.v81.{suffix}.gz"
    selected.append(
        (name, f"https://edc.dgfi.tum.de/pub/slr/products/orbits/TRF/lageos1/2024/240413/{name}")
    )
name = "ilrsa.pos+eop.240413.v80.snx.gz"
selected.append(
    (name, f"https://edc.dgfi.tum.de/pub/slr/products/pos+eop/weekly/2024/240413/{name}")
)


def retrieve(item):
    """Retain unchanged provider bytes and record both downloaded/decompressed identities."""
    name, url = item
    target = root / "2024-04" / name
    if target.exists():
        raise ValueError(f"Refusing to replace an existing frozen input: {name}")
    with urlopen(url, timeout=30) as response:
        raw = response.read()
    if name.endswith(".gz"):
        if not raw.startswith(b"\x1f\x8b"):
            raise ValueError(f"Not a provider gzip product: {name}")
        uncompressed = gzip.decompress(raw)
    else:
        if not raw[:20].lower().startswith(b"h1 crd"):
            raise ValueError(f"Not a provider daily CRD product: {name}")
        uncompressed = raw
    target.write_bytes(raw)
    return {
        "url": url,
        "path": f"2024-04/{name}",
        "bytes": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "uncompressed_bytes": len(uncompressed),
        "uncompressed_sha256": hashlib.sha256(uncompressed).hexdigest(),
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "rights": "Public ILRS/EDC scientific archive; retain provider citations. No explicit per-file redistribution licence observed.",
    }


records = list(ThreadPoolExecutor(max_workers=4).map(retrieve, selected))
manifest = {
    "selection_frozen_before_residuals": True,
    "week_start_utc": "2024-04-07",
    "week_end_inclusive_utc": "2024-04-13",
    "target": "LAGEOS1",
    "stations": {"7839": "11001S002", "7840": "13212S001"},
    "orbit": "ilrsa.orb.lageos1.240413.v81.sp3.gz",
    "weekly_eccentricity_metadata": "ilrsa.pos+eop.240413.v80.snx.gz",
    "selection_reason": "Next consecutive week; same combined-release conventions as first frozen week. No residual-based selection.",
    "products": records,
    "reused_metadata": [
        "metadata/SLRF2020_POS+VEL_2023.10.02.snx",
        "metadata/ILRS_Data_Handling_File_2024.02.13.snx",
        "metadata/com6.220915.tar.gz",
    ],
    "physical_selection_unchanged": "20 degree provider elevation cutoff; common8/10point interior; no residual clipping or fitted bias.",
    "residuals_computed": False,
}
(root / "download-manifest-week-20240407-13.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(
    json.dumps(
        {
            "file_count": len(records),
            "files": [{k: r[k] for k in ("path", "bytes", "sha256")} for r in records],
        },
        indent=2,
    )
)
