"""Pin original IMLS atmosphere/hydrology files for the frozen SLR window."""

from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import requests

root = Path("/home/dominic/Tudat/laser-ranging-work/research/slr-loading-202404")
start = datetime(2024, 3, 30)
stop = datetime(2024, 4, 15)
requests_to_make = []
for component, step in [("atm", 6), ("lws", 3)]:
    epoch = start
    while epoch <= stop:
        filename = f"{component}_merra2_{epoch:%Y%m%d_%H%M}.eph"
        requests_to_make.append(
            (
                component,
                filename,
                f"https://massloading.net/imsl/load_list/{component}/merra2/{filename}",
            )
        )
        epoch += timedelta(hours=step)


def retrieve(item):
    """Download an original single-epoch product, recording its bytes and site records."""
    component, filename, url = item
    path = root / filename
    if not path.exists():
        response = requests.get(url, timeout=35)
        response.raise_for_status()
        if not response.content.startswith(b"EPHEDISP"):
            raise ValueError(f"{filename} did not return EPHEDISP data")
        path.write_bytes(response.content)
    raw = path.read_bytes()
    lines = raw.decode("ascii").splitlines()
    site_lines = [
        line for line in lines if line.startswith("S ") and line.split()[1] in ("Grzl", "Herl")
    ]
    displacements = [
        line for line in lines if line.startswith("D ") and line.split()[5] in ("Grzl", "Herl")
    ]
    assert len(site_lines) == len(displacements) == 2, (filename, site_lines, displacements)
    begin = next(line for line in lines if line.startswith("T begin"))
    return {
        "component": component,
        "file": filename,
        "url": url,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "bytes": len(raw),
        "epoch_header": begin,
        "site_records": site_lines,
        "displacement_records": displacements,
    }


records = []
failures = []
with ThreadPoolExecutor(max_workers=4) as pool:
    futures = {pool.submit(retrieve, item): item for item in requests_to_make}
    for future in as_completed(futures):
        try:
            records.append(future.result())
        except Exception as error:
            failures.append({"file": futures[future][1], "error": str(error)})
report = {
    "retrieved_utc": datetime.now(timezone.utc).isoformat(),
    "window": "2024-03-30T00:00:00 through 2024-04-15T00:00:00; provider TAI epochs authoritative",
    "requested_files": len(requests_to_make),
    "records": sorted(records, key=lambda r: r["file"]),
    "failures": failures,
}
(root / "history_sources.json").write_text(json.dumps(report, indent=2) + "\n")
print(
    "retrieved", len(records), "failures", len(failures), "bytes", sum(r["bytes"] for r in records)
)
for item in failures:
    print(item)
