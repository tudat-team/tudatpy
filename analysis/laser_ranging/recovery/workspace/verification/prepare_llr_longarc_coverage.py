"""Retrieve fixed calendar-interval primary CRD files; inspect coverage only."""

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import datetime, hashlib, json, requests

WORK = Path("/home/dominic/Tudat/laser-ranging-work")
OUT = WORK / "research/longarc-20170304"
jobs = [
    (
        f"https://edc.dgfi.tum.de/pub/slr/data/npt_crd/{target}/2017/{target}_{month}.npt",
        OUT / "crd" / f"{target}_{month}.npt",
    )
    for target in ["apollo11", "apollo14", "apollo15", "luna17", "luna21"]
    for month in ["201703", "201704"]
]


def fetch(job):
    url, path = job
    r = requests.get(url, timeout=45)
    if r.status_code == 200:
        path.write_bytes(r.content)
    return {
        "url": url,
        "file": str(path),
        "status": r.status_code,
        "bytes": len(r.content),
        "sha256": hashlib.sha256(r.content).hexdigest(),
        "retrieved_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }


with ThreadPoolExecutor(max_workers=4) as pool:
    records = list(pool.map(fetch, jobs))
(OUT / "SOURCES.json").write_text(json.dumps(records, indent=2) + "\n")
for record in records:
    print(record["status"], record["bytes"], Path(record["file"]).name)
