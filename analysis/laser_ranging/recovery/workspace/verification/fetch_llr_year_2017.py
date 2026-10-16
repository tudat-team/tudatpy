"""Freeze the predeclared calendar2017 provider files without computing residuals."""

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import datetime, hashlib, json, urllib.request, shutil

ROOT = Path("/home/dominic/Tudat/laser-ranging-work/research/llr-year-2017")
OLD = Path("/home/dominic/Tudat/laser-ranging-work/research/longarc-20170304")
old = {Path(r["file"]).name: r for r in json.loads((OLD / "SOURCES.json").read_text())}


def get(target, month):
    name = f"{target}_2017{month:02d}.npt"
    url = f"https://edc.dgfi.tum.de/pub/slr/data/npt_crd/{target}/2017/{name}"
    path = ROOT / "crd" / name
    result = dict(
        url=url,
        file="crd/" + name,
        target=target,
        month=month,
        retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    )
    try:
        if name in old:
            payload = (OLD / "crd" / name).read_bytes()
            assert hashlib.sha256(payload).hexdigest() == old[name]["sha256"]
            result.update(status=200, reused_from_frozen_two_month_source=True)
        else:
            with urllib.request.urlopen(url, timeout=30) as response:
                payload = response.read()
                result["status"] = response.status
            if not payload.lstrip().lower().startswith(b"h1 crd"):
                raise ValueError("Non-CRD response")
        path.write_bytes(payload)
        result.update(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
    except Exception as error:
        result["error"] = str(error)
    return result


items = [
    (target, month)
    for target in ["apollo11", "apollo14", "apollo15", "luna17", "luna21"]
    for month in range(1, 13)
]
with ThreadPoolExecutor(max_workers=6) as pool:
    result = list(pool.map(lambda item: get(*item), items))
(ROOT / "SOURCES.json").write_text(json.dumps(result, indent=2) + "\n")
for target in ["apollo11", "apollo14", "apollo15", "luna17", "luna21"]:
    url = f"https://edc.dgfi.tum.de/pub/slr/data/npt_crd/{target}/2017/"
    try:
        data = urllib.request.urlopen(url, timeout=30).read()
        (ROOT / (target + "-2017-provider-index.html")).write_bytes(data)
    except Exception as error:
        print("INDEX_ERROR", target, str(error), flush=True)
print(
    json.dumps(
        {
            "retrieved": sum("sha256" in r for r in result),
            "reused": sum(r.get("reused_from_frozen_two_month_source", False) for r in result),
            "failures": [r for r in result if "error" in r],
        },
        indent=2,
    )
)
