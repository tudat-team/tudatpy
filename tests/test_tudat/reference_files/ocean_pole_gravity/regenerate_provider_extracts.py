"""Offline reproduction: python regenerate_provider_extracts.py CACHE [FIXTURE_DIRECTORY]."""

import gzip
import hashlib
import json
from pathlib import Path
import sys

fixtures = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).resolve().parent
cache = Path(sys.argv[1])
for name, record in json.loads((fixtures / "manifest.json").read_text()).items():
    source = cache / record["source"]
    payload = source.read_bytes()
    assert hashlib.sha256(payload).hexdigest() == record["source_sha256"], source
    if source.suffix == ".gz":
        payload = gzip.decompress(payload)
    lines = payload.splitlines(keepends=True)
    ranges = record["source_lines"]
    if isinstance(ranges[0], int):
        ranges = [ranges]
    extract = b"".join(b"".join(lines[first - 1 : last]) for first, last in ranges)
    assert record["edits"] == []
    assert hashlib.sha256(extract).hexdigest() == record["extract_sha256"], name
    (fixtures / name).write_bytes(extract)
