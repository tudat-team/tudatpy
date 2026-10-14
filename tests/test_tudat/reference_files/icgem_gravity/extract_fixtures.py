"""Reproduce unchanged static fixture extracts from pinned local source files."""

from pathlib import Path
import hashlib
import json
import sys

destination = Path(__file__).parent
cache = Path(sys.argv[1])
records = json.loads((destination / "sources.json").read_text())["fixtures"]
sources = [cache / "gfz-rl061-ddk4-202104.gfc", cache / "costg/models/EIGEN6-C4.gfc"]
for record, source in zip(records, sources):
    contents = source.read_bytes()
    assert hashlib.sha256(contents).hexdigest() == record["source_sha256"]
    lines = contents.splitlines(keepends=True)
    if source.name.startswith("gfz-"):
        end = next(
            i for i, line in enumerate(lines) if line.split() and line.split()[0] == b"end_of_head"
        )
        selected = lines[: end + 1] + [
            line
            for line in lines[end + 1 :]
            if line.split() and line.split()[0] == b"gfc" and int(line.split()[1]) <= 2
        ]
    else:
        selected = lines[:16]
    extract = b"".join(selected)
    assert hashlib.sha256(extract).hexdigest() == record["fixture_sha256"]
    (destination / record["file"]).write_bytes(extract)
