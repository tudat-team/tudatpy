#!/usr/bin/env python3
"""Extract pinned official IERS data locally; no downloads or production expressions."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_directory", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    sources = json.loads((root / "sources.json").read_text())
    for source in sources:
        path = args.reference_directory / source["file"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != source["sha256"]:
            raise RuntimeError(f"Source checksum mismatch: {path}")
    lines = (
        gzip.decompress((args.reference_directory / sources[0]["file"]).read_bytes())
        .decode()
        .splitlines()
    )
    headers = lines[:14]
    for i in (3, 7):
        headers[i] = headers[i].split("=")[0] + "= 2"
    for i, value in [(4, 232.25), (5, 232.75), (8, -43.75), (9, -43.25)]:
        headers[i] = headers[i].split("=")[0] + "= " + str(value)
    nodes = [
        line
        for line in lines[14:]
        if float(line.split()[0]) in (232.25, 232.75) and float(line.split()[1]) in (-43.75, -43.25)
    ]
    if len(nodes) != 4:
        raise RuntimeError("Missing reference grid nodes")
    products = {
        "iers_grid_extract.txt": ("\n".join(headers + nodes) + "\n").encode(),
        "iers_published_test.txt": (args.reference_directory / sources[1]["file"]).read_bytes(),
    }
    for name, contents in products.items():
        if args.check:
            if (root / name).read_bytes() != contents:
                raise RuntimeError(f"Frozen extract differs: {name}")
        else:
            (root / name).write_bytes(contents)
    print(
        "Verified pinned ocean pole-tide fixtures"
        if args.check
        else "Extracted pinned ocean pole-tide fixtures"
    )


if __name__ == "__main__":
    main()
