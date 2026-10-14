#!/usr/bin/env python3
"""Copy original provider rows after checking the pinned SHA256. No numerical model."""

import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_grid", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    metadata = json.loads((root / "sources.json").read_text())["s1_s2_def_ce.dat.source.json"]
    data = args.source_grid.read_bytes()
    if hashlib.sha256(data).hexdigest() != metadata["sha256"]:
        raise RuntimeError("Provider source checksum mismatch")
    lines = data.decode().splitlines()
    for name, longitudes, latitudes in [
        ("rp03_ce_grid_extract.txt", (47, 48), (-19, -20)),
        ("rp03_ce_polar_node_extract.txt", (0, 1), (90, 89)),
    ]:
        rows = [
            line
            for line in lines
            if not line.startswith("$")
            and line.strip()
            and float(line.split()[0]) in longitudes
            and float(line.split()[1]) in latitudes
        ]
        if len(rows) != 4:
            raise RuntimeError("Missing extract nodes")
        extracted = ("\n".join(lines[:3] + rows) + "\n").encode()
        target = root / name
        if args.check:
            if target.read_bytes() != extracted:
                raise RuntimeError("Frozen provider extract differs")
        else:
            target.write_bytes(extracted)
    print("Verified pinned RP03 grid extract" if args.check else "Extracted pinned RP03 grid rows")


if __name__ == "__main__":
    main()
