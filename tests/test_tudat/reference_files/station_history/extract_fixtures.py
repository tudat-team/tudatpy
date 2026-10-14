#!/usr/bin/env python3
"""Extract numerical station-product fixtures from the hash-pinned source cache."""

import argparse
import hashlib
import json
from pathlib import Path


def extract(cache, output, source, destination, codes, include_parameters=True):
    metadata = json.loads((cache / (source + ".source.json")).read_text())
    original = (cache / source).read_bytes()
    if hashlib.sha256(original).hexdigest() != metadata["sha256"]:
        raise ValueError(f"Source checksum mismatch: {source}")
    lines = original.decode().splitlines()
    blocks = []
    block = None
    for line_number, line in enumerate(lines[1:], 2):
        if line.startswith("+"):
            block = [line, []]
        elif line.startswith("-"):
            if block and block[1]:
                blocks.append([*block, line])
            block = None
        elif block and not line.startswith("*"):
            if block[0] == "+SOLUTION/ESTIMATE":
                selected = include_parameters and line[14:18] in codes
            elif block[0] in (
                "+SITE/ID",
                "+SOLUTION/EPOCHS",
                "+SOLUTION/DISCONTINUITY",
                "+SITE/RECEIVER",
                "+SITE/ANTENNA",
                "+SITE/ECCENTRICITY",
            ):
                selected = line[1:5] in codes
            else:
                selected = False
            if selected:
                block[1].append((line_number, line))
    count = sum(len(records) for name, records, end in blocks if name == "+SOLUTION/ESTIMATE")
    result = [lines[0][:60] + f"{count:5d}" + lines[0][65:]]
    mapping = []
    index = 0
    for name, records, end in blocks:
        result.append(name)
        for old_line, line in records:
            original_index = int(line[1:6]) if name == "+SOLUTION/ESTIMATE" else None
            if original_index is not None:
                index += 1
                line = line[:1] + f"{index:5d}" + line[6:]
            result.append(line)
            mapping.append(
                {
                    "fixture_line": len(result),
                    "source_line": old_line,
                    "source_parameter_index": original_index,
                }
            )
        result.append(end)
    result.append("%ENDSNX")
    data = ("\n".join(result) + "\n").encode()
    (output / destination).write_bytes(data)
    return {
        "fixture": destination,
        "fixture_sha256": hashlib.sha256(data).hexdigest(),
        "source": source,
        "source_metadata": metadata,
        "line_mapping": mapping,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_cache", type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    selections = [
        ("ITRF2020-IGS-TRF.SSC", "ab07_solutions.snx", {"AB07"}, True),
        ("ITRF2020-soln-gnss.snx", "ab07_discontinuities.snx", {"AB07"}, True),
        ("ITRF2020-psd-gnss.snx", "ab07_psd.snx", {"AB07"}, True),
        ("ITRF2020-Frequencies-ENU-CM.snx", "albh_seasonal.snx", {"ALBH"}, True),
        ("ITRF2020-psd-gnss.snx", "repeated_psd.snx", {"AB48", "ANTC"}, True),
        ("ITRF2020-Frequencies-ENU-CM.snx", "duplicate_sites.snx", {"GALA", "7894"}, False),
        ("igs.snx", "abmf_equipment.snx", {"abmf"}, False),
        ("igs.snx", "cedu_equipment.snx", {"cedu"}, False),
        ("igs.snx", "igs_malformed_site.snx", {"qs07"}, False),
    ]
    sources = [extract(args.reference_cache, args.output, *selection) for selection in selections]
    (args.output / "sources.json").write_text(json.dumps(sources, indent=2) + "\n")


if __name__ == "__main__":
    main()
