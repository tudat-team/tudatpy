"""Reduce the pinned GFZ RL07 daily product without changing coefficient text."""

import gzip
import hashlib
import json
from pathlib import Path
import re
import sys


def extract(source, destination):
    """Keep degree zero through two at the first three epochs, recording source lines."""
    # Reject another provider revision before deriving the offline fixture.
    expected = "bdb32469e4f9e9408886ebe2d4ccef9c9f56e90d27a8cd56a9ad5c6ba35d1b37"

    assert hashlib.sha256(source.read_bytes()).hexdigest() == expected
    with gzip.open(source, "rt", encoding="ascii") as stream:
        lines = stream.read().splitlines()

    # Preserve the original header separately; the extract's dimensional fields are rewritten below.
    header_end = lines.index("END OF HEADER") + 1

    original_header = lines[:header_end]

    destination.mkdir(parents=True, exist_ok=True)
    (destination / "original_provider_header.txt").write_text(
        "\n".join(original_header) + "\n", encoding="ascii"
    )

    # Retain all four components so their published coefficient-level sum can be checked.
    selected = []

    provenance = []

    index = header_end

    for block in range(12):
        match = re.fullmatch(
            r"DATA SET (\d+):\s+(\d+) COEFFICIENTS FOR (.+) OF TYPE (\w+)", lines[index]
        )

        assert match is not None
        selected.append(
            f"DATA SET {block + 1:02}:  6 COEFFICIENTS FOR {match[3]} OF TYPE {match[4]}"
        )
        provenance.append(
            {
                "source_first_line": index + 2,
                "extract_first_line": header_end + len(selected) + 1,
                "count": 6,
            }
        )
        for offset in range(1, 7):
            selected.append(lines[index + offset])
        index += int(match[2]) + 1

    # Keep the provider constants/revision and change only the extract's shape, range and filename.
    replacements = {
        "MAXIMUM DEGREE": "2",
        "NUMBER OF DATA RECORDS": "72",
        "NUMBER OF DATA SETS": "12",
        "TIME LAST OBS(SEC PAST EPOCH)": "820519200.000000 (2026-01-01 06:00:00)",
        "FILENAME": "AOD1B_2026-01-01_degree2_extract_07.asc",
    }

    header = [
        (
            line.split(":", 1)[0] + ": " + replacements[line.split(":", 1)[0].strip()]
            if ":" in line and line.split(":", 1)[0].strip() in replacements
            else line
        )
        for line in original_header
    ]

    # FILESIZE is an extract property, so compute it from the final ASCII content.
    size_index = next(i for i, line in enumerate(header) if line.startswith("FILESIZE (BYTES)"))

    for _ in range(3):
        content = "\n".join(header + selected) + "\n"
        header[size_index] = "FILESIZE (BYTES)              : " + str(len(content.encode("ascii")))
    content = "\n".join(header + selected) + "\n"
    assert int(header[size_index].split(":")[1]) == len(content.encode("ascii"))
    (destination / replacements["FILENAME"]).write_text(content, encoding="ascii")

    # Record every retained coefficient's original line and all intentional header changes.
    metadata = {
        "source_url": "https://isdc-data.gfz.de/grace-fo/Level-1B/GFZ/AOD/RL07/2026/AOD1B_2026-01-01_X_07.asc.gz",
        "source_sha256": expected,
        "source_bytes": source.stat().st_size,
        "dataset_doi": "10.5880/GFZ.1.3.2022.003",
        "license": "CC-BY-4.0",
        "license_metadata_url": "https://api.datacite.org/dois/10.5880/GFZ.1.3.2022.003",
        "changed_header_fields": sorted([*replacements, "FILESIZE (BYTES)"]),
        "coefficient_line_ranges": provenance,
        "fixture_sha256": hashlib.sha256(content.encode("ascii")).hexdigest(),
    }

    (destination / "sources.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    # The source remains an external reference; no download or dynamics import occurs here.
    extract(Path(sys.argv[1]), Path(__file__).resolve().parent)
