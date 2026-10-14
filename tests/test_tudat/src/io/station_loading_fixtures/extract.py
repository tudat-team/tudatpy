"""Extract three unchanged IMLS station/sample rows from the pinned provider file."""

import hashlib
import json
from pathlib import Path
import sys

source = Path(sys.argv[1]).read_bytes()
digest = hashlib.sha256(source).hexdigest()
assert digest == "0a9a2e507fbc6eb66901c2653fbb16ef06fe1acaaafa6c73466168c1c82e8ca8"
lines = source.splitlines(keepends=True)
identifiers = {b"00NA", b"weuc", b"yela"}
retained = []
records = []
count_change = None
for number, line in enumerate(lines, 1):
    keep = False
    replacement = line
    if line.startswith(b"EPHEDISP"):
        keep = True
    elif line.startswith(b"#") and number <= 170:
        keep = True
    elif line.startswith(b"P "):
        replacement = b"P T 3 S          3 E      1 D          3\n"
        count_change = {
            "source_line": number,
            "original": line.decode().rstrip("\n"),
            "extract": replacement.decode().rstrip("\n"),
        }
        keep = True
    elif line.startswith((b"A ", b"T ")):
        keep = True
    elif line.startswith(b"S "):
        keep = line[3:11].rstrip() in identifiers
    elif line.startswith(b"D "):
        keep = line[45:53].rstrip() in identifiers
    if keep:
        records.append(replacement)
        if replacement == line:
            retained.append(
                {
                    "source_line": number,
                    "extract_line": len(records),
                    "sha256": hashlib.sha256(line).hexdigest(),
                }
            )
output = b"".join(records)
root = Path(__file__).parent
(root / "imls_merra2_extract.eph").write_bytes(output)
manifest = {
    "url": "https://massloading.net/imsl/load_list/atm/merra2/atm_merra2_20200101_0000.eph",
    "source_sha256": digest,
    "source_bytes": len(source),
    "extract_sha256": hashlib.sha256(output).hexdigest(),
    "extract_bytes": len(output),
    "identifiers": sorted(x.decode() for x in identifiers),
    "changed_counts": count_change,
    "unchanged_records": retained,
}
duplicates = {}
for number, line in enumerate(lines, 1):
    if line.startswith(b"S "):
        duplicates.setdefault(line[3:11].rstrip().decode(), []).append(number)
manifest["source_duplicate_site_labels"] = {
    label: numbers for label, numbers in duplicates.items() if len(numbers) > 1
}
negative = []
negative_rows = []
for number, line in enumerate(lines, 1):
    keep = line.startswith(b"EPHEDISP") or line.startswith((b"A ", b"T "))
    if line.startswith(b"S "):
        keep = line[3:11].rstrip() == b"AGGO"
    if line.startswith(b"D "):
        keep = line[45:53].rstrip() == b"AGGO"
    if line.startswith(b"P "):
        negative.append(b"P T 3 S          2 E      1 D          2\n")
    elif keep:
        negative.append(line)
        negative_rows.append(
            {
                "source_line": number,
                "extract_line": len(negative),
                "sha256": hashlib.sha256(line).hexdigest(),
            }
        )
negative_bytes = b"".join(negative)
(root / "imls_duplicate_identifier.eph").write_bytes(negative_bytes)
manifest["negative_duplicate_fixture"] = {
    "sha256": hashlib.sha256(negative_bytes).hexdigest(),
    "bytes": len(negative_bytes),
    "unchanged_records": negative_rows,
    "changed_counts": "P T 3 S          2 E      1 D          2",
}
(root / "extract_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
