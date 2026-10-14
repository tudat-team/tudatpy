"""Extract unchanged IMLS source ranges, with only explicit P count edits."""

import hashlib
import json
from pathlib import Path
import sys

source = Path(sys.argv[1]).read_bytes()
digest = hashlib.sha256(source).hexdigest()
assert digest == "0a9a2e507fbc6eb66901c2653fbb16ef06fe1acaaafa6c73466168c1c82e8ca8"
lines = source.splitlines(keepends=True)
root = Path(__file__).parent


def source_ranges(numbers):
    ranges = []
    for number in numbers:
        if ranges and number == ranges[-1][1] + 1:
            ranges[-1][1] = number
        else:
            ranges.append([number, number])
    return ranges


def extract(identifiers, comments):
    selected = [
        (number, line)
        for number, line in enumerate(lines, 1)
        if line.startswith((b"EPHEDISP", b"P ", b"A ", b"T "))
        or (comments and number <= 170 and line.startswith(b"#"))
        or (line.startswith(b"S ") and line[3:11].rstrip() in identifiers)
        or (line.startswith(b"D ") and line[45:53].rstrip() in identifiers)
    ]
    sites = sum(line.startswith(b"S ") for _, line in selected)
    samples = sum(line.startswith(b"D ") for _, line in selected)
    counts = f"P T 3 S {sites:10d} E      1 D {samples:10d}\n".encode()
    output, retained, edits = [], [], []
    for extract_line, (source_line, line) in enumerate(selected, 1):
        if line.startswith(b"P "):
            edits.append(
                {
                    "source_line": source_line,
                    "extract_line": extract_line,
                    "original": line.decode(),
                    "replacement": counts.decode(),
                }
            )
            output.append(counts)
        else:
            retained.append(source_line)
            output.append(line)
    result = b"".join(output)
    return result, {
        "sha256": hashlib.sha256(result).hexdigest(),
        "bytes": len(result),
        "unchanged_source_ranges": source_ranges(retained),
        "edits": edits,
    }


positive, manifest = extract({b"00NA", b"weuc", b"yela"}, True)
(root / "imls_merra2_extract.eph").write_bytes(positive)
manifest = {
    "url": "https://massloading.net/imsl/load_list/atm/merra2/atm_merra2_20200101_0000.eph",
    "source_sha256": digest,
    "source_bytes": len(source),
    "extract_sha256": manifest.pop("sha256"),
    "extract_bytes": manifest.pop("bytes"),
    "identifiers": ["00NA", "weuc", "yela"],
    **manifest,
}
duplicates = {}
for number, line in enumerate(lines, 1):
    if line.startswith(b"S "):
        duplicates.setdefault(line[3:11].rstrip().decode(), []).append(number)
manifest["source_duplicate_site_labels"] = {
    label: numbers for label, numbers in duplicates.items() if len(numbers) > 1
}
negative, manifest["negative_duplicate_fixture"] = extract({b"AGGO"}, False)
(root / "imls_duplicate_identifier.eph").write_bytes(negative)
(root / "extract_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
