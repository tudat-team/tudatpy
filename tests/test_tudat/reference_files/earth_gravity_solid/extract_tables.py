"""Reproduce the fixed IERS frequency tables and compare the primary PDF text."""

import argparse
from pathlib import Path
import re

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("reference_directory", type=Path)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
root = Path(__file__).resolve().parents[4]
output = root / "src/tudat/astro/gravitation/iers2010SolidTideFrequencyTerms.inc"
primary = (args.reference_directory / "iers_chapter6_2018.txt").read_text()


def parse(text, suffix):
    result = []
    for line in text.splitlines():
        values = line.split()
        code = next((i for i, s in enumerate(values) if re.fullmatch(r"\d{2,3},\d{3}", s)), None)
        if code is None:
            continue
        tail = values[code + 1 :]
        multipliers = tail[:6] if suffix == "a" else tail[1:7]
        ip, op = (
            (tail[-2], tail[-1])
            if suffix == "a"
            else ((tail[-3], tail[-1]) if suffix == "b" else (tail[-1], "0.0"))
        )
        result.append((tuple(map(int, multipliers)), float(ip), float(op)))
    return result


lines = [
    "// IERS Conventions 2010 Tables 6.5a/b/c, transcribed from pinned Orekit 13.1.8 resources.",
    "// See reference_files/earth_gravity_solid/README.rst for extraction, hashes and Apache notices.",
]
for suffix, name, size in [
    ("a", "diurnalTerms", 48),
    ("b", "zonalTerms", 21),
    ("c", "semidiurnalTerms", 2),
]:
    rows = parse((args.reference_directory / f"orekit-tab6.5{suffix}.txt").read_text(), suffix)
    start = primary.index(f"Table 6.5{suffix}:")
    end = (
        primary.index(f"Table 6.5{chr(ord(suffix)+1)}:", start)
        if suffix != "c"
        else primary.index("6.2.2", start)
    )
    assert len(rows) == size and parse(primary[start:end], suffix) == rows, suffix
    lines += ["", f"constexpr std::array< FrequencyTerm, {size} > {name} = {{ {{"]
    lines += [
        "    { { " + ", ".join(map(str, m)) + f" }}, {ip:.1f}, {op:.1f} }}," for m, ip, op in rows
    ]
    lines += ["} };"]
text = "\n".join(lines) + "\n"
if args.check:
    assert output.read_text() == text, "Frozen table transcription changed."
else:
    output.write_text(text)
print("48/21/2 frequency rows match primary IERS tables and frozen transcription.")
