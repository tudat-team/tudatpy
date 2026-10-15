"""Reproduce degree-two TiME22 extracts and independent provider-formula references.

Usage: python regenerate.py /path/to/extracted/TiME22
No Tudat module or numerical backend is imported.
"""

import hashlib
import json
import sys
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = Path(sys.argv[1])
ARCHIVE_SHA = "79c9ef695b4c15d45f5a35dd706a7ed2e1df5acc4d7e509b5f6e404a71467ab3"
URL = "https://datapub.gfz-potsdam.de/download/10.5880.GFZ.1.3.2022.006-VENonf/2022-006_Sulzbach-et-al_time22.tar.gz"
D = Decimal
EPS = D(2) ** -52
PI = D("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628")


def sha(data):
    """Identify the exact provider/extract bytes."""
    return hashlib.sha256(data).hexdigest()


def line_ranges(numbers):
    """Record consecutive source lines compactly without per-line hashes."""
    ranges = []
    start = end = numbers[0]
    for number in numbers[1:]:
        if number == end + 1:
            end = number
        else:
            ranges.append(str(start) if start == end else f"{start}-{end}")
            start = end = number
    ranges.append(str(start) if start == end else f"{start}-{end}")
    return ",".join(ranges)


def identifier(multipliers):
    """Encode only the decimal/A/B/Z labels present in the supplied v1.1 table."""
    digits = [q + (0 if i == 0 else 5) for i, q in enumerate(multipliers)]
    text = "".join({10: "A", 11: "B", -1: "Z"}.get(q, str(q)) for q in digits)
    assert len(text) == 6
    return text[:3] + "." + text[3:]


def sincos(angle):
    """Use 80-digit Decimal Taylor sums, independent of the C++/libm evaluator."""
    angle %= 2 * PI
    if angle > PI:
        angle -= 2 * PI
    sine, cosine = angle, D(1)
    sterm, cterm = angle, D(1)
    for k in range(1, 180):
        sterm *= -(angle * angle) / D((2 * k) * (2 * k + 1))
        cterm *= -(angle * angle) / D((2 * k - 1) * (2 * k))
        sine += sterm
        cosine += cterm
    return sine, cosine


def gamma(count):
    """Conservative IEEE binary64 operation-count bound using epsilon=2^-52."""
    return count * EPS / (1 - count * EPS)


previous = ROOT / "manifest.json"
if previous.exists():
    for item in json.loads(previous.read_text())["files"]:
        assert sha((SOURCE / item["path"]).read_bytes()) == item["source_sha256"]

manifest = {
    "product": "TiME22 revision 1.1, 2025-03-22",
    "doi": "https://doi.org/10.5880/GFZ.1.3.2022.006",
    "license": "CC-BY-4.0",
    "authors": "Roman Sulzbach, Kyriakos Balidakis, Henryk Dobslaw, Maik Thomas",
    "archive_url": URL,
    "archive_bytes": 184944955,
    "archive_sha256": ARCHIVE_SHA,
    "provider_readme_sha256": "99dda6813e171fc3b1cd537ad30dbaa22ef43c71415ea19f0f27e4eccc78f206",
    "extraction": "GFC: retain selected original declaration/header lines and all original gfc rows with n<=2. No numeric line edits. Supplementary tables unchanged.",
    "files": [],
}
profiles = {}
cpp = [
    "// TiME22 revision 1.1 ordered inventories; reproduced by reference_files/time22-v1_1/regenerate.py.",
    "// clang-format off",
]
for component in ("ATM", "OCN"):
    prefix = f"TiME22_{component}"
    pair_lines = (SOURCE / component / f"{prefix}_001fileList.txt").read_text().splitlines()
    pairs = [line.split() for line in pair_lines]
    multipliers = [
        tuple(map(int, row.split()))
        for row in (SOURCE / component / f"{prefix}_002doodson.txt").read_text().splitlines()
    ]
    targets = [identifier(row) for row in multipliers]
    majors = [pair[0].split("_")[2] for pair in pairs]
    matrix = [
        [D(x) for x in row.split()]
        for row in (SOURCE / component / f"{prefix}_003admittance.txt").read_text().splitlines()
    ]
    assert len(pairs) == (16 if component == "ATM" else 57)
    assert len(targets) == (16 if component == "ATM" else 435)
    assert len(set(targets)) == len(targets)
    assert len(matrix) == len(pairs) and all(len(row) == len(targets) for row in matrix)
    declarations = {
        "modelname",
        "product_type",
        "earth_gravity_constant",
        "radius",
        "max_degree",
        "norm",
        "errors",
        "begin_of_head",
        "end_of_head",
    }
    coefficients = []
    for pair in pairs:
        fields = []
        for name in pair:
            relative = Path(component) / name
            data = (SOURCE / relative).read_bytes()
            kept = []
            terms = {}
            for number, line in enumerate(data.splitlines(keepends=True), 1):
                tokens = line.decode("ascii").split()
                if not tokens:
                    continue
                if tokens[0] == "gfc":
                    n, m = map(int, tokens[1:3])
                    if n <= 2:
                        kept.append((number, line))
                        terms[n, m] = (D(tokens[3]), D(tokens[4]))
                elif tokens[0] in declarations or any(
                    marker in line
                    for marker in (
                        b"--> ",
                        b"Converted to Stokes",
                        b"Contact:",
                        b"Full tidal potential",
                        b"Given in the Centre-of-figure",
                        b"Doodson-Warburg phase correction is already applied",
                    )
                ):
                    kept.append((number, line))
            output = b"".join(line for _, line in kept)
            destination = ROOT / relative
            destination.parent.mkdir(exist_ok=True)
            destination.write_bytes(output)
            manifest["files"].append(
                {
                    "path": str(relative),
                    "source_bytes": len(data),
                    "source_sha256": sha(data),
                    "extract_sha256": sha(output),
                    "retained_source_lines": line_ranges([i for i, _ in kept]),
                }
            )
            fields.append(terms)
        coefficients.append(fields)
    for suffix in ("001fileList", "002doodson", "003admittance"):
        relative = Path(component) / f"{prefix}_{suffix}.txt"
        data = (SOURCE / relative).read_bytes()
        (ROOT / relative).write_bytes(data)
        manifest["files"].append(
            {
                "path": str(relative),
                "source_bytes": len(data),
                "source_sha256": sha(data),
                "extract_sha256": sha(data),
                "retained_source_lines": "all; unchanged",
            }
        )
    profiles[component] = (majors, targets, multipliers, matrix, coefficients)
    label = "AtmosphericTargets" if component == "ATM" else "OceanTargets"
    cpp.append(f"const std::vector< std::string > time22{label} = {{")
    for start in range(0, len(targets), 8):
        cpp.append("    " + ", ".join('"' + x + '"' for x in targets[start : start + 8]) + ",")
    cpp.append("};")
    if component == "OCN":
        cpp.append("const std::vector< std::string > time22OceanInputs = {")
        for start in range(0, len(majors), 8):
            cpp.append("    " + ", ".join('"' + x + '"' for x in majors[start : start + 8]) + ",")
        cpp.append("};")
cpp.append("// clang-format on")
(ROOT.parents[3] / "src/tudat/simulation/environment_setup/time22Inventories.inc").write_text(
    "\n".join(cpp) + "\n"
)

# Read the exact binary64 images used by the native test, retaining the original COST-G angle file.
angle_path = ROOT.parent / "harmonic_potential/costg_90doodsonArguments.txt"
angles = [
    [D.from_float(float(x)) for x in row.split()[1:7]]
    for row in angle_path.read_text().splitlines()
]
reference = ["# component admittance sample n m C S absolute_C_budget absolute_S_budget"]
with localcontext() as context:
    context.prec = 80
    for component, (majors, targets, multipliers, matrix, coefficients) in profiles.items():
        for use_admittance in (False, True):
            phase_rows = (
                multipliers if use_admittance else [multipliers[targets.index(x)] for x in majors]
            )
            weights = (
                matrix
                if use_admittance
                else [[D(int(i == j)) for j in range(len(majors))] for i in range(len(majors))]
            )
            for sample, arguments in enumerate(angles):
                phases = [sum(D(q) * a for q, a in zip(row, arguments)) for row in phase_rows]
                trig = [sincos(x) for x in phases]
                for m in range(3):
                    values = []
                    budgets = []
                    for column in range(2):
                        # Provider formula: first contract the input-by-target map with target trigonometric phases, then with input Stokes fields.
                        value = D(0)
                        for i, pair in enumerate(coefficients):
                            cos_weight = sum(
                                weights[i][t] * trig[t][1] for t in range(len(phase_rows))
                            )
                            sin_weight = sum(
                                weights[i][t] * trig[t][0] for t in range(len(phase_rows))
                            )
                            value += (
                                pair[0][2, m][column] * cos_weight
                                + pair[1][2, m][column] * sin_weight
                            )
                        bound = D(0)
                        for t, row in enumerate(phase_rows):
                            magnitude = sum(
                                abs(weights[i][t] * pair[quadrature][2, m][column])
                                for i, pair in enumerate(coefficients)
                                for quadrature in range(2)
                            )
                            phase_error = gamma(12) * sum(
                                abs(D(q) * a) for q, a in zip(row, arguments)
                            )
                            bound += magnitude * (
                                phase_error
                                + 4 * EPS
                                + gamma(2 * len(majors) + 2)
                                + gamma(2 * len(phase_rows) + 4)
                            )
                        values.append(value)
                        budgets.append(bound)
                    reference.append(
                        f"{component} {int(use_admittance)} {sample} 2 {m} {values[0]:.35e} {values[1]:.35e} {budgets[0]:.17e} {budgets[1]:.17e}"
                    )
(ROOT / "provider-reference.txt").write_text("\n".join(reference) + "\n")
manifest["angles_sha256"] = sha(angle_path.read_bytes())
manifest["provider_reference_sha256"] = sha((ROOT / "provider-reference.txt").read_bytes())
manifest["inventory_sha256"] = sha(
    (ROOT.parents[3] / "src/tudat/simulation/environment_setup/time22Inventories.inc").read_bytes()
)
(ROOT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
