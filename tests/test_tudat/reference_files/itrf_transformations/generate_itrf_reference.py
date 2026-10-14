#!/usr/bin/env python3
"""Independent 80-digit scalar evaluation of pinned official linearized affine equations."""

from decimal import Decimal as D, getcontext
from pathlib import Path
import argparse
import hashlib
import json
import re

getcontext().prec = 80
YEAR = D(31557600)
PI = D("3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986")
INPUT = [
    D(x) for x in ["4075578.385", "931852.890", "4801570.154", ".012345", "-.023456", ".034567"]
]
EXPECTED = {
    2014: "6cc1c068d0d7e9ce53277752795dc2ef3ea5f8890c4b5d1878e1e9ccc05a0543",
    2020: "0d8dc872e7ef994966063e614adb5b1b0fe20582359c6e77f1d83259d1ff1a68",
}


def matvec(a, x):
    return [sum((a[i][j] * x[j] for j in range(3)), D(0)) for i in range(3)]


def solve(a, b):
    # Independent scalar Gaussian elimination; no Eigen or production inverse expression.
    work = [list(row) + [value] for row, value in zip(a, b)]
    for col in range(3):
        pivot = max(range(col, 3), key=lambda r: abs(work[r][col]))
        work[col], work[pivot] = work[pivot], work[col]
        divisor = work[col][col]
        work[col] = [value / divisor for value in work[col]]
        for row in range(3):
            if row != col:
                multiplier = work[row][col]
                work[row] = [x - multiplier * y for x, y in zip(work[row], work[col])]
    return [row[-1] for row in work]


def reference(values, rates, delta):
    factors = [D(".001")] * 3 + [D("1e-9")] + [PI / D(648000000)] * 3
    parameters = [(x + r * delta) * f for x, r, f in zip(values, rates, factors)]
    physical_rates = [r * f * D(".75") / YEAR for r, f in zip(rates, factors)]
    tx, ty, tz, s, rx, ry, rz = parameters
    dtx, dty, dtz, ds, drx, dry, drz = physical_rates
    a = [[1 + s, -rz, ry], [rz, 1 + s, -rx], [-ry, rx, 1 + s]]
    da = [[ds, -drz, dry], [drz, ds, -drx], [-dry, drx, ds]]
    x, v = INPUT[:3], INPUT[3:]
    forward = [t + y for t, y in zip([tx, ty, tz], matvec(a, x))]
    forward += [t + q + w for t, q, w in zip([dtx, dty, dtz], matvec(da, x), matvec(a, v))]
    inverse = solve(a, [y - t for y, t in zip(x, [tx, ty, tz])])
    inverse += solve(a, [w - t - q for w, t, q in zip(v, [dtx, dty, dtz], matvec(da, inverse))])
    return forward, inverse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    lines = [
        "// Independent 80-digit Decimal scalar reference of published IGN equations.",
        "// clang-format off",
        "struct ItrfReferenceFixture { int catalogue; const char* frame; double year; double yearRate; double forward[6]; double inverse[6]; };",
        "static const double itrfReferenceInputState[] = { "
        + ", ".join(str(float(x)) for x in INPUT)
        + " };",
        "static const ItrfReferenceFixture itrfReferenceFixtures[] = {",
    ]
    count = 0
    for catalogue, name in [
        (2014, "Transfo-ITRF2014_ITRFs.txt"),
        (2020, "Transfo-ITRF2020_TRFs.txt"),
    ]:
        path = root / name
        if hashlib.sha256(path.read_bytes()).hexdigest() != EXPECTED[catalogue]:
            raise RuntimeError("Pinned official input checksum mismatch: " + name)
        pinned_source = root.parents[3] / "src/tudat/io/pinnedItrfTransformationTables.h"
        embedded = re.search(
            r"pinnedItrf" + str(catalogue) + r'Table\[\] = R"ITRF\((.*?)\)ITRF";',
            pinned_source.read_text(),
            re.S,
        )
        if (
            embedded is None
            or hashlib.sha256(embedded.group(1).encode()).hexdigest() != EXPECTED[catalogue]
        ):
            raise RuntimeError("Production pinned raw input differs: " + name)
        pending = None
        for line in path.read_text().splitlines():
            parts = line.split()
            if parts and parts[0].startswith("ITRF") and len(parts) == 9:
                pending = (parts[0], [D(x) for x in parts[1:8]], D(parts[8]))
            elif parts and parts[0] == "rates" and len(parts) == 8:
                frame, values, epoch = pending
                rates = [D(x) for x in parts[1:]]
                for delta in [D("-10.25"), D(0), D("10.25")]:
                    forward, inverse = reference(values, rates, delta)
                    lines.append(
                        "    { "
                        + str(catalogue)
                        + ', "'
                        + frame
                        + '", '
                        + str(float(epoch + delta))
                        + ", "
                        + str(float(D(".75") / YEAR))
                        + ", { "
                        + ", ".join(str(float(x)) for x in forward)
                        + " }, { "
                        + ", ".join(str(float(x)) for x in inverse)
                        + " } },"
                    )
                    count += 1
    lines += ["};", "// clang-format on", ""]
    output = "\n".join(lines)
    target = root / "itrf_reference.h"
    if args.check:
        if target.read_text() != output:
            raise RuntimeError("Frozen independent ITRF fixtures differ")
    else:
        target.write_text(output)
    print(
        ("Verified" if args.check else "Generated")
        + f" {count} independent high-precision ITRF fixtures"
    )


if __name__ == "__main__":
    main()
