#!/usr/bin/env python3
"""Execute pinned, unmodified IERS references; never import or run Tudat."""

import argparse
import csv
import datetime
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import urllib.request

HERE = Path(__file__).resolve().parent


def header(rows):
    lines = [
        "// Generated from unmodified IERS software by generate.py.",
        "// clang-format off",
        "#ifndef TUDAT_SHORT_PERIOD_REFERENCE_H",
        "#define TUDAT_SHORT_PERIOD_REFERENCE_H",
        "#include <array>",
        "namespace short_period_reference {",
        "struct Fixture { double mjdUtc, ttSeconds; std::array<double,2> polarLibration; double ut1LibrationUs, lodLibrationUs; std::array<double,3> ortho; std::array<double,4> ocean; std::array<double,6> arguments; };",
        "static const std::array<Fixture, %d> fixtures = {{" % len(rows),
    ]
    for row in rows:
        values = list(row.values())
        lines.append(
            "    { %s, %s, {{ %s }}, %s, %s, {{ %s }}, {{ %s }}, {{ %s }} },"
            % (
                values[0],
                values[1],
                ", ".join(values[2:4]),
                values[4],
                values[5],
                ", ".join(values[6:9]),
                ", ".join(values[9:13]),
                ", ".join(values[13:19]),
            )
        )
    lines += ["}};", "}", "#endif", "// clang-format on", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    cache = args.reference_dir.resolve()
    repository = HERE.parents[5]
    if cache == repository or repository in cache.parents:
        parser.error("Reference downloads/builds must stay outside the repository")
    cache.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((HERE / "sources.json").read_text())
    for name, source in manifest["files"].items():
        file = cache / name
        if not file.exists():
            if not args.download:
                parser.error(str(file) + " is missing; use --download")
            with urllib.request.urlopen(source["url"], timeout=60) as response:
                file.write_bytes(response.read())
        if hashlib.sha256(file.read_bytes()).hexdigest() != source["sha256"]:
            raise RuntimeError("Reference checksum mismatch: " + name)
    sources = [str(cache / name) for name in manifest["files"] if name.endswith((".F", ".f"))]
    executable = cache / "reference_driver"
    subprocess.run(
        [
            args.compiler,
            "-O0",
            "-ffixed-line-length-none",
            "-o",
            str(executable),
            str(HERE / "reference_driver.f90"),
            *sources,
        ],
        check=True,
    )
    epochs = [row["utc_mjd"] for row in csv.DictReader((HERE / "epochs.csv").open())]
    calendar = [
        datetime.date(1858, 11, 17) + datetime.timedelta(days=int(float(epoch))) for epoch in epochs
    ]
    stdin = (
        str(len(epochs))
        + "\n"
        + "\n".join(
            f"{epoch} {day.year} {day.month} {day.day}" for epoch, day in zip(epochs, calendar)
        )
        + "\n"
    )
    result = subprocess.run(
        [str(executable)], input=stdin, capture_output=True, text=True, check=True
    )
    cols = [
        "utc_mjd",
        "tt_seconds",
        "pm_x_uas",
        "pm_y_uas",
        "libration_ut1_us",
        "libration_lod_us",
        "ortho_x_uas",
        "ortho_y_uas",
        "ortho_ut1_us",
        "ocean_x_arcsec",
        "ocean_y_arcsec",
        "ocean_ut1_s",
        "ocean_lod_s",
        "gmst_plus_pi",
        "l",
        "lp",
        "f",
        "d",
        "om",
    ]
    output = result.stdout.splitlines()
    if len(output) != len(epochs):
        raise RuntimeError("Unexpected reference output row count")
    rows = []
    for epoch, line in zip(epochs, output):
        values = line.split()
        if len(values) != len(cols) - 1:
            raise RuntimeError("Unexpected reference output column count")
        rows.append(dict(zip(cols, [epoch, *values])))
    checks = [
        (0, "libration_ut1_us", 2.441143834386761746),
        (0, "libration_lod_us", -14.78971247349449492),
        (1, "ortho_x_uas", -162.8386373279636530),
        (1, "ortho_ut1_us", -23.39092370609808214),
        (2, "pm_x_uas", 24.83144238273364834),
        (3, "libration_ut1_us", -2.655705844335680244),
    ]
    for index, column, expected in checks:
        if abs(float(rows[index][column]) - expected) > 2e-11:
            raise RuntimeError("Published reference benchmark failed: " + column)
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=cols, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    artifacts = {"expected.csv": buffer.getvalue(), "fixtures.h": header(rows)}
    # Coefficients for an independent UT1/LOD generic-engine check. Main UT1 tests
    # use production resources and explicitly account for the published Table5.1b
    # versus UTLIBR.F discrepancies, not these derived tables.
    text = (cache / "UTLIBR.F").read_text()
    block = text[text.index("     ./2, -2") : text.index("* Compute the harmonic")]
    coefficients = []
    for line in block.splitlines():
        numbers = re.findall(r"[-+]?\d+(?:\.\d*)?", line[7:])
        if len(numbers) == 11:
            coefficients.append([float(n) for n in numbers])
    if len(coefficients) != 11:
        raise RuntimeError("Unexpected UTLIBR source table shape")
    artifacts["libration_ut1_lod_amplitudes.txt"] = "".join(
        " ".join(str(v) for v in r[7:11]) + "\n" for r in coefficients
    )
    artifacts["libration_ut1_lod_multipliers.txt"] = "".join(
        " ".join(str(int(v)) for v in r[:6]) + "\n" for r in coefficients
    )
    for name, text in artifacts.items():
        if args.check:
            if (HERE / name).read_text() != text:
                raise RuntimeError("Reference artifact differs: " + name)
        else:
            (HERE / name).write_text(text)
    print(("Verified" if args.check else "Generated") + f" {len(rows)} independent EOP fixtures")


if __name__ == "__main__":
    main()
