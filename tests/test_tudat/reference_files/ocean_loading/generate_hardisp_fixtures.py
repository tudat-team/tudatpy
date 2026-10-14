#!/usr/bin/env python3
"""Generate C++ ocean-loading fixtures by running pinned original IERS Fortran.

Fortran is optional and used only here; production and regression tests are C++.
The generated executable is renamed, and only the output format is changed to
print 12 decimals. Numerical reference routines remain unchanged. Original
notices are retained in every source file. See reference_provenance.json.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile


def generate(compiler, check=False):
    root = Path(__file__).resolve().parent
    provenance = json.loads((root / "reference_provenance.json").read_text())
    original = root / "iers_hardisp"
    for name, expected in provenance["files"].items():
        actual = hashlib.sha256((original / name).read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError(f"Reference source checksum mismatch: {name}")
    stations = json.loads((root / "hardisp_station_inputs.json").read_text())
    epochs = [
        (1999, 3, 17, 13, 24, 51),
        (2004, 9, 2, 5, 6, 7),
        (2009, 6, 25, 1, 10, 45),
        (2013, 4, 12, 19, 17, 3),
        (2017, 11, 22, 17, 3, 15),
        (2024, 8, 19, 23, 41, 29),
    ]
    cases = []
    with tempfile.TemporaryDirectory(prefix="tudat-hardisp-reference-") as temporary:
        directory = Path(temporary)
        main = (original / "HARDISP.F").read_text()
        main = main.replace("PROGRAM HARDISP", "PROGRAM TUDAT_HARDISP_REFERENCE", 1)
        main = main.replace("FORMAT(3F14.6)", "FORMAT(3F20.12)", 1)
        main_path = directory / "TUDAT_HARDISP_REFERENCE.F"
        main_path.write_text(main)
        executable = directory / "hardisp-reference"
        sources = [str(main_path)] + [
            str(original / name) for name in sorted(provenance["files"]) if name != "HARDISP.F"
        ]
        command = [
            compiler,
            "-O0",
            "-std=legacy",
            "-ffixed-line-length-none",
            "-o",
            str(executable),
            *sources,
        ]
        subprocess.run(command, check=True)
        for station_index, station in enumerate(stations):
            coefficients = station["coefficients"]
            data = "".join(
                " " + "".join(f"{value:7.5f}" for value in row) + "\n" for row in coefficients[:3]
            )
            data += "".join(
                " " + "".join(f"{value:7.1f}" for value in row) + "\n" for row in coefficients[3:]
            )
            for epoch in epochs:
                result = subprocess.run(
                    [str(executable), *(str(value) for value in epoch), "1", "3600"],
                    input=data,
                    capture_output=True,
                    text=True,
                    check=True,
                )
                values = [float(value) for value in result.stdout.split()]
                if len(values) != 3:
                    raise RuntimeError(result.stdout + result.stderr)
                cases.append((station_index, epoch, values))
        provenance["generation_command"] = (
            f"python3 generate_hardisp_fixtures.py --compiler {compiler}"
        )
        provenance["reference_compile_options"] = ["-O0", "-std=legacy", "-ffixed-line-length-none"]
        provenance["compiler_version"] = subprocess.run(
            [compiler, "--version"], check=True, text=True, capture_output=True
        ).stdout.splitlines()[0]
    lines = [
        "// Generated independently by the pinned IERS HARDISP Fortran package.",
        "// See reference_provenance.json and generate_hardisp_fixtures.py; no runtime downloads.",
        "// clang-format off",
        "struct HardispReferenceStation { const char* identifier; double longitudeDegrees; double latitudeDegrees; double coefficients[6][11]; };",
        "static const HardispReferenceStation hardispReferenceStations[] = {",
    ]
    for station in stations:
        lines.append(
            '    { "'
            + station["identifier"]
            + '", '
            + repr(station["longitude_degrees"])
            + ", "
            + repr(station["latitude_degrees"])
            + ", {"
        )
        lines.extend(
            "        { " + ", ".join(repr(value) for value in row) + " },"
            for row in station["coefficients"]
        )
        lines.append("    } },")
    lines += [
        "};",
        "struct HardispReferenceCase { int stationIndex; int year, month, day, hour, minute, second; double upSouthWest[3]; };",
        "static const HardispReferenceCase hardispReferenceCases[] = {",
    ]
    for station_index, epoch, values in cases:
        lines.append(
            "    { "
            + ", ".join(map(str, (station_index, *epoch)))
            + ", { "
            + ", ".join(f"{value:.12f}" for value in values)
            + " } },"
        )
    lines += ["};", "// clang-format on", ""]
    generated = "\n".join(lines)
    if check:
        if (root / "hardisp_reference_cases.h").read_text() != generated:
            raise RuntimeError("Frozen HARDISP fixtures differ; regenerate and review the change")
        print(f"Verified {len(cases)} frozen independent cases.")
        return
    (root / "hardisp_reference_cases.h").write_text(generated)
    provenance["fixture_cases"] = len(cases)
    provenance["fixture_output_precision_metres"] = 1e-12
    provenance["test_tolerance_metres"] = 1e-6
    provenance["tolerance_reason"] = (
        "IERS Fortran uses single-precision admittance/synthesis; C++ promotes to double. One micrometre covers numerical precision for the explicit high-load inputs, below HARDISP's stated model precision."
    )
    (root / "reference_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Generated {len(cases)} independent cases across {len(stations)} stations and 25 years.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default=os.environ.get("FC", "gfortran"))
    parser.add_argument(
        "--check", action="store_true", help="Verify frozen fixtures without modifying files"
    )
    args = parser.parse_args()
    generate(args.compiler, args.check)
