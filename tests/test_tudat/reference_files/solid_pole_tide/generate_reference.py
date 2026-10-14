#!/usr/bin/env python3
"""Independent scalar colatitude equations; reads fixtures, never Tudat model code."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("primary_pdf", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    directory = Path(__file__).parent
    provenance = json.loads((directory / "source.json").read_text())
    if hashlib.sha256(args.primary_pdf.read_bytes()).hexdigest() != provenance["sha256"]:
        raise ValueError("Expected the exact pinned 1 February 2018 IERS chapter PDF")
    cpp = directory.parents[1] / "src/astro/ground_stations/unitTestPoleTideDisplacement.cpp"
    block = cpp.read_text().split("const PoleTideReferenceCase cases[] = {", 1)[1].split("};", 1)[0]
    values = [float(s) for s in re.findall(r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?", block)]
    if len(values) != 8 * 11:
        raise ValueError("Expected the original eight input/ENU/Cartesian fixture rows")
    maximum = 0.0
    for index in range(8):
        year, latitude, longitude, xp, yp, *expected = values[index * 11 : (index + 1) * 11]
        # Eq.21: t is years of 365.25 days. Eq.25: mas -> arcsec, m2 sign.
        m1 = (xp - (55.0 + 1.677 * year)) / 1000.0
        m2 = -(yp - (320.5 + 3.460 * year)) / 1000.0
        theta = math.pi / 2 - math.radians(latitude)
        longitude = math.radians(longitude)
        c, s = math.cos(longitude), math.sin(longitude)
        # Eq.26 in its original south/east/radial colatitude form, mm -> m.
        south = -0.009 * math.cos(2 * theta) * (m1 * c + m2 * s)
        east = 0.009 * math.cos(theta) * (m1 * s - m2 * c)
        radial = -0.033 * math.sin(2 * theta) * (m1 * c + m2 * s)
        # Eq.27/28: R transpose; no Tudat or ENU implementation is invoked.
        ct, st = math.cos(theta), math.sin(theta)
        x = ct * c * south - s * east + st * c * radial
        y = ct * s * south + c * east + st * s * radial
        z = -st * south + ct * radial
        generated = (east, -south, radial, x, y, z)
        maximum = max(maximum, *(abs(a - b) for a, b in zip(generated, expected)))
        if not args.check:
            print(index, *(format(value, ".17g") for value in generated), sep=",")
    if maximum > provenance["tolerance_metres"]:
        raise AssertionError(f"Frozen fixture mismatch: {maximum:.17g} metres")
    print(f"8 unchanged fixtures verified; maximum scalar difference {maximum:.3g} m")


if __name__ == "__main__":
    main()
