"""Freeze/check four published COST-G samples and geometric DE421 vectors offline."""

import argparse
from pathlib import Path
import hashlib
import numpy as np
import spiceypy as spice

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("reference_directory", type=Path)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
ref = args.reference_directory
assert hashlib.sha256((ref / "de421.bsp").read_bytes()).hexdigest() == (
    "08b20db2ae22488650641c5a9033e5bfda4b1c4b440cfeaf20f621cfa18ecdb3"
)


def read(name, skip):
    return np.loadtxt(ref / "costg" / name, skiprows=skip)


orbit = read("00orbit_icrf.txt", 5)
fixed = read("00orbit_itrf.txt", 5)
rotation = read("01earthRotation_rotaryMatrix.txt", 5)
eop = read("01earthRotation_interpolatedEOP.txt", 2)
arguments = read("90doodsonArguments.txt", 5)
solid = read("04solidEarthTide_icrf.txt", 5)
pole = read("05poleTide_icrf.txt", 5)
spice.furnsh(str(ref / "de421.bsp"))
# Standard NAIF DELTET constants: TDT -> ET needs no civil leap-second lookup.
spice.pdpool("DELTET/DELTA_T_A", [32.184])
spice.pdpool("DELTET/K", [1.657e-3])
spice.pdpool("DELTET/EB", [1.671e-2])
spice.pdpool("DELTET/M", [6.239996, 1.99096871e-7])
lines = [
    "# index mjdGPS tdb tt R[9 row-major] r[3] Doodson[6] pole[2 radians] Sun[3] Moon[3] solidA[3] poleA[3]"
]
try:
    for i in [0, 720, 1440, 2879]:
        gps = (orbit[i, 0] - 51544.5) * 86400
        tt = gps + 51.184
        et = spice.unitim(tt, "TDT", "TDB")
        matrix = rotation[i, 1:].reshape(3, 3)
        assert np.linalg.norm(matrix @ orbit[i, 1:4] - fixed[i, 1:4]) < 1e-7
        bodies = [
            (spice.spkpos(body, et, "J2000", "NONE", "EARTH")[0] * 1000).tolist()
            for body in ["SUN", "MOON"]
        ]
        numbers = (
            [i, orbit[i, 0], et, tt]
            + matrix.reshape(9).tolist()
            + orbit[i, 1:4].tolist()
            + arguments[i, 1:7].tolist()
            + eop[i, 1:3].tolist()
            + bodies[0]
            + bodies[1]
            + solid[i, 1:4].tolist()
            + pole[i, 1:4].tolist()
        )
        lines.append(" ".join(format(value, ".17g") for value in numbers))
finally:
    spice.kclear()
text = "\n".join(lines) + "\n"
output = Path(__file__).with_name("costg_samples.txt")
if args.check:
    assert output.read_text() == text, "Frozen COST-G/DE421 inputs changed."
else:
    output.write_text(text)
print("Four frozen samples reproduced; matrix direction agrees with published ITRF positions.")
