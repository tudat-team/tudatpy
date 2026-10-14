"""Independent scalar oracle; takes the original NASA IVS ASCII supplement.

No Tudat import or production parser is used. Phases follow the IERS 2003
Delaunay polynomials and SOFA GMST06/ERA00 equations; rates differentiate
those polynomials analytically along dTT/dt=1, dUT1/dt=0.99999998.
"""

import hashlib
import math
from decimal import Decimal as D, getcontext
from pathlib import Path
import sys

getcontext().prec = 50
source = Path(sys.argv[1]).read_bytes()
assert (
    hashlib.sha256(source).hexdigest()
    == "53398e40b70ccfa411d292e9a7f752a01cbac1be07d4e7fb18906d6481e57dd8"
)
terms = []
for line in source.decode().splitlines():
    if line.startswith("#") or not line.strip():
        continue
    fields = line.split()
    assert len(fields) == 17
    terms.append(([int(x) for x in fields[1:7]], [float(x) for x in fields[9:17]]))
assert len(terms) == 159
pi = D("3.1415926535897932384626433832795028841971693993751")
arcsec = pi / D(648000)
century = D(36525 * 86400)
polynomials = [
    ["485868.249036", "1717915923.2178", "31.8792", ".051635", "-.00024470"],
    ["1287104.793048", "129596581.0481", "-.5532", ".000136", "-.00001149"],
    ["335779.526232", "1739527262.8478", "-12.7512", "-.001037", ".00000417"],
    ["1072260.703692", "1602961601.2090", "-6.3706", ".006593", "-.00003169"],
    ["450160.398036", "-6962890.5431", "7.4722", ".007702", "-.00005939"],
]
gmst = [".014506", "4612.156534", "1.3915817", "-.00000044", "-.000029956", "-.0000000368"]


def polynomial(coefficients, t):
    return sum(D(c) * (t ** i if i else D(1)) for i, c in enumerate(coefficients))


def derivative(coefficients, t):
    return sum(
        i * D(c) * (t ** (i - 1) if i > 1 else D(1)) for i, c in enumerate(coefficients) if i
    )


def evaluate(tt):
    ut1 = tt - D("69.584")
    t = tt / century
    earth_rate = 2 * pi * D("1.00273781191135448") / D(86400)
    gamma = 2 * pi * (D(".7790572732640") + D("1.00273781191135448") * ut1 / D(86400))
    gamma += polynomial(gmst, t) * arcsec + pi
    arguments = [gamma % (2 * pi)] + [(polynomial(p, t) * arcsec) % (2 * pi) for p in polynomials]
    rates = [earth_rate * D(".99999998") + derivative(gmst, t) * arcsec / century]
    rates += [derivative(p, t) * arcsec / century for p in polynomials]
    values, slopes = [], []
    for component in range(4):
        contributions, tangents = [], []
        for multipliers, coefficients in terms:
            phase = float(sum(m * a for m, a in zip(multipliers, arguments)))
            velocity = float(sum(m * r for m, r in zip(multipliers, rates)))
            s, c = coefficients[2 * component : 2 * component + 2]
            contributions.append(s * math.sin(phase) + c * math.cos(phase))
            tangents.append((s * math.cos(phase) - c * math.sin(phase)) * velocity)
        values.append(math.fsum(contributions))
        slopes.append(math.fsum(tangents))
    return float(tt), float(ut1), [float(a) for a in arguments], values, slopes


epochs = [
    D(x)
    for x in [
        "-1199167165.96986127",
        "-631186508.816",
        "-384004744.816",
        "0",
        "64.184",
        "10864.184",
        "43264.184",
        "75764.184",
        "94683664.184",
        "284061666.184",
        "536500869.184",
        "586982229.184",
        "662752869.184",
        "788983269.184",
    ]
]
output = [
    "// Generated independently from the pinned NASA IVS supplement; see README.md.",
    "// clang-format off",
    "#ifndef TUDAT_DESAI_SIBOIS_REFERENCE_H",
    "#define TUDAT_DESAI_SIBOIS_REFERENCE_H",
    "#include <array>",
    "namespace desai_sibois_reference {",
    "struct Fixture { double tt, ut1; std::array<double,6> arguments; std::array<double,4> values, rates; };",
    f"static const std::array<Fixture, {len(epochs)}> fixtures = {{{{",
]
format_array = lambda values: "{{ " + ", ".join(f"{x:.17e}" for x in values) + " }}"
for epoch in epochs:
    tt, ut1, arguments, values, rates = evaluate(epoch)
    output.append(
        f"    {{ {tt:.17e}, {ut1:.17e}, {format_array(arguments)}, {format_array(values)}, {format_array(rates)} }},"
    )
output += ["}};", "}", "#endif", "// clang-format on", ""]
Path(__file__).with_name("expected.h").write_text("\n".join(output))
