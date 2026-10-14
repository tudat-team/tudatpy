#!/usr/bin/env python3
"""Original 70-digit IERS 11.18/11.19 fixture calculation; no Tudat imports."""

import argparse
import json
from decimal import Decimal, localcontext
from pathlib import Path

FROZEN = {
    "unscaled": (
        "-0.05927009731058432782438507239103079827",
        "0.00067500388626950807702023348174785247",
        "-0.05176757493587349398402470505696181319",
    ),
    "tt_to_tdb": (
        "-0.11962235703801889715288507239103079827",
        "-0.01312412407358704507197976651825214753",
        "-0.12287051583384520679542470505696181319",
    ),
    "rate": (
        "-7.23953384523679299906582658942596120526e-9",
        "-1.65526948537310939009488779777092310547e-9",
        "-8.52912798048462520908934078398968329889e-9",
    ),
}


def generate():
    with localcontext() as context:
        context.prec = 70
        r = tuple(map(Decimal, ("4075578.385", "931852.890", "4801570.154")))
        v = tuple(map(Decimal, ("27000", "-14000", "6200")))
        c2 = Decimal("299792458") ** 2
        solar_gm, lunar_gm = Decimal("1.32712440018e20"), Decimal("4.9028e12")
        solar_distance, lunar_distance = Decimal("149597870700"), Decimal("384400000")
        potential = solar_gm / solar_distance + lunar_gm / lunar_distance
        potential_rate = solar_gm * v[0] / solar_distance**2 + lunar_gm * v[1] / lunar_distance**2
        dot = sum(a * b for a, b in zip(r, v))
        result = {}
        for name, scale in (("unscaled", "0"), ("tt_to_tdb", "1.48082686741e-8")):
            result[name] = tuple(
                -(potential / c2 + Decimal(scale)) * x - y * dot / (2 * c2) for x, y in zip(r, v)
            )
        result["rate"] = tuple(-potential_rate * x / c2 for x in r)
        for name, expected in FROZEN.items():
            for generated, frozen in zip(result[name], expected):
                target = Decimal(frozen)
                if generated.quantize(Decimal(1).scaleb(target.as_tuple().exponent)) != target:
                    raise AssertionError(f"Frozen {name} value changed: {frozen} != {generated}")
        return {name: tuple(str(value) for value in values) for name, values in result.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true", help="Check frozen values and their C++ fixture literals."
    )
    args = parser.parse_args()
    generated = generate()
    if args.check:
        source = (
            Path(__file__).resolve().parents[2]
            / "src/astro/ground_stations/unitTestRelativisticStationPosition.cpp"
        )
        text = source.read_text()
        for values in FROZEN.values():
            for literal in values:
                if literal not in text:
                    raise AssertionError(f"C++ fixture literal changed or disappeared: {literal}")
        print("All 9 frozen IERS spatial correction/rate literals reproduce at 70-digit precision.")
    else:
        print(json.dumps(generated, indent=2))


if __name__ == "__main__":
    main()
