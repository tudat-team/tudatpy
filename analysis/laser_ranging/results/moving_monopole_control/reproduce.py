"""Reproduce the saved fixed-endpoint Hees2014 comparison without Tudat or new solves."""

import csv
import json
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).parent

mp.mp.dps = 80

inputs = json.loads((ROOT / "endpoint_inputs.json").read_text())

C = mp.mpf(inputs["speed_of_light_m_s"])


def dot(a, b):
    """Retain arbitrary precision during short finite-path vector operations."""

    return sum(x * y for x, y in zip(a, b))


def norm(a):
    """Evaluate an ordinary Euclidean vector norm."""

    return mp.sqrt(dot(a, a))


def vec(values):
    """Represent the exact supplied binary doubles at the selected decimal precision."""

    return [mp.mpf(float(v)) for v in values]


def static(xa, xb, centre, gm):
    """Evaluate the native static monopole logarithm as a length, gammaPPN=1."""

    a = norm([x - y for x, y in zip(xa, centre)])

    b = norm([x - y for x, y in zip(xb, centre)])

    r = norm([x - y for x, y in zip(xb, xa)])

    return 2 * gm / C**2 * mp.log((a + b + r) / (a + b - r))


def moving(xa, xb, centre, velocity, tb_minus_t0, gm, quadrature=False):
    """Use Hees Eqs26/30/45 and, independently, the null-segment potential integral23."""

    displacement = [b - a for a, b in zip(xa, xb)]

    r = norm(displacement)

    n = [x / r for x in displacement]

    beta = [v / C for v in velocity]

    lorentz = 1 / mp.sqrt(1 - dot(beta, beta))

    factor = lorentz * (1 - dot(beta, n))

    k = lorentz**2 / (1 + lorentz)

    vectors = []

    for x in [xa, xb]:
        relative = [a - b for a, b in zip(x, centre)]

        vectors.append(
            [
                relative[i]
                + k * beta[i] * dot(beta, relative)
                - lorentz * velocity[i] * tb_minus_t0
                for i in range(3)
            ]
        )
    transformed_a = [vectors[0][i] + lorentz * beta[i] * r for i in range(3)]

    transformed_b = vectors[1]

    ar = norm(transformed_a)

    br = norm(transformed_b)

    transformed_r = r * factor

    result = 2 * gm / C**2 * factor * mp.log((ar + br + transformed_r) / (ar + br - transformed_r))

    integral_error = None

    if quadrature:
        g = [b - a for a, b in zip(transformed_a, transformed_b)]

        integral = mp.quad(
            lambda s: 1 / norm([transformed_b[i] - s * g[i] for i in range(3)]),
            [0, mp.mpf(".01"), mp.mpf(".1"), mp.mpf(".5"), mp.mpf(".9"), mp.mpf(".99"), 1],
        )

        direct = 2 * r * gm / C**2 * factor**2 * integral

        integral_error = direct - result

    return result, factor, integral_error


expected = list(csv.DictReader((ROOT / "body_delay_comparison.csv").open()))

errors = []

for leg in inputs["legs"]:
    xa, xb = [vec(x) for x in leg["endpoint_bcrs_m"]]

    for name, state in leg["body_states"].items():
        centre = vec(state["body_reference_bcrs_state"][:3])

        velocity = vec(state["body_reference_bcrs_state"][3:])

        delta_t = mp.mpf(leg["endpoint_epochs_tdb_s"][1]) - mp.mpf(
            state["body_reference_epoch_tdb_s"]
        )

        gm = mp.mpf(state["gm_si"])

        frozen = static(xa, xb, centre, gm)

        value, _, error = moving(xa, xb, centre, velocity, delta_t, gm, True)

        row = next(
            r
            for r in expected
            if r["station"] == leg["station"]
            and int(r["source_line"]) == leg["source_line"]
            and r["leg"] == leg["leg"]
            and r["body"] == name
        )

        errors.append(abs(float(value - frozen) - float(row["moving_minus_static_m"])))
        print(
            leg["station"],
            leg["source_line"],
            leg["leg"],
            name,
            "moving_minus_static_m",
            mp.nstr(value - frozen, 20),
            "quadrature_error_m",
            mp.nstr(error, 5),
        )
print("max_difference_from_recorded_binary_double_m", max(errors))
