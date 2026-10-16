"""Compare old/stable cubic evaluation with an independent 80-digit natural spline."""

from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from tudatpy import kernel
from tudatpy.astro import time_representation as time
from tudatpy.data_input import resource_paths
from tudatpy.dynamics import environment_setup
from tudatpy.dynamics.environment_setup import ephemeris
from tudatpy.interface import spice
from tudatpy.math import interpolators

work = Path("/home/dominic/Tudat/laser-ranging-work")
assert (
    hashlib.sha256(Path(kernel.__file__).read_bytes()).hexdigest()
    == "bd8d893574bd29281188934ff5a62647e4c08494eb4bfb846c17cfc06a062968"
)
kernel_dir = Path(resource_paths.get_spice_kernel_path())
spice.clear_kernels()
for name in (
    "pck00010.tpc",
    "inpop19a_TDB_m100_p100_spice.tpc",
    "inpop19a_TDB_m100_p100_spice.bsp",
    "naif0012.tls",
):
    spice.load_kernel(str(kernel_dir / name))
start, end, step = 765028800, 765806400, 30
settings = ephemeris.interpolated_spice(
    start,
    end,
    step,
    "SSB",
    "J2000",
    interpolators.cubic_spline_interpolation(
        boundary_interpolation=interpolators.throw_exception_at_boundary
    ),
    "Earth",
)
native = environment_setup.create_body_ephemeris(settings, "Earth")
epochs = list(range(start, end + 1, step))
knots = np.array(
    [np.asarray(native.cartesian_state(time.Time(epoch))).reshape(6) for epoch in epochs]
)
assert all(
    np.array_equal(
        knots[i],
        spice.get_body_cartesian_state_at_epoch("Earth", "SSB", "J2000", "NONE", epochs[i]),
    )
    for i in range(0, len(epochs), 997)
)
trace = json.loads((work / "logs/slr-ocean-harpos-first-point-compact-attempt3.json").read_text())
query_rows = trace["last_receiver_epochs"][-2:]
query_rows = sorted(query_rows, key=lambda row: row["epoch_fraction_seconds"])
whole = query_rows[0]["epoch_whole_seconds"]
interval = (whole - start) // step

# Reproduce the existing Thomas solve with long-double scalar factors and double matrix elements.
count = len(knots) - 2
curvature_native = np.zeros_like(knots)
rhs = 6.0 * ((knots[2:] - knots[1:-1]) / 30.0 - (knots[1:-1] - knots[:-2]) / 30.0)
intermediate = np.zeros(count, dtype=np.longdouble)
solution = np.zeros((count, 6))
scale = np.longdouble(120)
solution[0] = rhs[0] / float(scale)
for j in range(1, count):
    intermediate[j] = np.longdouble(30) / scale
    scale = np.longdouble(120) - np.longdouble(30) * intermediate[j]
    solution[j] = (rhs[j] - 30.0 * solution[j - 1]) / float(scale)
for j in range(count - 2, -1, -1):
    solution[j] -= float(intermediate[j + 1]) * solution[j + 1]
curvature_native[1:-1] = solution

high_precision = []
with localcontext() as context:
    context.prec = 80
    h = Decimal(30)
    for component in range(6):
        values = [Decimal.from_float(float(value)) for value in knots[:, component]]
        source = [
            Decimal(6) * ((values[i + 2] - values[i + 1]) / h - (values[i + 1] - values[i]) / h)
            for i in range(count)
        ]
        factors = [Decimal(0)] * count
        solved = [Decimal(0)] * count
        scale = Decimal(120)
        solved[0] = source[0] / scale
        for j in range(1, count):
            factors[j] = h / scale
            scale = Decimal(120) - h * factors[j]
            solved[j] = (source[j] - h * solved[j - 1]) / scale
        for j in range(count - 2, -1, -1):
            solved[j] -= factors[j + 1] * solved[j + 1]
        curvatures = [Decimal(0), *solved, Decimal(0)]
        component_results = []
        for row in query_rows:
            x = Decimal(whole - epochs[interval]) + Decimal.from_float(
                row["epoch_fraction_seconds"]
            )
            b = x / h
            a = Decimal(1) - b
            c = (a * a * a - a) * h * h / Decimal(6)
            d = (b * b * b - b) * h * h / Decimal(6)
            value = (
                a * values[interval]
                + b * values[interval + 1]
                + c * curvatures[interval]
                + d * curvatures[interval + 1]
            )
            rounded = float(value)
            lower = (
                Decimal.from_float(math.nextafter(rounded, -math.inf)) + Decimal.from_float(rounded)
            ) / 2
            upper = (
                Decimal.from_float(math.nextafter(rounded, math.inf)) + Decimal.from_float(rounded)
            ) / 2
            component_results.append(
                {
                    "value_80_digits": str(value),
                    "correctly_rounded_double": rounded,
                    "correctly_rounded_hex": rounded.hex(),
                    "distance_to_nearest_double_rounding_midpoint": str(
                        min(value - lower, upper - value)
                    ),
                }
            )
        high_precision.append(
            {
                "component": component,
                "left_curvature_80_digits": str(curvatures[interval]),
                "right_curvature_80_digits": str(curvatures[interval + 1]),
                "evaluations": component_results,
                "high_precision_delta": str(
                    Decimal(component_results[1]["value_80_digits"])
                    - Decimal(component_results[0]["value_80_digits"])
                ),
            }
        )

evaluations = []
for row in query_rows:
    fraction = row["epoch_fraction_seconds"]
    epoch = time.Time(whole) + fraction
    actual = np.asarray(native.cartesian_state(epoch)).reshape(6)
    assert np.array_equal(actual[:3], row["earth_bcrs_position_m"])
    a = (np.longdouble(epochs[interval + 1] - whole) - np.longdouble(fraction)) / np.longdouble(30)
    b = np.longdouble(1) - a
    c = (a * a * a - a) / np.longdouble(6) * np.longdouble(900)
    d = (b * b * b - b) / np.longdouble(6) * np.longdouble(900)
    a, b, c, d = map(float, (a, b, c, d))
    y0, y1 = knots[interval : interval + 2]
    m0, m1 = curvature_native[interval : interval + 2]
    old = ((a * y0 + b * y1) + c * m0) + d * m1
    stable = y0 + ((b * (y1 - y0) + c * m0) + d * m1)
    evaluations.append(
        {
            "epoch_fraction_seconds": fraction,
            "actual_native_state": actual.tolist(),
            "actual_native_hex": [float(x).hex() for x in actual],
            "reconstructed_old_state": old.tolist(),
            "reconstructed_old_matches_native_every_component": np.array_equal(old, actual),
            "stable_rearranged_state": stable.tolist(),
            "stable_rearranged_hex": [float(x).hex() for x in stable],
            "stable_minus_correctly_rounded_m": [
                (
                    float(stable[i])
                    - high_precision[i]["evaluations"][len(evaluations)]["correctly_rounded_double"]
                )
                for i in range(6)
            ],
        }
    )

result = {
    "kind": "Read-only same-polynomial evaluation; no production replacement or convergence change",
    "kernel_sha256": hashlib.sha256(Path(kernel.__file__).read_bytes()).hexdigest(),
    "diagnostic_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "knots_sha256_float64_native_bytes": hashlib.sha256(knots.tobytes()).hexdigest(),
    "natural_spline": {
        "start": start,
        "end": end,
        "step": step,
        "node_count": len(knots),
        "curvature_end_conditions": "zero at both endpoints",
        "independent_decimal_precision": 80,
        "exact_native_knots_match_direct_spice_selected_nodes": True,
        "interval_lower": epochs[interval],
        "interval_upper": epochs[interval + 1],
        "left_knot": knots[interval].tolist(),
        "right_knot": knots[interval + 1].tolist(),
        "left_curvature_reproduced_double": curvature_native[interval].tolist(),
        "right_curvature_reproduced_double": curvature_native[interval + 1].tolist(),
    },
    "evaluations": evaluations,
    "high_precision": high_precision,
    "high_precision_both_epochs_round_to_same_position": all(
        hp["evaluations"][0]["correctly_rounded_hex"]
        == hp["evaluations"][1]["correctly_rounded_hex"]
        for hp in high_precision[:3]
    ),
    "stable_both_epochs_same_position": evaluations[0]["stable_rearranged_hex"][:3]
    == evaluations[1]["stable_rearranged_hex"][:3],
}
path = work / "logs/slr-cubic-same-polynomial-probe-attempt1.json"
path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
print(json.dumps(result, indent=2))
