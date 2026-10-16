"""Read frozen DE440 SPK/PCK coefficients and compare independent 80-digit evaluation.

No observation model is run and no native kernel is imported or changed.
"""

from pathlib import Path
import hashlib
import json
import math
import numpy as np
import mpmath as mp
import spiceypy as spice
from jplephem.spk import SPK
from jplephem.pck import PCK

mp.mp.dps = 80
ROOT = Path("/home/dominic/Tudat/laser-ranging-work")
KERNELS = ROOT / "data/llr/kernels"
manifest = json.loads((KERNELS / "SOURCES.json").read_text())
spice.kclear()
for record in manifest["files"]:
    path = KERNELS / record["file"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"]
    spice.furnsh(str(path))
spk = SPK.open(str(KERNELS / "de440.bsp"))
pck = PCK.open(str(KERNELS / "moon_pa_de440_200625.bpc"))


def exact(value):
    """Retain the exact binary value stored in a coefficient or exposed Time fraction."""
    return mp.mpf(float(value))


def polynomial(segment, epoch):
    """Evaluate the stored polynomial with independent Clenshaw arithmetic."""
    initial, interval, coefficients = (
        segment._data if hasattr(segment, "target") else segment._load()
    )
    index = int(mp.floor((epoch - exact(initial)) / exact(interval)))
    offset = epoch - exact(initial) - index * exact(interval)
    argument = 2 * offset / exact(interval) - 1
    result = []
    for component in range(3):
        previous, previous_previous = mp.mpf(0), mp.mpf(0)
        for coefficient in coefficients[:-1, component, index]:
            current = exact(coefficient) + 2 * argument * previous - previous_previous
            previous_previous, previous = previous, current
        result.append(
            exact(coefficients[-1, component, index]) + argument * previous - previous_previous
        )
    return mp.matrix(result)


def state_position(body, epoch):
    """Compose the actual DE440 EMB and geocentric segments before final rounding."""
    emb = polynomial(spk[0, 3], epoch)
    if body == "Earth":
        return 1000 * (emb + polynomial(spk[3, 399], epoch))
    if body == "Moon":
        return 1000 * (emb + polynomial(spk[3, 301], epoch))
    return 1000 * polynomial(spk[0, 10], epoch)


def norm(vector):
    """Return a metre-valued norm after all high-precision differences are formed."""
    return float(mp.sqrt(sum(value * value for value in vector)))


def mp_rotation(angle, axis):
    """Use the SPICE passive 3-1-3 Euler convention of the stored PCK angles."""
    c, s = mp.cos(angle), mp.sin(angle)
    if axis == 3:
        return mp.matrix([[c, s, 0], [-s, c, 0], [0, 0, 1]])
    return mp.matrix([[1, 0, 0], [0, c, s], [0, -s, c]])


summary = {
    "sources": {
        "SPK": "https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Type%202:%20Chebyshev%20(position%20only)",
        "PCK": "https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html#Type%202:%20Chebyshev%20(Angles%20only)",
    },
    "coefficient_reference": "Exact double coefficients in the hash-pinned SPK/PCK, independently evaluated with80decimaldigits; not original JPL ASCII coefficients, which are not cached.",
    "kernel_load_order": [
        spice.kdata(i, "ALL")[0].split("/")[-1] for i in range(spice.ktotal("ALL"))
    ],
    "spk_segments": [],
    "pck_segments": [],
    "cases": {},
}
for segment in spk.segments:
    initial, interval, coefficients = segment._data
    summary["spk_segments"].append(
        {
            "center": segment.center,
            "target": segment.target,
            "type": segment.data_type,
            "degree": len(coefficients) - 1,
            "interval_seconds": interval,
            "coverage_et_seconds": [segment.start_second, segment.end_second],
        }
    )
for segment in pck.segments:
    initial, interval, coefficients = segment._load()
    summary["pck_segments"].append(
        {
            "frame_class_id": segment.body,
            "inertial_frame_id": segment.frame,
            "type": segment.data_type,
            "degree": len(coefficients) - 1,
            "interval_seconds": interval,
            "coverage_et_seconds": [segment.initial_second, segment.final_second],
        }
    )

for case, directory in [
    ("apollo", "apollo-20170309-survey-stable-cubic30"),
    ("matera", "matera-20220821-stable-cubic30"),
]:
    result = []
    events = json.loads((ROOT / "results" / directory / "event_diagnostics.json").read_text())
    for event in events:
        hours, seconds = event["launch_tdb_time_split"]
        precise = mp.mpf(hours) * 3600 + exact(seconds)
        rounded = exact(float(precise))
        bodies = {}
        for body in ["Earth", "Moon", "Sun"]:
            reference = state_position(body, precise)
            rounded_reference = state_position(body, rounded)
            actual = mp.matrix(
                [
                    exact(value) * 1000
                    for value in spice.spkezr(body, float(rounded), "J2000", "NONE", "SSB")[0][:3]
                ]
            )
            bodies[body] = {
                "time_rounding_seconds": float(rounded - precise),
                "time_rounding_position_error_m": norm(rounded_reference - reference),
                "double_spice_evaluation_error_at_rounded_epoch_m": norm(
                    actual - rounded_reference
                ),
                "total_direct_error_at_precise_epoch_m": norm(actual - reference),
            }
        direct_relative = (
            spice.spkezr("Moon", float(rounded), "J2000", "NONE", "Earth")[0][:3] * 1000
        )
        subtraction = (
            spice.spkezr("Moon", float(rounded), "J2000", "NONE", "SSB")[0][:3]
            - spice.spkezr("Earth", float(rounded), "J2000", "NONE", "SSB")[0][:3]
        ) * 1000
        subtraction_metres = (
            spice.spkezr("Moon", float(rounded), "J2000", "NONE", "SSB")[0][:3] * 1000
            - spice.spkezr("Earth", float(rounded), "J2000", "NONE", "SSB")[0][:3] * 1000
        )
        relative_reference = 1000 * (
            polynomial(spk[3, 301], rounded) - polynomial(spk[3, 399], rounded)
        )
        angles = polynomial(pck.segments[0], precise)
        angles_at_double = polynomial(pck.segments[0], rounded)
        rotation = (
            mp_rotation(angles[2], 3) * mp_rotation(angles[1], 1) * mp_rotation(angles[0], 3)
        ).T
        rotation_at_double = (
            mp_rotation(angles_at_double[2], 3)
            * mp_rotation(angles_at_double[1], 1)
            * mp_rotation(angles_at_double[0], 3)
        ).T
        # A15 is only a fixed lever arm here, not a new reflector selection.
        reflector = mp.matrix([mp.mpf("1554678.104"), mp.mpf("98094.498"), mp.mpf("765005.863")])
        actual_rotation = mp.matrix(spice.pxform("MOON_PA_DE440", "J2000", float(rounded)).tolist())
        result.append(
            {
                "index": event["index"],
                "launch_tdb_split": [hours, seconds],
                "bodies": bodies,
                "same_epoch_relative_error_m": {
                    "direct_moon_earth": norm(
                        mp.matrix(direct_relative.tolist()) - relative_reference
                    ),
                    "bcrs_km_subtraction_then_si": norm(
                        mp.matrix(subtraction.tolist()) - relative_reference
                    ),
                    "bcrs_si_subtraction": norm(
                        mp.matrix(subtraction_metres.tolist()) - relative_reference
                    ),
                },
                "a15_pck_lever_arm_error_m": {
                    "time_rounding": norm((rotation_at_double - rotation) * reflector),
                    "double_pck_evaluation_at_rounded_epoch": norm(
                        (actual_rotation - rotation_at_double) * reflector
                    ),
                    "total_direct_at_precise_epoch": norm((actual_rotation - rotation) * reflector),
                },
            }
        )
    summary["cases"][case] = result
(ROOT / "logs/llr-kernel-coefficient-audit.json").write_text(json.dumps(summary, indent=2) + "\n")
for case, rows in summary["cases"].items():
    print(case)
    for body in ["Earth", "Moon", "Sun"]:
        print(
            body,
            {key: max(r["bodies"][body][key] for r in rows) for key in rows[0]["bodies"][body]},
        )
    print(
        "relative maxima",
        {
            key: max(r["same_epoch_relative_error_m"][key] for r in rows)
            for key in rows[0]["same_epoch_relative_error_m"]
        },
    )
    print(
        "PCK maxima",
        {
            key: max(r["a15_pck_lever_arm_error_m"][key] for r in rows)
            for key in rows[0]["a15_pck_lever_arm_error_m"]
        },
    )
