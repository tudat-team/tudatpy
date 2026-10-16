"""Measure existing native Time cubic tables at frozen LLR events; no residual model."""

from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import mpmath as mp
from jplephem.spk import SPK
from tudatpy import kernel
from tudatpy.astro import time_representation as time
from tudatpy.dynamics import environment_setup as setup
from tudatpy.interface import spice
from tudatpy.math import interpolators

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--kernel-sha256", required=True)
args = parser.parse_args()
assert hashlib.sha256(Path(kernel.__file__).read_bytes()).hexdigest() == args.kernel_sha256
mp.mp.dps = 80
root = Path("/home/dominic/Tudat/laser-ranging-work")
source = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-data/analysis/laser_ranging")
manifest = json.loads((source / "llr_dataset.json").read_text())
spice.clear_kernels()
for record in manifest["kernels"]:
    spice.load_kernel(str(root / "data/llr/kernels" / record["file"]))
spice.load_kernel(str(source / "icrf_spice_alias.tf"))
spk = SPK.open(str(root / "data/llr/kernels/de440.bsp"))


def polynomial(segment, epoch):
    """Evaluate exact stored binary coefficients independently at80decimaldigits."""
    initial, duration, coefficients = segment._data
    initial, duration = mp.mpf(float(initial)), mp.mpf(float(duration))
    interval = int(mp.floor((epoch - initial) / duration))
    x = 2 * (epoch - initial - interval * duration) / duration - 1
    values = []
    for axis in range(3):
        previous = previous_previous = mp.mpf(0)
        for coefficient in coefficients[:-1, axis, interval]:
            current = mp.mpf(float(coefficient)) + 2 * x * previous - previous_previous
            previous_previous, previous = previous, current
        values.append(
            mp.mpf(float(coefficients[-1, axis, interval])) + x * previous - previous_previous
        )
    return mp.matrix(values)


def epoch_mp(epoch):
    """Reconstruct the exposed Time split without adding a large double epoch."""
    hours, seconds = epoch.__getstate__()
    return mp.mpf(hours) * 3600 + mp.mpf(seconds)


def norm(vector):
    """Form high-precision position differences before converting their norm to float."""
    return float(mp.sqrt(sum(value * value for value in vector)))


report = {
    "kernel_binary_sha256": args.kernel_sha256,
    "coefficient_reference": "Independent80digit evaluation of the exact pinned DE440 SPK coefficients, not original decimal ASCII",
    "sample_selection": "Every frozen precise launch and original baseline bounce/reception diagnostic double; no residual calculation or selection",
    "cases": {},
}
for case, directory in [
    ("apollo", "apollo-20170309-survey-stable-cubic30"),
    ("matera", "matera-20220821-stable-cubic30"),
]:
    events = json.loads((root / "results" / directory / "event_diagnostics.json").read_text())
    # UTC versus TDB differs by about69seconds; these rounded bounds reproduce the script's actual tables.
    import sys

    sys.path.insert(0, str(source))
    if case == "apollo":
        import llr_residuals

        points = llr_residuals.read_points(
            root / "data/llr" / manifest["normal_points"]["file"],
            manifest["normal_points"]["sha256"],
        )
    else:
        import matera_case

        data = json.loads((source / "matera_dataset.json").read_text())
        points = matera_case.read_points(root / "data/llr" / data["normal_points"]["file"])
    start = 60 * int(float(points[0]["utc"]) // 60) - 3600
    stop = 60 * int(float(points[-1]["utc"]) // 60 + 1) + 3600
    samples = []
    for event in events:
        launch = time.Time(*event["launch_tdb_time_split"])
        samples.extend(
            [
                launch,
                time.Time(event["complete"]["link_end_epochs_tdb_s"][1]),
                time.Time(event["complete"]["link_end_epochs_tdb_s"][3]),
            ]
        )
    references = []
    for epoch in samples:
        exact = epoch_mp(epoch)
        emb = polynomial(spk[0, 3], exact)
        earth = polynomial(spk[3, 399], exact)
        moon = polynomial(spk[3, 301], exact)
        references.append(
            {
                "Earth": 1000 * (emb + earth),
                "Moon": 1000 * (emb + moon),
                "Sun": 1000 * polynomial(spk[0, 10], exact),
                "Moon-Earth": 1000 * (moon - earth),
            }
        )
    configurations = {}
    for step in [0, 60, 30, 15, 10]:
        ephemerides = {}
        for name in ["Earth", "Moon", "Sun", "Moon-Earth"]:
            body = "Moon" if name == "Moon-Earth" else name
            origin = "Earth" if name == "Moon-Earth" else "SSB"
            settings = (
                setup.ephemeris.interpolated_spice(
                    start,
                    stop,
                    step,
                    origin,
                    "ICRF",
                    interpolators.cubic_spline_interpolation(
                        boundary_interpolation=interpolators.throw_exception_at_boundary
                    ),
                    body,
                )
                if step
                else setup.ephemeris.direct_spice(origin, "ICRF", body_name_to_use=body)
            )
            ephemerides[name] = setup.create_body_ephemeris(settings, body)
        errors = {name: [] for name in ephemerides}
        for epoch, reference in zip(samples, references):
            for name, ephemeris in ephemerides.items():
                actual = mp.matrix(
                    np.asarray(ephemeris.cartesian_state(epoch)).reshape(6)[:3].tolist()
                )
                errors[name].append(norm(actual - reference[name]))
        configurations[str(step)] = {
            name: {
                "maximum_position_error_m": max(values),
                "rms_position_error_m": float(np.sqrt(np.mean(np.square(values)))),
            }
            for name, values in errors.items()
        }
    report["cases"][case] = {
        "sample_count": len(samples),
        "table_bounds_tdb_seconds": [start, stop],
        "configurations_by_step_seconds_zero_is_direct": configurations,
    }
(root / "logs/llr-native-interpolation-audit.json").write_text(json.dumps(report, indent=2) + "\n")
for case, result in report["cases"].items():
    print(
        case,
        json.dumps(result["configurations_by_step_seconds_zero_is_direct"], indent=2),
        flush=True,
    )
