"""Record callback states at the first unmodified HER HARPOS light-time failure."""

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from tudatpy import kernel
from tudatpy.astro import time_representation as time

root = Path("/home/dominic/Tudat/laser-ranging-work")
source = Path(
    "/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-data/analysis/laser_ranging/fixed_orbit_slr.py"
)
expected_sha = "bd8d893574bd29281188934ff5a62647e4c08494eb4bfb846c17cfc06a062968"
assert hashlib.sha256(Path(kernel.__file__).read_bytes()).hexdigest() == expected_sha
sys.path.insert(0, str(source.parent))
spec = importlib.util.spec_from_file_location("fixed_orbit_slr", source)
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)
original_state = analysis.FixedOrbit.state
original_loading = analysis.got410c_ocean_loading
state_calls = []
loading_calls = []
loading_models = {}


def traced_state(self, epoch):
    """Save exact local epoch differences and the returned double BCRS state."""
    value = np.asarray(original_state(self, epoch)).reshape(6)
    state_calls.append((epoch, value.copy()))
    return value


def traced_loading(*arguments):
    """Record original native station callback evaluations without extra queries."""
    model, metadata = original_loading(*arguments)
    site = metadata["site"]
    loading_models[site] = model

    def evaluate(epoch):
        """Keep the original reported harmonic displacement and callback epoch."""
        value = np.asarray(model(epoch)).reshape(3)
        loading_calls.append((site, float(epoch), value.copy()))
        return value

    return evaluate, metadata


analysis.FixedOrbit.state = traced_state
analysis.got410c_ocean_loading = traced_loading
sys.argv = [
    str(source),
    "--data-root",
    str(root / "data/slr"),
    "--output",
    str(root / "results/slr-ocean-harpos-first-point-trace-attempt1"),
    "--got410c-harpos",
    str(root / "research/ocean-loading/toc_got410c_harmod.hps"),
]
error = None

try:
    analysis.main()
except Exception as caught:
    error = str(caught)

first_epoch = state_calls[0][0]
trace = []
for index, (epoch, state) in enumerate(state_calls):
    trace.append(
        {
            "index": index,
            "epoch_split": epoch.__getstate__() if isinstance(epoch, time.Time) else None,
            "epoch_double_hex": float(epoch).hex(),
            "epoch_minus_first_seconds": float(epoch - first_epoch),
            "dt_previous_seconds": float(epoch - state_calls[index - 1][0]) if index else None,
            "dt_two_calls_back_seconds": (
                float(epoch - state_calls[index - 2][0]) if index > 1 else None
            ),
            "bcrs_state": state.tolist(),
            "bcrs_state_hex": [float(value).hex() for value in state],
            "position_minus_previous_m": (
                (state[:3] - state_calls[index - 1][1][:3]).tolist() if index else None
            ),
            "position_minus_two_calls_back_m": (
                (state[:3] - state_calls[index - 2][1][:3]).tolist() if index > 1 else None
            ),
        }
    )

loading_trace = [
    {"site": site, "epoch_double_hex": epoch.hex(), "displacement_xyz_m": value.tolist()}
    for site, epoch, value in loading_calls[-120:]
]
last_epoch = next(epoch for site, epoch, _ in reversed(loading_calls) if site == "Herstmon")
smoothness = []
base_displacement = loading_models["Herstmon"](last_epoch)
for offset in (-1e-3, -1e-6, -120e-9, -60e-9, 0.0, 60e-9, 120e-9, 1e-6, 1e-3):
    actual_epoch = last_epoch + offset
    delta = loading_models["Herstmon"](actual_epoch) - base_displacement
    smoothness.append(
        {
            "requested_offset_seconds": offset,
            "actual_double_offset_seconds": actual_epoch - last_epoch,
            "displacement_delta_xyz_m": delta.tolist(),
            "displacement_delta_norm_m": float(np.linalg.norm(delta)),
        }
    )

result = {
    "kernel_sha256": expected_sha,
    "analysis_script_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "diagnostic_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "convergence_unchanged": {"maximum_iterations": 50, "fractional_tolerance": 1e-12},
    "error": error,
    "callback_state_count": len(state_calls),
    "loading_callback_count": len(loading_calls),
    "last_satellite_calls": trace[-120:],
    "last_loading_calls": loading_trace,
    "loading_smoothness_at_last_Herstmon_evaluation": smoothness,
}
path = root / "logs/slr-ocean-harpos-first-point-trace-attempt1.json"
path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
print(
    json.dumps(
        {
            "error": error,
            "callback_state_count": len(state_calls),
            "loading_callback_count": len(loading_calls),
            "trace": str(path),
        },
        indent=2,
    )
)
