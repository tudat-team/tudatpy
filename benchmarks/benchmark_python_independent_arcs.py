"""Strong scaling through the public native-arc API; JSONL on stdout or --output.

The default is a short smoke test. --preset representative explicitly selects
8 arcs of 2 days at degree/order 50 with three repetitions. No speedup threshold
is a correctness assertion. See PYTHON_SCALING.md for timing/provenance limits.
"""

import argparse
from contextlib import nullcontext
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
from time import perf_counter

import numpy as np

REPOSITORY = Path(__file__).resolve().parents[1]
ATOL = 1e-14
RTOL = 1e-11
PRESETS = {
    "smoke": dict(arcs=8, degree=20, days_per_arc=600 / 86400, repeats=1),
    "representative": dict(arcs=8, degree=50, days_per_arc=2.0, repeats=3),
}


@dataclass(frozen=True)
class Configuration:
    preset: str
    arcs: int
    degree: int
    order: int
    days_per_arc: float
    step_size: float
    workers: tuple
    repeats: int

    @property
    def duration(self):
        return self.days_per_arc * 86400.0


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preset", choices=PRESETS, default="smoke")
    parser.add_argument("--arcs", type=int)
    parser.add_argument("--degree", type=int)
    parser.add_argument("--order", type=int, help="defaults to the selected degree")
    parser.add_argument("--days-per-arc", type=float)
    parser.add_argument("--step-size", type=float, default=30.0)
    parser.add_argument("--workers", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--repeats", type=int)
    parser.add_argument("--output", type=Path, help="new JSONL file; never overwrite")
    args = parser.parse_args(argv)
    values = PRESETS[args.preset].copy()
    for name in values:
        if getattr(args, name) is not None:
            values[name] = getattr(args, name)
    config = Configuration(
        preset=args.preset,
        **values,
        order=values["degree"] if args.order is None else args.order,
        step_size=args.step_size,
        workers=tuple(args.workers),
    )
    if config.arcs < 1 or config.repeats < 1:
        parser.error("arcs and repeats must be positive")
    if not 2 <= config.degree <= 120 or not 0 <= config.order <= config.degree:
        parser.error("require 2 <= degree <= 120 and 0 <= order <= degree")
    if (
        not math.isfinite(config.duration)
        or not math.isfinite(config.step_size)
        or config.duration <= 0
        or config.step_size <= 0
        or config.duration / config.step_size < 8
    ):
        parser.error("finite positive duration/step with at least eight steps required")
    if (
        1 not in config.workers
        or len(set(config.workers)) != len(config.workers)
        or min(config.workers) < 1
        or max(config.workers) > config.arcs
    ):
        parser.error("workers must be unique, include 1, and lie between 1 and arcs")
    return config, args.output


def command_output(arguments, cwd=REPOSITORY):
    try:
        return subprocess.check_output(
            arguments, cwd=cwd, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def collect_metadata(config, starts, states, model):
    import tudatpy
    import tudatpy.kernel
    from tudatpy.data_input import resource_paths

    module = Path(tudatpy.kernel.__file__).resolve()
    build = next((p for p in module.parents if (p / "CMakeCache.txt").is_file()), None)
    cache = {}
    if build:
        for line in (build / "CMakeCache.txt").read_text().splitlines():
            if not line.startswith(("#", "//")) and ":" in line and "=" in line:
                key_type, value = line.split("=", 1)
                key = key_type.split(":", 1)[0]
                if (
                    key.startswith(("TUDAT_", "CMAKE_CXX_", "CEREAL_", "Cereal_"))
                    or key == "CMAKE_BUILD_TYPE"
                ):
                    cache[key] = value
    distributions = {}
    for name in ("numpy", "pybind11"):
        try:
            distributions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            distributions[name] = None
    packages = {}
    for file in (Path(sys.prefix) / "conda-meta").glob("*.json"):
        package = json.loads(file.read_text())
        if package["name"] in {
            "eigen",
            "libboost-headers",
            "cspice",
            "tudat-resources",
            "cereal",
            "sofa-cmake",
            "nrlmsise-00",
            "pybind11",
            "numpy",
        }:
            packages[package["name"]] = {key: package.get(key) for key in ("version", "build")}
    cpu = platform.processor()
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        cpu = next(
            (
                line.split(":", 1)[1].strip()
                for line in cpuinfo.read_text().splitlines()
                if line.startswith("model name")
            ),
            cpu,
        )
    gravity_file = Path(resource_paths.get_gravity_models_path()) / "Mars/jgmro120d.txt"
    return {
        "kind": "metadata",
        "schema_version": 1,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "command_line": [sys.executable, *sys.argv],
        "configuration": asdict(config),
        "input": {
            "first_epoch": 778118400.0,
            "arc_gap_seconds": 3600.0,
            "starts": starts.tolist(),
            "initial_states_m_m_per_s": states.tolist(),
            "reference_radius_m": model.reference_radius,
            "mu_m3_per_s2": model.gravitational_parameter,
            "keplerian_recipe": "a=R+300000, e=.001, i=93deg, omega=.2, Omega=.3, nu=fmod(.4+.7*arc,2pi)",
            "gravity_file": str(gravity_file),
            "gravity_file_sha256": sha256(gravity_file),
        },
        "code": {
            "head": command_output(["git", "rev-parse", "HEAD"]),
            "branch": command_output(["git", "branch", "--show-current"]),
            "status_porcelain": command_output(["git", "status", "--porcelain=v1"]),
            "driver_sha256": sha256(Path(__file__)),
            "tracked_diff_sha256": hashlib.sha256(
                (command_output(["git", "diff", "HEAD"]) or "").encode()
            ).hexdigest(),
        },
        "environment": {
            "python": sys.version,
            "executable": sys.executable,
            "tudatpy": tudatpy.__version__,
            "kernel_path": str(module),
            "kernel_size_bytes": module.stat().st_size,
            "kernel_mtime_ns": module.stat().st_mtime_ns,
            "platform": platform.platform(),
            "cpu": cpu,
            "logical_cpus": os.cpu_count(),
            "cpu_affinity": (
                sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None
            ),
            "thread_environment": {
                name: os.environ.get(name)
                for name in (
                    "OPENBLAS_NUM_THREADS",
                    "OMP_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "BLIS_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS",
                    "NUMEXPR_NUM_THREADS",
                )
            },
            "python_distributions": distributions,
            "conda_packages": packages,
            "cmake_cache": cache,
            "compiler_version": (
                command_output([cache["CMAKE_CXX_COMPILER"], "--version"])
                if cache.get("CMAKE_CXX_COMPILER")
                else None
            ),
        },
        "timing_contract": {
            "setup": "fresh native model descriptor and facade; input generation excluded",
            "integration": "integrate_full_equations, including serial worker factory and native postprocessing",
            "total": "setup + integration",
            "excluded": [
                "imports",
                "SPICE kernel loading",
                "warmup",
                "sequential reference",
                "validation/result extraction",
                "destruction",
                "JSON output",
            ],
        },
        "historical_native_comparison": {
            "directly_comparable": False,
            "matching_recipe_source": "benchmarks/cpp/marsArcBenchmark.h:makeSettings",
            "differences": [
                "Python binding uses tudat::Time; historical native benchmark uses double time",
                "native benchmark also saves relative-position dependent variables and updates their interpolator",
                "historical timed native one-worker run uses isolated factory; Python one-worker run uses true sequential default",
                "Python setup also constructs a model descriptor; initialization/timing boundaries differ",
            ],
        },
    }


def make_inputs(config, model):
    from tudatpy.astro import element_conversion

    starts = 778118400.0 + np.arange(config.arcs) * (config.duration + 3600.0)
    states = np.array(
        [
            element_conversion.keplerian_to_cartesian(
                np.array(
                    [
                        model.reference_radius + 300000.0,
                        0.001,
                        93.0 * math.pi / 180.0,
                        0.2,
                        0.3,
                        math.fmod(0.4 + 0.7 * arc, 2 * math.pi),
                    ]
                ),
                model.gravitational_parameter,
            )
            for arc in range(config.arcs)
        ]
    )
    return starts, states


def integrate(config, starts, states, workers):
    from tudatpy.dynamics import simulator

    start = perf_counter()
    model = simulator.MarsGravityArcModel(config.degree, config.order)
    options = {} if workers == 1 else {"number_of_workers": workers}
    simulation = simulator.IndependentArcVariationalSimulator(
        model=model,
        arc_start_times=starts,
        initial_states=states,
        arc_duration=config.duration,
        step_size=config.step_size,
        **options,
    )
    setup_end = perf_counter()
    simulation.integrate_full_equations()
    integration_end = perf_counter()
    used = simulation.last_integration_used_independent_arcs
    if used != (workers > 1):
        raise RuntimeError(f"workers={workers}: unexpected independent route {used}")
    return simulation, {
        "setup_seconds": setup_end - start,
        "integration_seconds": integration_end - setup_end,
        "total_seconds": integration_end - start,
        "used_independent_arcs": used,
    }


def checked_difference(actual, expected, label):
    actual, expected = np.asarray(actual), np.asarray(expected)
    if (
        actual.shape != expected.shape
        or not np.isfinite(actual).all()
        or not np.isfinite(expected).all()
    ):
        raise ValueError(f"{label}: shape mismatch or non-finite values")
    absolute = np.abs(actual - expected)
    scaled = absolute / (ATOL + RTOL * np.abs(expected))
    maximum = float(np.max(scaled))
    if maximum > 1:
        raise ValueError(f"{label}: maximum tolerance ratio {maximum:.6g} exceeds 1")
    return {"max_absolute_difference": float(np.max(absolute)), "max_tolerance_ratio": maximum}


def snapshot(simulation, config, starts):
    dynamics = simulation.propagation_results.single_arc_results
    variations = simulation.variational_propagation_results.single_arc_results
    if len(dynamics) != config.arcs or len(variations) != config.arcs:
        raise ValueError("unexpected native arc count")
    arcs = []
    for arc, (state_result, variation_result) in enumerate(zip(dynamics, variations)):
        histories = [
            state_result.state_history,
            variation_result.state_transition_matrix_history,
            variation_result.sensitivity_matrix_history,
        ]
        times = np.array(list(histories[0]))
        arrays = [np.array(list(history.values())) for history in histories]
        if len(times) < 8 or not np.all(np.diff(times) > 0):
            raise ValueError("insufficient or unordered epochs")
        for history, values, shape in zip(histories, arrays, [(6,), (6, 6), (6, 1)]):
            if (
                not np.array_equal(times, np.array(list(history)))
                or values.shape != (len(times), *shape)
                or not np.isfinite(values).all()
            ):
                raise ValueError("incompatible or non-finite native histories")
        if times[0] != starts[arc] or times[-1] != starts[arc] + config.duration:
            raise ValueError("unexpected integration interval")
        interface = simulation.state_transition_interface
        probes = []
        # At stored epochs, check scatter mapping against the local histories.
        for index in sorted({0, len(times) // 2, len(times) - 1}):
            full = interface.full_state_transition_sensitivity_at_epoch(float(times[index]))
            expected = np.zeros((6, 6 * config.arcs + 1))
            expected[:, 6 * arc : 6 * (arc + 1)] = arrays[1][index]
            expected[:, -1:] = arrays[2][index]
            checked_difference(full, expected, "global STM/shared-mu mapping")
            probes.append((float(times[index]), full))
        # Also compare interpolation between saved epochs against the reference.
        epoch = float((times[1] + times[2]) / 2)
        full = interface.full_state_transition_sensitivity_at_epoch(epoch)
        if full.shape != (6, 6 * config.arcs + 1) or not np.isfinite(full).all():
            raise ValueError("invalid interpolated global matrix")
        probes.append((epoch, full))
        if not np.any(arrays[2][len(times) // 2]):
            raise ValueError("shared-mu sensitivity is identically zero at midpoint")
        arcs.append(
            {
                "times": times,
                "states": arrays[0],
                "stms": arrays[1],
                "sensitivities": arrays[2],
                "global_probes": probes,
            }
        )
    return arcs


def compare_snapshots(actual, expected):
    if len(actual) != len(expected):
        raise ValueError("arc count mismatch")
    metrics = {
        name: {"max_absolute_difference": 0.0, "max_tolerance_ratio": 0.0}
        for name in ("states", "stms", "sensitivities", "global_matrices")
    }
    epochs = 0
    probes = 0
    for arc, (first, second) in enumerate(zip(actual, expected)):
        if not np.array_equal(first["times"], second["times"]):
            raise ValueError(f"arc {arc}: epoch mismatch")
        epochs += len(first["times"])
        for name in ("states", "stms", "sensitivities"):
            result = checked_difference(first[name], second[name], f"arc {arc} {name}")
            for key, value in result.items():
                metrics[name][key] = max(metrics[name][key], value)
        if len(first["global_probes"]) != len(second["global_probes"]):
            raise ValueError("global probe count mismatch")
        for (time, matrix), (other_time, other_matrix) in zip(
            first["global_probes"], second["global_probes"]
        ):
            if time != other_time:
                raise ValueError("global probe epoch mismatch")
            result = checked_difference(matrix, other_matrix, f"arc {arc} global matrix")
            for key, value in result.items():
                metrics["global_matrices"][key] = max(metrics["global_matrices"][key], value)
            probes += 1
    return {
        "status": "passed",
        "checked_epochs": epochs,
        "global_matrix_probes": probes,
        "atol": ATOL,
        "rtol": RTOL,
        "differences": metrics,
    }


def worker_schedule(workers, repeats):
    return [
        workers[index % len(workers) :] + workers[: index % len(workers)]
        for index in range(repeats)
    ]


def summarize(records, workers):
    result = []
    for count in workers:
        selected = [row for row in records if row["workers"] == count]
        summary = {"workers": count, "repetitions": len(selected)}
        for name in ("setup_seconds", "integration_seconds", "total_seconds"):
            values = [row[name] for row in selected]
            summary[name] = {
                "median": statistics.median(values),
                "min": min(values),
                "max": max(values),
            }
        result.append(summary)
    baseline = next(row for row in result if row["workers"] == 1)
    for row in result:
        row["total_speedup"] = baseline["total_seconds"]["median"] / row["total_seconds"]["median"]
        row["integration_speedup"] = (
            baseline["integration_seconds"]["median"] / row["integration_seconds"]["median"]
        )
    return {"kind": "summary", "workers": result, "validation": "passed"}


def execute(config, emit):
    from tudatpy.data_input.environment_data import spice
    from tudatpy.dynamics import simulator

    spice.load_standard_kernels()
    model = simulator.MarsGravityArcModel(config.degree, config.order)
    starts, states = make_inputs(config, model)
    emit(collect_metadata(config, starts, states, model))
    # Same short warmup on every invocation; it is never included in timings.
    warmup = Configuration("warmup", 1, 2, 2, 600 / 86400, 30.0, (1,), 1)
    warm_starts, warm_states = make_inputs(warmup, model)
    warm_simulation, _ = integrate(warmup, warm_starts, warm_states, 1)
    del warm_simulation
    print("Building untimed sequential reference...", file=sys.stderr, flush=True)
    reference_simulation, _ = integrate(config, starts, states, 1)
    reference = snapshot(reference_simulation, config, starts)
    emit({"kind": "reference", "workers": 1, "validation": compare_snapshots(reference, reference)})
    del reference_simulation
    records = []
    for repetition, order in enumerate(worker_schedule(config.workers, config.repeats), start=1):
        for position, workers in enumerate(order, start=1):
            print(
                f"Repetition {repetition}/{config.repeats}, workers={workers}",
                file=sys.stderr,
                flush=True,
            )
            simulation, timing = integrate(config, starts, states, workers)
            validation_start = perf_counter()
            actual = snapshot(simulation, config, starts)
            validation = compare_snapshots(actual, reference)
            record = {
                "kind": "run",
                "repetition": repetition,
                "position": position,
                "workers": workers,
                **timing,
                "validation": validation,
                "validation_seconds": perf_counter() - validation_start,
            }
            emit(record)
            records.append(record)
            del actual, simulation
    summary = summarize(records, config.workers)
    summary["completed_utc"] = datetime.now(timezone.utc).isoformat()
    emit(summary)
    for row in summary["workers"]:
        print(
            f'workers={row["workers"]}: total median {row["total_seconds"]["median"]:.6f}s '
            f'(range {row["total_seconds"]["min"]:.6f}..{row["total_seconds"]["max"]:.6f}), '
            f'speedup {row["total_speedup"]:.2f}x',
            file=sys.stderr,
        )


def main(argv=None):
    config, output = parse_arguments(argv)
    try:
        context = output.open("x") if output else nullcontext(sys.stdout)
        with context as stream:

            def emit(record):
                print(json.dumps(record, allow_nan=False), file=stream, flush=True)

            try:
                execute(config, emit)
            except Exception as error:
                emit({"kind": "error", "exception": type(error).__name__, "message": str(error)})
                raise
    except (OSError, ValueError, RuntimeError) as error:
        print(f"Benchmark failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
