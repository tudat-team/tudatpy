"""Benchmark reporting and fail-loud validation, including a real short CLI run."""

import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

DRIVER = Path(__file__).resolve().parents[2] / "benchmarks/benchmark_python_independent_arcs.py"
SPEC = importlib.util.spec_from_file_location("python_arc_benchmark", DRIVER)
benchmark = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = benchmark
SPEC.loader.exec_module(benchmark)


@pytest.mark.parametrize(
    "arguments",
    [
        ["--arcs", "0"],
        ["--arcs", "4"],
        ["--workers", "2", "4"],
        ["--workers", "1", "2", "2"],
        ["--workers", "1", "0"],
        ["--repeats", "0"],
        ["--degree", "121"],
        ["--degree", "2", "--order", "3"],
        ["--days-per-arc", "nan"],
        ["--step-size", "inf"],
        ["--step-size", "100"],
        ["--days-per-arc", "-1"],
    ],
)
def test_invalid_workloads_are_rejected(arguments):
    with pytest.raises(SystemExit) as error:
        benchmark.parse_arguments(arguments)
    assert error.value.code == 2


def test_presets_and_rotating_schedule():
    smoke, _ = benchmark.parse_arguments([])
    heavy, _ = benchmark.parse_arguments(["--preset", "representative"])
    assert smoke.arcs == 8 and smoke.duration == 600 and smoke.repeats == 1
    assert heavy.arcs == 8 and heavy.degree == heavy.order == 50
    assert heavy.duration == 172800 and heavy.repeats == 3
    assert benchmark.worker_schedule(smoke.workers, 3) == [
        (1, 2, 4, 8),
        (2, 4, 8, 1),
        (4, 8, 1, 2),
    ]


@pytest.mark.parametrize(
    "actual,reference",
    [
        ([np.nan], [np.nan]),
        ([np.inf], [1.0]),
        ([0.0], [np.inf]),
        ([[1.0]], [1.0]),
        ([2e-14], [0.0]),
        ([1.000000001], [1.0]),
    ],
)
def test_nonfinite_mismatched_or_out_of_tolerance_data_fails(actual, reference):
    with pytest.raises(ValueError):
        benchmark.checked_difference(actual, reference, "injected corruption")


def test_validation_metrics_do_not_hide_near_zero_errors():
    metrics = benchmark.checked_difference([5e-15, 1.0], [0.0, 1.0], "small values")
    assert metrics["max_absolute_difference"] == pytest.approx(5e-15)
    assert metrics["max_tolerance_ratio"] == pytest.approx(0.5)


def synthetic_snapshot():
    return [
        {
            "times": np.arange(8.0),
            "states": np.ones((8, 6)),
            "stms": np.ones((8, 6, 6)),
            "sensitivities": np.ones((8, 6, 1)),
            "global_probes": [(3.0, np.ones((6, 7)))],
        }
    ]


@pytest.mark.parametrize("name", ["states", "stms", "sensitivities", "global_matrices", "times"])
def test_corruption_in_every_compared_quantity_is_detected(name):
    expected, actual = synthetic_snapshot(), synthetic_snapshot()
    if name == "global_matrices":
        actual[0]["global_probes"][0][1][0, -1] += 1.0
    else:
        actual[0][name].flat[0] += 1.0
    with pytest.raises(ValueError):
        benchmark.compare_snapshots(actual, expected)


def test_summary_uses_median_and_retains_spread():
    records = []
    for workers, totals in [(1, [4.0, 5.0, 6.0]), (2, [1.0, 2.0, 8.0])]:
        for total in totals:
            records.append(
                dict(
                    workers=workers,
                    setup_seconds=0.1,
                    integration_seconds=total - 0.1,
                    total_seconds=total,
                )
            )
    summary = benchmark.summarize(records, (1, 2))
    parallel = summary["workers"][1]
    assert parallel["total_seconds"] == dict(median=2.0, min=1.0, max=8.0)
    assert parallel["total_speedup"] == pytest.approx(2.5)
    # An outlier remains visible rather than being silently discarded.
    assert parallel["repetitions"] == 3


def test_existing_output_is_never_overwritten(tmp_path):
    output = tmp_path / "existing.jsonl"
    output.write_text("preserve this file\n")
    assert benchmark.main(["--output", str(output)]) == 1
    assert output.read_text() == "preserve this file\n"


def test_real_cli_sequential_parallel_and_metadata(tmp_path):
    import tudatpy.kernel

    output = tmp_path / "smoke.jsonl"
    subprocess.run(
        [sys.executable, str(DRIVER), "--degree", "2", "--repeats", "2", "--output", str(output)],
        check=True,
        capture_output=True,
        text=True,
        timeout=90,
    )
    rows = [json.loads(line) for line in output.read_text().splitlines()]
    assert [row["kind"] for row in rows] == ["metadata", "reference"] + ["run"] * 8 + ["summary"]
    meta = rows[0]
    assert meta["environment"]["kernel_path"] == str(Path(tudatpy.kernel.__file__).resolve())
    assert meta["configuration"]["workers"] == [1, 2, 4, 8]
    assert meta["configuration"]["arcs"] == 8
    assert len(meta["input"]["gravity_file_sha256"]) == 64
    assert len(meta["code"]["driver_sha256"]) == 64
    assert meta["historical_native_comparison"]["directly_comparable"] is False
    runs = rows[2:-1]
    assert [row["workers"] for row in runs] == [1, 2, 4, 8, 2, 4, 8, 1]
    for row in runs:
        assert row["used_independent_arcs"] == (row["workers"] > 1)
        assert row["validation"]["status"] == "passed"
        assert row["validation"]["checked_epochs"] == 168
        assert row["validation"]["global_matrix_probes"] == 32
        assert row["total_seconds"] == pytest.approx(
            row["setup_seconds"] + row["integration_seconds"]
        )
        assert all(
            metric["max_tolerance_ratio"] <= 1
            for metric in row["validation"]["differences"].values()
        )
    # Initial-state recipe matches the native benchmark's Keplerian semimajor axis.
    state = np.array(meta["input"]["initial_states_m_m_per_s"][0])
    mu = meta["input"]["mu_m3_per_s2"]
    energy = 0.5 * np.dot(state[3:], state[3:]) - mu / np.linalg.norm(state[:3])
    assert energy == pytest.approx(
        -mu / (2 * (meta["input"]["reference_radius_m"] + 300000)), rel=1e-12
    )
    assert rows[-1]["validation"] == "passed"
