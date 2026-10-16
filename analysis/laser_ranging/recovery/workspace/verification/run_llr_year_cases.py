"""Record the authorized three annual cases using separate native interpreters."""

from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
import time

WORK = Path("/home/dominic/Tudat/laser-ranging-work")
SOURCE = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-signature")
PYTHON = Path("/home/dominic/miniconda3/envs/tudatpy-dev/bin/python")
KERNEL = Path("/home/dominic/Tudat/laser-ranging-build/src/tudatpy/kernel.so")
EXPECTED = "cab49e738cfcf0ea5f251cb140ab8919d5372dbc3e5ba09314d0d5be4eedfcef"
PHASE = sys.argv[1]
assert PHASE in ("boundary", "full")
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
assert sha(KERNEL) == EXPECTED
manifest = SOURCE / "analysis/laser_ranging/llr_year_2017_dataset.json"
assert sha(manifest) == "4140847ee51866fed8d3d2a92d6ead4d03113a7e67885adcbfc52958d9174446"
env = os.environ.copy()
env.update(
    PYTHONPATH="/home/dominic/Tudat/laser-ranging-build/src",
    LD_LIBRARY_PATH="/home/dominic/miniconda3/envs/tudatpy-dev/lib",
)
for name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
):
    env[name] = "1"
INTEGRATION = "9082ce261"
for relative in (
    "llr_residuals.py",
    "llr_long_arc.py",
    "llr_year_2017_dataset.json",
    "inpop17a_dataset.json",
    "plot_llr_residuals.py",
):
    archived = subprocess.check_output(
        ["git", "show", f"{INTEGRATION}:analysis/laser_ranging/{relative}"], cwd=SOURCE
    )
    assert hashlib.sha256(archived).hexdigest() == sha(SOURCE / "analysis/laser_ranging" / relative)
revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=SOURCE, text=True).strip()
assert not subprocess.check_output(["git", "status", "--porcelain"], cwd=SOURCE, text=True)
cases = [
    ("de440", "de440", "jpl-nh"),
    ("inpop17a_jpl_nh", "inpop17a", "jpl-nh"),
    ("inpop17a_no_separate_glass", "inpop17a", "none"),
]
if len(sys.argv) > 2:
    cases = [case for case in cases if case[0] == sys.argv[2]]
    assert len(cases) == 1
pending = []
records = []


def start_case(name, ephemeris, glass):
    """Launch a single process with immutable inputs and a fresh output directory."""
    output = WORK / "results/llr-year-2017-attempt2" / PHASE / name
    assert not output.exists()
    command = [
        str(PYTHON),
        str(SOURCE / "analysis/laser_ranging/llr_long_arc.py"),
        "--dataset-manifest",
        str(manifest),
        "--data-root",
        str(WORK / "research/llr-year-2017"),
        "--kernel-dir",
        str(WORK / ("data/llr/kernels" if ephemeris == "de440" else "data/llr/inpop17a-kernels")),
        "--ephemeris",
        ephemeris,
        "--glass-convention",
        glass,
        "--output",
        str(output),
    ]
    if PHASE == "boundary":
        command.append("--first-last-control")
    log = WORK / "logs" / f"llr-year-2017-attempt2-{PHASE}-{name}.log"
    stream = log.open("w")
    began = time.perf_counter()
    process = subprocess.Popen(
        command, cwd=SOURCE, env=env, stdout=stream, stderr=subprocess.STDOUT
    )
    record = dict(
        case=name,
        phase=PHASE,
        command=command,
        source_revision=revision,
        integration_source_revision=INTEGRATION,
        started_utc=datetime.now(timezone.utc).isoformat(),
        kernel_sha256=EXPECTED,
        manifest_sha256=sha(manifest),
        source_sha256={
            path.name: sha(path)
            for path in (
                SOURCE / "analysis/laser_ranging/llr_long_arc.py",
                SOURCE / "analysis/laser_ranging/llr_residuals.py",
            )
        },
        numerical_threads=1,
        process_id=process.pid,
        output=str(output),
        log=str(log),
    )
    return process, stream, began, record


def finish_case(item):
    """Record completed timings and validate the declared retained/source counts."""
    process, stream, began, record = item
    status = process.wait()
    stream.close()
    record.update(
        returncode=status,
        external_wall_seconds=time.perf_counter() - began,
        ended_utc=datetime.now(timezone.utc).isoformat(),
    )
    failed_path = WORK / "logs" / f"llr-year-2017-attempt2-{PHASE}-{record['case']}-provenance.json"
    failed_path.write_text(json.dumps(record, indent=2) + "\n")
    assert status == 0, (record["case"], status)
    output = Path(record["output"])
    summary = json.loads((output / "summary.json").read_text())
    assert summary["kernel_binary_sha256"] == EXPECTED
    assert summary["count"] == (4 if PHASE == "boundary" else 1579)
    if PHASE == "full":
        assert (
            summary["source_count"] == 1582
            and summary["excluded_count"] == 3
            and summary["control_unselected_count"] == 0
        )
        assert summary["station_statistics"]["APOLLO"]["count"] == 311
        assert summary["station_statistics"]["GRASSE"]["count"] == 1268
    rows = list(csv.DictReader((output / "residuals.csv").open()))
    assert len(rows) == summary["count"] and all(
        math.isfinite(float(row["one_way_residual_m"])) for row in rows
    )
    record.update(
        csv_sha256=sha(output / "residuals.csv"),
        summary_sha256=sha(output / "summary.json"),
        count=summary["count"],
        setup_seconds=summary["setup_seconds"],
        solver_seconds=summary["solver_seconds"],
        total_seconds=summary["total_seconds"],
    )
    path = WORK / "logs" / f"llr-year-2017-attempt2-{PHASE}-{record['case']}-provenance.json"
    path.write_text(json.dumps(record, indent=2) + "\n")
    records.append(record)
    print(json.dumps(record), flush=True)


for case in cases:
    item = start_case(*case)
    if PHASE == "boundary":
        finish_case(item)
    else:
        pending.append(item)
for item in pending:
    finish_case(item)
assert sha(KERNEL) == EXPECTED
(
    WORK
    / "logs"
    / f'llr-year-2017-attempt2-{PHASE}-{sys.argv[2] if len(sys.argv)>2 else "all"}-summary.json'
).write_text(json.dumps(records, indent=2) + "\n")
