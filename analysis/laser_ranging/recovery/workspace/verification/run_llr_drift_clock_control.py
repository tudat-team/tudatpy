"""Run one authorized all-306 drift-clock control with pinned kernel and serial numerical threads."""

from pathlib import Path
import datetime, hashlib, json, os, subprocess, sys, time

SOURCE = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-signature")
WORK = Path("/home/dominic/Tudat/laser-ranging-work")
label = sys.argv[1]
assert label in ["baseline", "candidate", "baseline-repeat", "candidate-repeat"]
out = WORK / "results" / ("llr-drift-clock-" + label + "-306")
assert not out.exists()
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
kernel = Path("/home/dominic/Tudat/laser-ranging-build/src/tudatpy/kernel.so")
assert sha(kernel) == "cab49e738cfcf0ea5f251cb140ab8919d5372dbc3e5ba09314d0d5be4eedfcef"
env = os.environ.copy()
env["PYTHONPATH"] = "/home/dominic/Tudat/laser-ranging-build/src"
env["LD_LIBRARY_PATH"] = "/home/dominic/miniconda3/envs/tudatpy-dev/lib"
for name in [
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
]:
    env[name] = "1"
command = [
    "/home/dominic/miniconda3/envs/tudatpy-dev/bin/python",
    str(SOURCE / "analysis/laser_ranging/llr_long_arc.py"),
    "--data-root",
    str(WORK / "research/longarc-20170304"),
    "--kernel-dir",
    str(WORK / "data/llr/kernels"),
    "--output",
    str(out),
]
record = {
    "label": label,
    "command": command,
    "source_revision": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE, text=True
    ).strip(),
    "source_dirty": subprocess.check_output(["git", "status", "--short"], cwd=SOURCE, text=True),
    "source_sha256": {
        p.name: sha(p)
        for p in (
            SOURCE / "analysis/laser_ranging/llr_long_arc.py",
            SOURCE / "analysis/laser_ranging/llr_residuals.py",
        )
    },
    "kernel_sha256": sha(kernel),
    "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "numerical_threads": 1,
}
started = time.perf_counter()
with (WORK / "logs" / ("llr-drift-clock-" + label + "-306.log")).open("w") as log:
    r = subprocess.run(command, cwd=SOURCE, env=env, stdout=log, stderr=subprocess.STDOUT)
record.update(
    returncode=r.returncode,
    external_wall_seconds=time.perf_counter() - started,
    ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
)
assert record["kernel_sha256"] == sha(kernel)
if r.returncode == 0:
    summary = json.loads((out / "summary.json").read_text())
    record.update(
        {k: summary[k] for k in ["count", "setup_seconds", "solver_seconds", "total_seconds"]}
    )
    record["csv_sha256"] = sha(out / "residuals.csv")
(WORK / "logs" / ("llr-drift-clock-" + label + "-306-provenance.json")).write_text(
    json.dumps(record, indent=2) + "\n"
)
print(json.dumps(record, indent=2))
raise SystemExit(r.returncode)
