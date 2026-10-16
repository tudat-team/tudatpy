"""Record the authorized existing environment tests and the new kernel provenance."""

from pathlib import Path
import datetime, hashlib, json, os, subprocess, sys
import tudatpy
import tudatpy.kernel

work = Path("/home/dominic/Tudat/laser-ranging-work")
source = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging")
logs = work / "logs"
kernel = Path(tudatpy.kernel.__file__)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
record = {
    "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "native_source_revision": subprocess.check_output(
        ["git", "rev-parse", "7d52716f6"], cwd=source, text=True
    ).strip(),
    "integration_revision_at_import": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=source, text=True
    ).strip(),
    "kernel_origin": str(kernel),
    "package_origin": tudatpy.__file__,
    "interpreter": sys.executable,
    "kernel_sha256": sha(kernel),
    "kernel_bytes": kernel.stat().st_size,
    "kernel_mtime_utc": datetime.datetime.fromtimestamp(
        kernel.stat().st_mtime, datetime.timezone.utc
    ).isoformat(),
    "previous_kernel_sha256": "10f8307ba799847f11d9614d2c3f3da3db23962c46a739d73f536f9dc31ae5d5",
    "compiler_job_limit": 8,
    "build_log": str(logs / "selected-clock-time-rotation-build.log"),
    "build_command": [
        "cmake",
        "--build",
        "/home/dominic/Tudat/laser-ranging-build",
        "--target",
        "test_observation_models_NWayRangeObservationModel",
        "kernel",
        "--parallel",
        "8",
    ],
    "environment": {
        "PYTHONPATH": os.environ.get("PYTHONPATH"),
        "LD_LIBRARY_PATH": os.environ.get("LD_LIBRARY_PATH"),
        "numerical_threads": 1,
    },
}
assert record["kernel_sha256"] != record["previous_kernel_sha256"]
command = [
    sys.executable,
    "-m",
    "pytest",
    str(source / "tests/test_tudatpy/test_environment_setup.py"),
    "-q",
    "--junitxml=" + str(logs / "selected-clock-environment-python.xml"),
]
with (logs / "selected-clock-environment-python.log").open("w") as output:
    completed = subprocess.run(
        command, cwd=source, env=os.environ.copy(), stdout=output, stderr=subprocess.STDOUT
    )
record["python_test_command"] = command
record["python_test_returncode"] = completed.returncode
assert record["kernel_sha256"] == sha(kernel)
(logs / "selected-clock-time-rotation-kernel-provenance.json").write_text(
    json.dumps(record, indent=2) + "\n"
)
print(json.dumps(record, indent=2))
raise SystemExit(completed.returncode)
