"""Record the authorized focused n-way native validation without broadening suites."""

from pathlib import Path
import datetime, hashlib, json, os, subprocess, xml.etree.ElementTree as ET

build = Path("/home/dominic/Tudat/laser-ranging-build")
source = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging")
logs = Path("/home/dominic/Tudat/laser-ranging-work/logs")
target = "test_observation_models_NWayRangeObservationModel"
exe = build / "tests" / target
report = logs / (target + "-selected-clock-framefix.xml")
log = logs / (target + "-selected-clock-framefix.log")
env = os.environ.copy()
env["LD_LIBRARY_PATH"] = "/home/dominic/miniconda3/envs/tudatpy-dev/lib"
for key in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
]:
    env[key] = "1"
command = [
    str(exe),
    "--report_level=detailed",
    "--report_format=XML",
    "--report_sink=" + str(report),
    "--log_level=test_suite",
]
record = {
    "target": target,
    "native_source_revision": "7d52716f6",
    "integration_revision_at_test": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=source, text=True
    ).strip(),
    "executable_sha256": hashlib.sha256(exe.read_bytes()).hexdigest(),
    "command": command,
    "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "log": str(log),
    "report": str(report),
}
with log.open("w") as output:
    completed = subprocess.run(
        command,
        cwd=build / "tests/test_tudat/src/astro/observation_models",
        env=env,
        stdout=output,
        stderr=subprocess.STDOUT,
    )
record.update(
    returncode=completed.returncode,
    ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
)
if report.exists():
    record["boost_report"] = ET.parse(report).getroot().find("TestSuite").attrib
(logs / "selected-clock-nway-native-framefix-result.json").write_text(
    json.dumps(record, indent=2) + "\n"
)
print(json.dumps(record, indent=2))
raise SystemExit(completed.returncode)
