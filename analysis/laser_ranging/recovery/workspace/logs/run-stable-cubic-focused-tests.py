"""Run only the cubic interpolation and relevant n-way range suites, serially."""

from pathlib import Path
import datetime
import hashlib
import json
import os
import subprocess
import xml.etree.ElementTree as ET

BUILD = Path("/home/dominic/Tudat/laser-ranging-build")
LOGS = Path("/home/dominic/Tudat/laser-ranging-work/logs")
env = os.environ.copy()
env["LD_LIBRARY_PATH"] = "/home/dominic/miniconda3/envs/tudatpy-dev/lib"
for key in [
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
]:
    env[key] = "1"
results = []
for target, relative in [
    ("test_interpolators_CubicSplineInterpolator", "math/interpolators"),
    ("test_observation_models_NWayRangeObservationModel", "astro/observation_models"),
]:
    exe = BUILD / "tests" / target
    report = LOGS / (target + "-stable-cubic.xml")
    log = LOGS / (target + "-stable-cubic.log")
    command = [
        str(exe),
        "--report_level=detailed",
        "--report_format=XML",
        "--report_sink=" + str(report),
        "--log_level=test_suite",
    ]
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with log.open("w") as out:
        completed = subprocess.run(
            command,
            cwd=BUILD / "tests/test_tudat/src" / relative,
            env=env,
            stdout=out,
            stderr=subprocess.STDOUT,
        )
    root = ET.parse(report).getroot()
    suite = root.find("TestSuite") if root.tag == "TestResult" else root
    result = {
        "target": target,
        "returncode": completed.returncode,
        "native_source_revision": "b1aa772503070347e8ac0b263c8de9278d390c17",
        "executable_sha256": hashlib.sha256(exe.read_bytes()).hexdigest(),
        "command": command,
        "report": str(report),
        "log": str(log),
        "started_utc": started,
        "completed_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "boost_report": suite.attrib,
    }
    results.append(result)
    (LOGS / "stable-cubic-focused-native-results.json").write_text(
        json.dumps({"results": results}, indent=2) + "\n"
    )
    print(json.dumps(result), flush=True)
    if completed.returncode:
        break
raise SystemExit(any(result["returncode"] for result in results))
