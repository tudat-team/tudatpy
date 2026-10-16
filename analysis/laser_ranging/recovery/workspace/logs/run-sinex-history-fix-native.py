"""Run the three affected SINEX reader suites serially; retain unrelated passes."""

from pathlib import Path
import datetime
import hashlib
import json
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

build = Path("/home/dominic/Tudat/laser-ranging-build")
source = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging")
logs = Path("/home/dominic/Tudat/laser-ranging-work/logs")
suffix = os.environ.get("TUDAT_READER_CHECK_SUFFIX", "sinex-fix")
targets = ["test_io_LaserRangingData", "test_io_CrdAndSinexReaders", "test_io_SinexStationHistory"]
if len(sys.argv) > 1:
    targets = sys.argv[1:]
env = os.environ.copy()
env["LD_LIBRARY_PATH"] = "/home/dominic/miniconda3/envs/tudatpy-dev/lib"
for name in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
]:
    env[name] = "1"
revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()
results = []
for target in targets:
    executable = build / "tests" / target
    working_directory = build / "tests/test_tudat/src/io"
    report = logs / (target + "-" + suffix + ".xml")
    command = [
        str(executable),
        "--report_level=detailed",
        "--report_format=XML",
        "--report_sink=" + str(report),
        "--log_level=test_suite",
    ]
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with (logs / (target + "-" + suffix + ".log")).open("w") as output:
        completed = subprocess.run(
            command, cwd=working_directory, env=env, stdout=output, stderr=subprocess.STDOUT
        )
    result = {
        "target": target,
        "source_revision": revision,
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "command": command,
        "working_directory": str(working_directory),
        "started": started,
        "ended": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "returncode": completed.returncode,
        "report": str(report),
        "log": str(logs / (target + "-" + suffix + ".log")),
    }
    if report.exists():
        root = ET.parse(report).getroot()
        result["boost_report"] = (
            root.find("TestSuite").attrib if root.tag == "TestResult" else root.attrib
        )
        result["test_cases"] = [case.attrib for case in root.iter("TestCase")]
    results.append(result)
    (logs / (suffix + "-native-results.json")).write_text(
        json.dumps({"source_revision": revision, "results": results}, indent=2) + "\n"
    )
    print(
        target,
        "exit",
        completed.returncode,
        "cases",
        len(result.get("test_cases", [])),
        "report",
        result.get("boost_report"),
        flush=True,
    )
raise SystemExit(any(result["returncode"] for result in results))
