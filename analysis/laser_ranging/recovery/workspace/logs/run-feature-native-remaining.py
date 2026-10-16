"""Run only the five selected laser-ranging native suites, serially."""

from pathlib import Path
import datetime
import hashlib
import json
import os
import subprocess
import xml.etree.ElementTree as ET

build = Path("/home/dominic/Tudat/laser-ranging-build")
source = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging")
logs = Path("/home/dominic/Tudat/laser-ranging-work/logs")
targets = [
    ("io", "CrdAndSinexReaders"),
    ("observation_models", "NWayRangeObservationModel"),
    ("observation_models", "AtmosphereCorrection"),
]
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
for group, name in targets:
    target = "test_" + group + "_" + name
    executable = build / "tests" / target
    working_directory = (
        build / "tests/test_tudat/src" / ("io" if group == "io" else "astro/observation_models")
    )
    report = logs / (target + ".xml")
    command = [
        str(executable),
        "--report_level=detailed",
        "--report_format=XML",
        "--report_sink=" + str(report),
        "--log_level=test_suite",
    ]
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with (logs / (target + ".log")).open("w") as output:
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
        "log": str(logs / (target + ".log")),
    }
    if report.exists():
        root = ET.parse(report).getroot()
        result["boost_report"] = root.attrib
        result["test_cases"] = [case.attrib for case in root.iter("TestCase")]
    results.append(result)
    (logs / "feature-native-results-remaining.json").write_text(
        json.dumps({"source_revision": revision, "results": results}, indent=2) + "\n"
    )
    print(
        target, "exit", completed.returncode, "cases", len(result.get("test_cases", [])), flush=True
    )
    if completed.returncode:
        break
raise SystemExit(any(result["returncode"] for result in results))
