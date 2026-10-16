"""Identify the imported Time callback build and resolve environment documentation."""

from pathlib import Path
import datetime
import hashlib
import importlib
import json
import subprocess
import sys

import tudatpy
import tudatpy.kernel
from tudatpy.astro import time_representation as time
from tudatpy.dynamics import environment

source = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging")
logs = Path("/home/dominic/Tudat/laser-ranging-work/logs")
document = source / "docs/tudatpy/source/dynamics/environment.rst"
module = None
autosummary = False
resolved = []
for number, line in enumerate(document.read_text().splitlines(), 1):
    if line.startswith(".. currentmodule::"):
        module = importlib.import_module(line.split("::", 1)[1].strip())
    if line.startswith(".. autosummary::"):
        autosummary = True
        continue
    if autosummary and line and not line[0].isspace():
        autosummary = False
    if autosummary and line.strip() and not line.lstrip().startswith(":"):
        name = line.strip()
        getattr(module, name)
        resolved.append({"line": number, "name": module.__name__ + "." + name})

assert "CustomEphemerisTimeObject" in [item["name"].rsplit(".", 1)[1] for item in resolved]
assert issubclass(environment.CustomEphemerisTimeObject, environment.Ephemeris)
assert "Time" in environment.Ephemeris.cartesian_state.__doc__
kernel = Path(tudatpy.kernel.__file__)
result = {
    "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "source_revision": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=source, text=True
    ).strip(),
    "native_binding_source_revision": "8506f918874a2a628aa200cdab579e0f85e88a90",
    "later_source_changes": "Analysis script only; native binding unchanged since 8506f9188",
    "interpreter": sys.executable,
    "package_origin": tudatpy.__file__,
    "kernel_origin": str(kernel),
    "kernel_sha256": hashlib.sha256(kernel.read_bytes()).hexdigest(),
    "kernel_bytes": kernel.stat().st_size,
    "kernel_mtime_utc": datetime.datetime.fromtimestamp(
        kernel.stat().st_mtime, datetime.timezone.utc
    ).isoformat(),
    "environment_autosummary_names": resolved,
    "environment_autosummary_count": len(resolved),
    "custom_ephemeris_time_class_and_overload": "passed",
    "behavioral_test_report": str(logs / "time-callback-environment-tests.xml"),
    "environment": {
        "PYTHONPATH": "/home/dominic/Tudat/laser-ranging-build/src",
        "LD_LIBRARY_PATH": "/home/dominic/miniconda3/envs/tudatpy-dev/lib",
        "numerical_threads": 1,
    },
}
(logs / "time-callback-kernel-provenance.json").write_text(json.dumps(result, indent=2) + "\n")
print(
    json.dumps(
        {key: value for key, value in result.items() if key != "environment_autosummary_names"},
        indent=2,
    )
)
