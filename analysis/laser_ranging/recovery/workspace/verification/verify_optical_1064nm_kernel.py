"""Check actual Grasse wavelength settings and record the focused kernel provenance."""

from pathlib import Path
import datetime, hashlib, json, os, subprocess, sys
import tudatpy
import tudatpy.kernel
from tudatpy.estimation.observable_models_setup import light_time_corrections as ltc

WORK = Path("/home/dominic/Tudat/laser-ranging-work")
SOURCE = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging")
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
factory = ltc.mendes_pavlis_tropospheric_light_time_correction
assert "1.0645e-6" in factory.__doc__
wavelengths = set()
for entry in json.loads((WORK / "research/longarc-20170304/SOURCES.json").read_text()):
    path = Path(entry["file"])
    assert sha(path) == entry["sha256"]
    for line in path.read_text().splitlines():
        fields = line.split()
        if fields and fields[0].lower() == "c0":
            wavelengths.add(float(fields[2]))
assert wavelengths == {532.1, 1064.2}
checks = []
for nm in sorted(wavelengths):
    settings = factory(nm * 1e-9)
    assert isinstance(settings, ltc.LightTimeCorrectionSettings)
    checks.append(
        {
            "source_vacuum_wavelength_nm": nm,
            "supplied_vacuum_wavelength_m": nm * 1e-9,
            "settings_type": type(settings).__name__,
            "result": "accepted",
        }
    )
kernel = Path(tudatpy.kernel.__file__)
record = {
    "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "native_source_revision": subprocess.check_output(
        ["git", "rev-parse", "a366f436c"], cwd=SOURCE, text=True
    ).strip(),
    "integration_revision_at_import": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE, text=True
    ).strip(),
    "interpreter": sys.executable,
    "package_origin": tudatpy.__file__,
    "kernel_origin": str(kernel),
    "kernel_sha256": sha(kernel),
    "kernel_bytes": kernel.stat().st_size,
    "kernel_mtime_utc": datetime.datetime.fromtimestamp(
        kernel.stat().st_mtime, datetime.timezone.utc
    ).isoformat(),
    "previous_kernel_sha256": "80b924287707124a8efea66455b26233b5511b7856983d35c05d21f82734f215",
    "build_log": str(WORK / "logs/optical-1064_2nm-focused-build.log"),
    "compiler_job_limit": 8,
    "python_check": {
        "checks": checks,
        "docstring_acceptance_convention": "passed",
        "source_crd_hashes_verified": 10,
    },
    "environment": {
        "PYTHONPATH": os.environ.get("PYTHONPATH"),
        "LD_LIBRARY_PATH": os.environ.get("LD_LIBRARY_PATH"),
        "numerical_threads": 1,
    },
    "scope": "Settings constructors for actual source wavelengths; existing native optical case independently checks1064.2nm formula against unmodified officialIERS Fortran. No input rounding, formula, tolerance or optimization change.",
}
assert record["kernel_sha256"] != record["previous_kernel_sha256"]
(WORK / "logs/optical-1064_2nm-kernel-provenance.json").write_text(
    json.dumps(record, indent=2) + "\n"
)
print(json.dumps(record, indent=2))
