"""Run the six predeclared fixed-observation LLR controls on one verified kernel."""

from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
from tudatpy import kernel

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--kernel-sha256", required=True)
args = parser.parse_args()
work = Path("/home/dominic/Tudat/laser-ranging-work")
source = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-data")
script = source / "analysis/laser_ranging/llr_residuals.py"
assert hashlib.sha256(Path(kernel.__file__).read_bytes()).hexdigest() == args.kernel_sha256
assert (
    subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()
    == subprocess.check_output(["git", "rev-parse", "cb39aa79d"], cwd=source, text=True).strip()
)
runs = [
    ("de440-inpop17a", "apollo20170309", "inpop17a-gir", "de440"),
    ("de440-matera", "matera20220821", "survey2005", "de440"),
    ("de440-hofmann-apollo", "apollo20170309", "hofmann2017-apollo-standard", "de440"),
    ("de440-hofmann-matera", "matera20220821", "hofmann2017-matera-standard", "de440"),
    ("de421-inpop17a", "apollo20170309", "inpop17a-gir", "de421-mer"),
    ("de421-matera", "matera20220821", "survey2005", "de421-mer"),
]
manifest = {
    "kernel_binary_sha256": args.kernel_sha256,
    "native_source_revision": "7d52716f6",
    "analysis_source_revision": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=source, text=True
    ).strip(),
    "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
    "numerical_threads": {
        name: os.environ.get(name)
        for name in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"]
    },
    "runs": [],
}
keys = {}
for label, case, station, ephemeris in runs:
    output = work / "results" / f"llr-six-{label}-attempt1"
    assert not output.exists()
    kernel_directory = work / (
        "data/llr/kernels" if ephemeris == "de440" else "research/de421-coherent"
    )
    command = [
        sys.executable,
        str(script),
        "--data-dir",
        str(work / "data/llr"),
        "--kernel-dir",
        str(kernel_directory),
        "--output-dir",
        str(output),
        "--case",
        case,
        "--station-profile",
        station,
        "--ephemeris",
        ephemeris,
    ]
    assert hashlib.sha256(Path(kernel.__file__).read_bytes()).hexdigest() == args.kernel_sha256
    with (work / "logs" / f"llr-six-{label}-attempt1.log").open("w") as stream:
        subprocess.run(command, cwd=source, stdout=stream, stderr=subprocess.STDOUT, check=True)
    summary = json.loads((output / "summary.json").read_text())
    assert summary["kernel_binary_sha256"] == args.kernel_sha256
    assert summary["script_sha256"] == manifest["script_sha256"]
    assert summary["count"] == (11 if case == "apollo20170309" else 2)
    assert summary["ocean_loading"] is None
    assert summary["inertial_axes"]["environment_frame"] == "ICRF"
    with (output / "normal_points.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    current = [(row["index"], row["utc_launch"], row["reflector"]) for row in rows]
    if case in keys:
        assert current == keys[case]
    else:
        keys[case] = current
    record = {
        "label": label,
        "case": case,
        "station_profile": summary["station_profile"],
        "ephemeris": ephemeris,
        "output": str(output),
        "count": summary["count"],
        "mean_m": summary["mean_m"],
        "rms_about_zero_m": summary["rms_about_zero_m"],
        "scatter_about_mean_m": summary["scatter_about_mean_m"],
    }
    manifest["runs"].append(record)
    (work / "logs/llr-six-controls-provenance.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(json.dumps(record), flush=True)
