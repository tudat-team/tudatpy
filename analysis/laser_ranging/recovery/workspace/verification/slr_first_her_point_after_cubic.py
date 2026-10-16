"""Evaluate only the original failing HER point with the exact historical source/helper."""

import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

from tudatpy import kernel

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--kernel-sha256", required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
assert hashlib.sha256(Path(kernel.__file__).read_bytes()).hexdigest() == args.kernel_sha256
root = Path("/home/dominic/Tudat/laser-ranging-work")
source = root / "verification/slr-radial-historical/fixed_orbit_slr.py"
sys.path.insert(0, str(source.parent))
spec = importlib.util.spec_from_file_location("historical_slr", source)
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)
original_reduce = analysis.reduce_session
original_extract = analysis.slr.extract_normal_point_measurements


def first_point_only(session, index, *remaining):
    """Keep the original first HER observable and its weather/physics unchanged."""
    if index != 29:
        return [], []
    analysis.slr.extract_normal_point_measurements = lambda selected: original_extract(selected)[:1]
    try:
        return original_reduce(session, index, *remaining)
    finally:
        analysis.slr.extract_normal_point_measurements = original_extract


analysis.reduce_session = first_point_only
sys.argv = [
    str(source),
    "--data-root",
    str(root / "data/slr"),
    "--output",
    str(args.output),
    "--got410c-harpos",
    str(root / "research/ocean-loading/toc_got410c_harmod.hps"),
]
analysis.main()
with (args.output / "residuals.csv").open(newline="") as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 1 and (rows[0]["station"], rows[0]["session"], rows[0]["normal_point"]) == (
    "7840",
    "29",
    "0",
)
print(
    json.dumps(
        {
            "diagnostic": "Original first HER failure now evaluated without changing50/1e-12 or physical inputs",
            "kernel_sha256": args.kernel_sha256,
            "historical_script_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "point": rows[0],
        },
        indent=2,
    )
)
