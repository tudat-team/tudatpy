#!/usr/bin/env bash
set -eu
export PYTHONPATH=/home/dominic/Tudat/laser-ranging-build/src
export LD_LIBRARY_PATH=/home/dominic/miniconda3/envs/tudatpy-dev/lib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
python_executable=/home/dominic/miniconda3/envs/tudatpy-dev/bin/python
"$python_executable" - <<'PY' > /home/dominic/Tudat/laser-ranging-work/logs/slr-runtime-provenance-attempt2.json
import hashlib, json, sys
from pathlib import Path
import tudatpy
from tudatpy import kernel
expected = Path('/home/dominic/Tudat/laser-ranging-build/src/tudatpy/kernel.so')
assert Path(kernel.__file__).resolve() == expected
assert Path(tudatpy.__file__).resolve().parent == expected.parent
sha = hashlib.sha256(expected.read_bytes()).hexdigest()
assert sha == '1864ed04530210fa056e78128bdbbfe1ad536cae6a9d3ade4f23295fbeda1744'
print(json.dumps({'interpreter':sys.executable,'package':tudatpy.__file__,'kernel':str(expected),'sha256':sha,'source':'2fe4f1cae','script_source':'2d32e3c34','numerical_threads':1},indent=2))
PY
"$python_executable" analysis/laser_ranging/fixed_orbit_slr.py --data-root /home/dominic/Tudat/laser-ranging-work/data/slr --output /home/dominic/Tudat/laser-ranging-work/results/slr-baseline-combined-slrf2020-20261016-attempt2 > /home/dominic/Tudat/laser-ranging-work/logs/slr-baseline-combined-slrf2020-attempt2.log 2>&1
