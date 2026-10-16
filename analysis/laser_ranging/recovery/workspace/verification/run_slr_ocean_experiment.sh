#!/usr/bin/env bash
set -eu
export PYTHONPATH=/home/dominic/Tudat/laser-ranging-build/src
export LD_LIBRARY_PATH=/home/dominic/miniconda3/envs/tudatpy-dev/lib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
slr_python=/home/dominic/miniconda3/envs/tudatpy-dev/bin/python
"$slr_python" - <<'PY' > /home/dominic/Tudat/laser-ranging-work/logs/slr-ocean-runtime-provenance-attempt1.json
import hashlib,json,sys
from pathlib import Path
import tudatpy
from tudatpy import kernel
expected=Path('/home/dominic/Tudat/laser-ranging-build/src/tudatpy/kernel.so')
assert Path(kernel.__file__).resolve()==expected
sha=hashlib.sha256(expected.read_bytes()).hexdigest()
assert sha=='bd8d893574bd29281188934ff5a62647e4c08494eb4bfb846c17cfc06a062968'
print(json.dumps({'interpreter':sys.executable,'package':tudatpy.__file__,'kernel':str(expected),'sha256':sha,'binding_source':'8506f918874a2a628aa200cdab579e0f85e88a90','script_source':'74c207ff3','numerical_threads':1},indent=2))
PY
"$slr_python" analysis/laser_ranging/fixed_orbit_slr.py --data-root /home/dominic/Tudat/laser-ranging-work/data/slr --output /home/dominic/Tudat/laser-ranging-work/results/slr-ocean-noload-control-attempt1 > /home/dominic/Tudat/laser-ranging-work/logs/slr-ocean-noload-control-attempt1.log 2>&1
printf 'no-loading control completed\n'
"$slr_python" analysis/laser_ranging/fixed_orbit_slr.py --data-root /home/dominic/Tudat/laser-ranging-work/data/slr --output /home/dominic/Tudat/laser-ranging-work/results/slr-ocean-harpos-attempt1 --got410c-harpos /home/dominic/Tudat/laser-ranging-work/research/ocean-loading/toc_got410c_harmod.hps > /home/dominic/Tudat/laser-ranging-work/logs/slr-ocean-harpos-attempt1.log 2>&1
printf 'HARPOS companion completed\n'
