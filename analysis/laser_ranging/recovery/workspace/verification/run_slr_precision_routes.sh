#!/usr/bin/env bash
set -eu
export PYTHONPATH=/home/dominic/Tudat/laser-ranging-build/src
export LD_LIBRARY_PATH=/home/dominic/miniconda3/envs/tudatpy-dev/lib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
slr_python=/home/dominic/miniconda3/envs/tudatpy-dev/bin/python
"$slr_python" - <<'PY' > /home/dominic/Tudat/laser-ranging-work/logs/slr-precise-runtime-provenance-attempt1.json
import hashlib,json,sys
from pathlib import Path
import tudatpy
from tudatpy import kernel
expected=Path('/home/dominic/Tudat/laser-ranging-build/src/tudatpy/kernel.so')
assert Path(kernel.__file__).resolve()==expected
sha=hashlib.sha256(expected.read_bytes()).hexdigest()
assert sha=='bd8d893574bd29281188934ff5a62647e4c08494eb4bfb846c17cfc06a062968'
print(json.dumps({'interpreter':sys.executable,'package':tudatpy.__file__,'kernel':str(expected),'sha256':sha,'binding_source':'8506f918874a2a628aa200cdab579e0f85e88a90','script_source':'0a8f805a5','numerical_threads':1},indent=2))
PY
for slr_route in legacy cubic30-double cubic30-time cubic60-time; do
 case "$slr_route" in
  legacy) slr_args=(--body-ephemeris direct --satellite-epoch double) ;;
  cubic30-double) slr_args=(--body-step 30 --satellite-epoch double) ;;
  cubic30-time) slr_args=(--body-step 30 --satellite-epoch time) ;;
  cubic60-time) slr_args=(--body-step 60 --satellite-epoch time) ;;
 esac
 "$slr_python" analysis/laser_ranging/fixed_orbit_slr.py --data-root /home/dominic/Tudat/laser-ranging-work/data/slr --output "/home/dominic/Tudat/laser-ranging-work/results/slr-precise-$slr_route-attempt1" "${slr_args[@]}" > "/home/dominic/Tudat/laser-ranging-work/logs/slr-precise-$slr_route-attempt1.log" 2>&1
 printf '%s completed\n' "$slr_route"
done
