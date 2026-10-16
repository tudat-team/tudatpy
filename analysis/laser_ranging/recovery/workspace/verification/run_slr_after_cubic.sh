#!/usr/bin/env bash
set -eu
slr_expected_kernel_sha="$1"
export PYTHONPATH=/home/dominic/Tudat/laser-ranging-build/src
export LD_LIBRARY_PATH=/home/dominic/miniconda3/envs/tudatpy-dev/lib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
slr_python=/home/dominic/miniconda3/envs/tudatpy-dev/bin/python
slr_work=/home/dominic/Tudat/laser-ranging-work
slr_historical="$slr_work/verification/slr-radial-historical/fixed_orbit_slr.py"
"$slr_python" - "$slr_expected_kernel_sha" <<'PY' > /home/dominic/Tudat/laser-ranging-work/logs/slr-after-cubic-runtime-provenance-attempt1.json
import hashlib,json,sys
from pathlib import Path
from tudatpy import kernel
path=Path(kernel.__file__).resolve()
assert path==Path('/home/dominic/Tudat/laser-ranging-build/src/tudatpy/kernel.so')
sha=hashlib.sha256(path.read_bytes()).hexdigest()
assert sha==sys.argv[1]
print(json.dumps({'kernel':str(path),'sha256':sha,'interpreter':sys.executable,'script_source':'9a912e85ee38d07e235f93bffcf8578a395b91d5','numerical_threads':1},indent=2))
PY
"$slr_python" "$slr_work/verification/slr_first_her_point_after_cubic.py" --kernel-sha256 "$slr_expected_kernel_sha" --output "$slr_work/results/slr-after-cubic-first-her-radial-attempt1" > "$slr_work/logs/slr-after-cubic-first-her-radial-attempt1.log" 2>&1
printf 'original HER point completed\n'
"$slr_python" "$slr_historical" --data-root "$slr_work/data/slr" --output "$slr_work/results/slr-after-cubic-historical-noload-attempt1" > "$slr_work/logs/slr-after-cubic-historical-noload-attempt1.log" 2>&1
printf 'historical no-load completed\n'
"$slr_python" "$slr_historical" --data-root "$slr_work/data/slr" --output "$slr_work/results/slr-after-cubic-historical-radial-attempt1" --got410c-harpos "$slr_work/research/ocean-loading/toc_got410c_harmod.hps" > "$slr_work/logs/slr-after-cubic-historical-radial-attempt1.log" 2>&1
printf 'historical radial completed\n'
"$slr_python" analysis/laser_ranging/fixed_orbit_slr.py --data-root "$slr_work/data/slr" --output "$slr_work/results/slr-after-cubic-current-noload-attempt1" > "$slr_work/logs/slr-after-cubic-current-noload-attempt1.log" 2>&1
printf 'current no-load completed\n'
"$slr_python" analysis/laser_ranging/fixed_orbit_slr.py --data-root "$slr_work/data/slr" --output "$slr_work/results/slr-after-cubic-current-geodetic-attempt1" --got410c-harpos "$slr_work/research/ocean-loading/toc_got410c_harmod.hps" > "$slr_work/logs/slr-after-cubic-current-geodetic-attempt1.log" 2>&1
printf 'current geodetic completed\n'
"$slr_python" - <<'PY'
import csv,json
from pathlib import Path
root=Path('/home/dominic/Tudat/laser-ranging-work/results')
def read(name):
 with (root/name/'residuals.csv').open(newline='') as stream:rows=list(csv.DictReader(stream))
 return rows,json.loads((root/name/'summary.json').read_text())
reference,old=read('slr-precise-cubic30-time-attempt1')
keys=lambda rows:[(r['station'],r['session'],r['normal_point']) for r in rows]
for name in ('historical-noload','historical-radial','current-noload','current-geodetic'):
 rows,s=read('slr-after-cubic-'+name+'-attempt1')
 assert len(rows)==23 and keys(rows)==keys(reference) and s['excluded']==old['excluded']
historical,_=read('slr-after-cubic-historical-noload-attempt1')
current,_=read('slr-after-cubic-current-noload-attempt1')
assert historical==current
print('All first-week routes retain23 original keys and6 exclusions; date parameterization reproduces same-kernel historical no-load CSV exactly.')
PY
"$slr_python" analysis/laser_ranging/fixed_orbit_slr.py --week-ending 2024-04-13 --data-root "$slr_work/data/slr-week-20240407-13" --output "$slr_work/results/slr-week2-noload-attempt1" > "$slr_work/logs/slr-week2-noload-attempt1.log" 2>&1
printf 'second frozen week no-load completed\n'
"$slr_python" analysis/laser_ranging/fixed_orbit_slr.py --week-ending 2024-04-13 --data-root "$slr_work/data/slr-week-20240407-13" --output "$slr_work/results/slr-week2-geodetic-attempt1" --got410c-harpos "$slr_work/research/ocean-loading/toc_got410c_harmod.hps" > "$slr_work/logs/slr-week2-geodetic-attempt1.log" 2>&1
printf 'second frozen week geodetic completed\n'
