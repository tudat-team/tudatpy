# Spherical-harmonic gravity through the Python API

The driver `benchmark_python_spherical_harmonics.py` uses the ordinary Python
acceleration factory and native solvers. There are no replacement objects,
experimental acceleration subclasses or global gravity switches. Each case
owns an independent environment, acceleration model, cache and solver. Buffers
are reused within that model. Team counts are requested upper bounds, not
measurements of actual OpenMP team sizes.

## Run

Build the current native libraries and Python kernel consistently with
`TUDAT_BUILD_WITH_OPENMP=ON`; activate the configured Python environment and
select that build's package, for example:

```bash
export PYTHONPATH="$PWD/build-gravity-native-on/src"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
python benchmarks/benchmark_python_spherical_harmonics.py \
  --workers 1 2 4 8 --degree 100 --duration 86400 --repeats 3 \
  --output local-profiles/gravity-dynamics.jsonl
python benchmarks/benchmark_python_spherical_harmonics.py \
  --workers 1 2 4 8 --degree 100 --duration 21600 --mode variational --repeats 3 \
  --output local-profiles/gravity-variational.jsonl
```

Use new output paths; existing results are not overwritten. A short API example
is the same command with `--workers 1 4 --degree 20 --duration 300 --repeats 1`.
This short case is for learning/validation, not evidence of speedup.

Requires standard SPICE kernels and JGMRO120D. The Mars-only environment does not
require the missing historical Earth-orientation resource. Its physical setup
matches the earlier experiment: degree/order 100, 300-km near-circular Mars
orbit, inclination 93 degrees, RKDP87 fixed 30-second steps, start epoch
778118400. Variational propagation estimates the six initial-state components
and Mars's gravitational parameter. Hessian summation remains sequential.

Setup includes creating each fresh environment and solver. Integration is timed
separately; result extraction, finite/shape/epoch validation, comparison and
teardown are outside the measured interval. A one-thread reference is untimed.
Order rotates across three repetitions. There is no fixed speedup test assertion.

## Measurements, 15 November 2026

AMD Ryzen 9 7900X, 12 physical cores / 24 logical CPUs, GCC 14.4, Release `-O3`,
Python 3.14.6, no CPU pinning. Timing build: OpenMP ON, Cereal OFF, extended
precision OFF. No compiler or other benchmark was running during these cases.
All speedups compare medians within the same build and workload.

| Workload | Requested threads | Integration median (range), s | Speedup | Setup + integration median, s |
|---|---:|---|---:|---:|
| Dynamics, one day | 1 | 3.947 (3.933–4.002) | 1.00 | 3.956 |
| Dynamics, one day | 2 | 2.673 (2.624–2.678) | 1.48 | 2.682 |
| Dynamics, one day | 4 | 1.640 (1.616–1.731) | 2.41 | 1.649 |
| Dynamics, one day | 8 | 1.710 (1.626–1.897) | 2.31 | 1.719 |
| Variational, six hours | 1 | 1.769 (1.767–1.773) | 1.00 | 1.779 |
| Variational, six hours | 2 | 1.422 (1.414–1.427) | 1.24 | 1.431 |
| Variational, six hours | 4 | 1.140 (1.138–1.184) | 1.55 | 1.149 |
| Variational, six hours | 8 | 1.145 (1.138–1.181) | 1.55 | 1.154 |

All 24 timed propagations pass comparison at all 2,881 dynamics epochs or 721
variational epochs, including the complete 6x6 STM and 6x1 sensitivity. All
observed maximum absolute differences are zero. The unchanged comparison is
`abs(actual-reference) <= 1e-14 + 1e-11*abs(reference)`; finite values, interval,
epoch count and shapes are checked separately.

An older complete Python kernel was run separately, with the new thread keyword
removed but all physical/integration settings retained. Three sequential medians
were 4.005 s for dynamics and 1.799 s for variational propagation. Its full
histories also match the current build with zero observed difference. This is
a historical control from a different existing build, not an otherwise identical
fresh build of the parent commit. It offers no evidence that the new speedup is
caused by a slowed sequential reference, but is not a portable regression bound.

The earlier standalone C++ timings used a different executable/build layout;
their 2.17x factor and absolute seconds should not be combined with these values.
Four threads are a sensible measured choice here; eight add no clear benefit.
The short degree-20 smoke case was slower at four threads (about 1.81 ms versus
1.21 ms), reinforcing explicit opt-in rather than automatic parallel defaults.

Raw JSONL, run metadata, build/test logs, historical-control script and result
snapshots are retained locally in `local-profiles/20261115/python-gravity/`.
The related install/API documentation is
[installation_openmp.rst](../docs/tudatpy/source/installation_openmp.rst).
