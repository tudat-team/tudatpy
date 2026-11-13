# C++ prototype: parallel independent arcs

The native C++ integration is now available separately; see
[independent_arc_execution.md](../../docs/independent_arc_execution.md) and
`nativeMultiArc.cpp`. The standalone prototype below remains a historical
baseline. Its shared Mars setup helpers have moved to `marsArcBenchmark.h`;
the executor is now a library header with a benchmark compatibility include.

`parallelMultiArc.cpp` is an experimental executable, not a change to Tudat's
public solver API. It compares the existing native sequential multi-arc solver
with a C++ thread executor running independently constructed single-arc solvers.
It uses the same Mars configuration and arc epochs as `profile_multi_arc.py`.
Each arc has six local initial-state parameters and the same Mars mu value.

## Isolation and scope

All environments and solvers are built serially before launching threads.
Every arc owns fresh bodies, force models, gravity/Legendre caches, parameter
objects and result containers. There is no shallow copy of SystemOfBodies.
The coefficient files are read independently during setup, which is conservative
and makes setup/memory costs visible. This prototype allocates one environment
per arc, rather than reusing one per worker.

Mars ephemeris/rotation evaluations delegate to the original models under one
shared mutex. This serializes calls into SPICE; loading/unloading kernels and
model construction happen only before threads launch. These wrappers cover this
restricted Mars-only configuration, **not arbitrary Tudat models** or all future
SPICE entry points. The reference uses the original unwrapped models.

`isolatedArcExecutor.h` schedules complete arc integrations with a dynamic task
index. Results retain their arc index regardless of completion order. Task
exceptions are collected and rethrown after every thread has joined. The lowest
failed arc index is reported deterministically. Dependent arcs/state transfer,
hybrid arcs and observation/estimation processing are outside this prototype.

## Merge and validation

For each saved epoch, the local 6x6 STM is scattered into columns 6i..6i+5
of a 6x(6N+1) global matrix, and local mu sensitivity into the last column.
Other arc initial-state columns are zero. Full global histories are materialized
in the serial merge phase. This is intentionally simple and memory intensive;
a production implementation should preserve sparse/arc-local storage and expose
the existing multi-arc interpolation interface rather than duplicate zero blocks.

Every saved state and combined variational matrix is compared against the native
sequential multi-arc solver at the same epoch. Endpoint matrices are additionally
compared against its actual global state-transition/sensitivity interface.
Tolerance is componentwise `1e-14 + 1e-11 * abs(reference)`; finiteness, dimensions,
history lengths and end epochs are checked. Validation is outside timed work.
The executor also checks task coverage and injected exception propagation.

The standalone executor test can be built independently of Tudat:

```bash
c++ -std=c++17 -O2 -pthread -I include -fsanitize=address,undefined \
    benchmarks/cpp/testIsolatedArcExecutor.cpp -o /tmp/tudat-test-arc-executor
ASAN_OPTIONS=detect_leaks=0 /tmp/tudat-test-arc-executor
```

AddressSanitizer/UBSan checks passed. Leak detection is disabled for this run
because LeakSanitizer cannot operate under the execution environment's ptrace
supervision; this is not a leak check or a ThreadSanitizer test of Tudat.

This does not yet transplant results into a production MultiArcSimulationResults,
update a caller's ephemerides, perform a joint fit, or expose Python controls.
It is a measured implementation of the proposed isolated-arc task boundary.

## Build

The local build helper reuses compiler, defines, includes and link libraries from
the configured `test_propagators_MultiArcVariationalEquations` CMake target.
That target and its Tudat libraries must have been built from compatible source.
No source-library ABI compatibility can be inferred just from a successful link.

For this session a separate `build-parallel-arcs` directory is used because the
original `build` libraries predate the current headers. Configuration uses the
`tudatpy-dev` compiler/dependencies, Release, serialization off and extended
precision off. Cereal 1.3.2 was downloaded and installed only under `/tmp` to
provide the missing header dependency. No installation into the Conda environment
is performed. For a fresh configure, supply your own dependency prefix and
compiler, for example:

```bash
cmake -S . -B build-parallel-arcs \
    -DCMAKE_CXX_COMPILER="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-c++" \
    -DCMAKE_C_COMPILER="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-cc" \
    -DCMAKE_PREFIX_PATH="$CONDA_PREFIX;/tmp/tudat-cereal-prefix" \
    -DCMAKE_BUILD_TYPE=Release -DTUDAT_BUILD_TESTS=ON \
    -DTUDAT_BUILD_WITH_SERIALIZATION=OFF \
    -DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=OFF
cmake --build build-parallel-arcs --target test_propagators_MultiArcVariationalEquations -j 6
```

The named test target pulls in the library set needed by the prototype. The
executable explicitly uses `double` time and state types, with a native reference
using the same types; absolute timings should not be equated with those of a
different installed Python build.

```bash
python benchmarks/cpp/build_prototype.py --build-dir build-parallel-arcs
# If Cereal isn't on the compiler's normal include path:
python benchmarks/cpp/build_prototype.py --build-dir build-parallel-arcs \
    --cereal-include /tmp/cereal-1.3.2/include
```

The helper writes only `<build-dir>/benchmarks/parallelMultiArc.o` and
`<build-dir>/benchmarks/parallel_multi_arc`; it does not install a package or replace
the Python environment. Compilation uses the configured optimization flags,
debug information, and pthread support. Cereal is a dependency of current Tudat
headers, not of the thread scheduler itself.

## Run

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
build-parallel-arcs/benchmarks/parallel_multi_arc --self-test
build-parallel-arcs/benchmarks/parallel_multi_arc --arcs 4 --degree 20 --days-per-arc 0.02

# Fixed workload; measures 1, 2, 4 and 8 workers in the same executable.
build-parallel-arcs/benchmarks/parallel_multi_arc --arcs 8 --degree 50 --days-per-arc 2
```

JSON output separates serial setup, threaded integration, serial merge and total
wall time. Total time excludes validation, reference data copies, console output
and destruction. Keep this definition in mind: it is not whole-process time.
The native reference time includes its own setup, integration and native
postprocessing; it does not materialize the prototype's global matrix histories.
For scaling, compare prototype worker counts with each other, not just with the
native reference. Repeat runs on an otherwise quiet machine before claiming a
speedup, and monitor memory as arc count grows.

Use `--workers 8,4,2,1` to vary run order, or a single worker count to profile a
particular configuration. At most min(arcs, workers) threads are started.

## Measurements, 2026-11-08

The first complete run, using eight arcs of two days each, degree/order 50 and
step 30 s, gave the following times. BLAS/OpenMP were set to one thread and
compilation was stopped before measurement. These are prototype total times
with the definition above, not whole-process timings.

| Workers | Setup (s) | Integration (s) | Merge (s) | Total (s) | Speedup vs 1 |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.056 | 29.202 | 0.056 | 29.314 | 1.00 |
| 2 | 0.054 | 14.729 | 0.015 | 14.799 | 1.98 |
| 4 | 0.054 | 7.457 | 0.015 | 7.526 | 3.90 |
| 8 | 0.055 | 4.502 | 0.014 | 4.571 | 6.41 |

The native sequential reference took 29.467 s including setup/native
postprocessing. All 46,088 stored arc epochs passed the state and variational
matrix comparisons at every worker count. Endpoint global matrices also
matched the native multi-arc interface. Raw output: `scaling-forward.jsonl`.
The shorter four-arc, degree-20 check also passed for 1, 2, 4 and 8 requested
workers (the last uses only four active threads).

A second run in reverse worker order also passed every comparison. Totals were
29.175 s (1 worker), 14.681 s (2), 7.490 s (4), and 4.467 s (8), giving 6.53x
speedup at eight workers relative to that run's one-worker result. Its native
reference took 29.201 s. Raw output: `scaling-reverse.jsonl`. Two runs support
this preliminary scaling result; they do not establish general mission-analysis
performance or thread safety for other models.

The local source/library mismatch required a fresh library build; no result
from the incompatible original build is used. All required libraries built
successfully in `build-parallel-arcs`. After that the additional native unit-test
compilation was stopped to avoid competing with the benchmark. That existing
unit test was not executed; the prototype performs its own full-history native
reference comparison. No ThreadSanitizer-instrumented Tudat build was run.
