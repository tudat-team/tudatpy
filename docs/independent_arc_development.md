# Building and validating independent arcs

This guide covers the native independent-arc route and its restricted Python
facade. Read the [native design](independent_arc_execution.md) and
[Python API contract](python_independent_arcs_api.md) before extending either.
Keep sequential defaults, isolated per-arc models, native result interfaces and
sequential fallback for coupled arcs.

## Environment and resources

Use a development environment with the dependencies in `environment.yaml`.
That file is not a version lock. Check CMake, C/C++ compilers, Eigen, Boost,
Cereal, CSPICE, SOFA, NRLMSISE, Tudat Resources, pybind11, NumPy and pytest.
Cereal headers are also needed by native test support when serialization is OFF.
Public Python imports require the declared `atdf2ascii` dependency and its
transitive dependencies. Do not substitute an older installed kernel for a
fresh build.

The example needs standard SPICE kernels and the real
`gravity_models/Mars/jgmro120d.txt` resource. Load kernels before creating models;
workers must not load/unload kernels. Do not replace gravity data with synthetic
coefficients for acceptance checks.

For reference, validation used Python 3.14.6, GCC 14.4.0, pybind11 3.0.3,
Cereal 1.3.2 and Tudat Resources 2.4. Native reproduction also used Eigen 5.0.1,
Boost 1.85.0 and CSPICE 67. These are tested versions, not a portability guarantee.
Install dependencies in a durable prefix; no temporary local prefix is required.

## Separate serialization-OFF and ON builds

From the repository root, activate the development environment. These compiler
paths assume a Linux Conda toolchain; select equivalent compilers elsewhere.
Use fresh build directories and adjust jobs to available memory. Do not copy
libraries between ON and OFF configurations.

```bash
conda activate tudatpy-dev
cmake -S . -B build-python-arcs -G 'Unix Makefiles' \
  -DCMAKE_CXX_COMPILER="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-c++" \
  -DCMAKE_C_COMPILER="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-cc" \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DCMAKE_BUILD_TYPE=Release -DTUDAT_BUILD_TESTS=ON \
  -DTUDAT_BUILD_WITH_SERIALIZATION=OFF \
  -DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=OFF \
  -DTUDAT_BUILD_EXPLICIT_INSTANTIATIONS=ON \
  -DCCACHE_PROGRAM=CCACHE_PROGRAM-NOTFOUND
cmake --build build-python-arcs --target kernel -j 2

cmake -S . -B build-python-arcs-serialization-on -G 'Unix Makefiles' \
  -DCMAKE_CXX_COMPILER="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-c++" \
  -DCMAKE_C_COMPILER="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-cc" \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DCMAKE_BUILD_TYPE=Release -DTUDAT_BUILD_TESTS=ON \
  -DTUDAT_BUILD_WITH_SERIALIZATION=ON \
  -DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=OFF \
  -DTUDAT_BUILD_EXPLICIT_INSTANTIATIONS=ON \
  -DCCACHE_PROGRAM=CCACHE_PROGRAM-NOTFOUND
cmake --build build-python-arcs-serialization-on --target kernel -j 2
```

Check the actual loaded module and compiled feature flag in separate Python
processes. Do not load both kernels in one interpreter.

```bash
PYTHONPATH=build-python-arcs/src python -c \
  'from tudatpy import kernel; assert kernel._build_with_serialization is False; print(kernel.__file__)'
PYTHONPATH=build-python-arcs-serialization-on/src python -c \
  'from tudatpy import kernel; assert kernel._build_with_serialization is True; print(kernel.__file__)'
```

## Python checks

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
PYTHONPATH=build-python-arcs/src python benchmarks/python_independent_arcs.py
PYTHONPATH=build-python-arcs/src python -m pytest \
  tests/test_tudatpy/test_independent_arc_execution.py \
  tests/test_tudatpy/test_python_arc_benchmark.py --no-remote-data -q
PYTHONPATH=build-python-arcs/src python -m pytest tests/test_tudatpy --no-remote-data -q

PYTHONPATH=build-python-arcs-serialization-on/src python -m pytest \
  tests/test_tudatpy/test_serialization_file_io.py \
  tests/test_tudatpy/test_serialization_pickle.py \
  tests/test_tudatpy/test_serialization_build_policy.py --no-remote-data -q
PYTHONPATH=build-python-arcs-serialization-on/src python -m pytest \
  tests/test_tudatpy --no-remote-data -q
PYTHONPATH=build-python-arcs-serialization-on/src python benchmarks/python_independent_arcs.py
```

The short example should report independent execution, 21 saved states in each
of two arcs and a `(6, 13)` global matrix. The loaded binary's explicit
`_build_with_serialization` flag controls Cereal test selection; missing or
non-boolean flags are configuration errors. `Time` tuple-pickle tests remain
active in OFF builds. Network tests, including GRAIL/MRO downloads, respect
`--no-remote-data`. Do not use manual module or name exclusions. See the
[Python test guide](../tests/test_tudatpy/README.md).

A policy test with an artificial ON flag does not exercise Cereal: actual
roundtrips require the ON binary. Test settings and results through JSON,
binary and pickle routes, including histories, polymorphism and invalid loaders.
For polymorphic file tests, write and read through the same declared root type.
Some native derived types have no dedicated Python binding; check virtual C++
equality in both directions rather than assuming identical Python wrappers.
Explicitly bound result types should also retain type and history checks.

Result equality must compare owned settings/interfaces by value, with null
handling, rather than pointer identity after deserialization. Default print/save
intervals use NaN for "unset": only those configuration sentinels treat two
unset values as equal. Numerical histories still require finite values. Native
`test_DefaultProcessingSettingsEquality` checks defaults and finite-vs-unset
comparisons in both directions. Serialization fixtures use Tudat's own SPICE
module and the current public propagation API, not an unrelated `spice` package.
Two observation-collection file tests remain explicitly skipped for refactoring;
the native observation-collection target below is exercised separately.

## Native acceptance and serialization checks

```bash
cmake --build build-python-arcs --target copy_test_data \
  test_propagators_IndependentArcExecution test_integrators_RungeKuttaCoefficients \
  test_integrators_RungeKutta87DormandPrinceIntegrator -j 2
ctest --test-dir build-python-arcs \
  -R '^test_(propagators_IndependentArcExecution|integrators_RungeKuttaCoefficients|integrators_RungeKutta87DormandPrinceIntegrator)$' \
  --output-on-failure
ctest --test-dir build-python-arcs \
  -R '^test_propagators_IndependentArcExecution$' --repeat until-fail:30
ctest --test-dir build-python-arcs \
  -R '^test_integrators_RungeKuttaCoefficients$' --repeat until-fail:20

cmake --build build-python-arcs-serialization-on --target copy_test_data \
  test_io_SettingsSerialization test_io_ObservationCollectionSerialization \
  test_io_PropagationResultsSerialization test_io_PodInputOutputTypesSerialization -j 2
ctest --test-dir build-python-arcs-serialization-on \
  -R '^test_io_(SettingsSerialization|ObservationCollectionSerialization|PropagationResultsSerialization|PodInputOutputTypesSerialization)$' \
  --output-on-failure
```

`copy_test_data` provides the integrator reference files. Native acceptance uses
`double` time, whereas the actual Python bindings use `tudat::Time`; native tests
alone do not replace Python tests. Repeated acceptance checks are useful for
concurrency regression testing, but do not establish absence of races. A full
ThreadSanitizer run and free-threaded Python validation remain outstanding.

## Validation evidence and limits

On 12 November 2026, the separate ON build passed 235 focused Python tests
(two refactoring skips), all four native serialization targets, the short
example and the full offline suite (915 passed, 108 skipped). The preserved
OFF binary passed the full offline suite (694 passed, 329 skipped); its file
fingerprint was unchanged. The new C++ equality fixes were compiled in ON;
the retained OFF binary was not rebuilt with those fixes. Both full-suite
totals include the 28 scaling-driver tests. The fixes and benchmark are recorded
in `eba2dca47` and `bc682fefc`. Counts are historical evidence, not fixed test
assertions or a claim that all native Tudat tests were run.

Native arc/RK checks were also reproduced on an Intel Core i7-8550U on
9 November 2026, before Python bindings existed: three targets passed, with
30 acceptance and 20 coefficient-test repetitions. The small four-arc timing
run does not establish scaling. The [native design](independent_arc_execution.md)
and [Python benchmark guide](../benchmarks/PYTHON_SCALING.md) retain the separate
performance evidence and its measurement boundaries.

Successful roundtrips within one build/schema family do not establish arbitrary
cross-version or cross-platform archive compatibility. Retained results are live
objects: follow the [serialization contract](python_independent_arcs_api.md#serialization-contract-for-independent-arcs)
and save serially only after successful integration and worker joins. They are
not resumable checkpoints. Automated checks and these AI-assisted instructions
do not replace independent scientific or code review.
