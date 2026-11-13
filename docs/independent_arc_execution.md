# Independent arcs in the native C++ multi-arc solver

The native `MultiArcDynamicsSimulator` now has an opt-in independent-arc path.
`MultiArcVariationalEquationsSolver` exposes the corresponding variational
configuration. Default construction and integration remain sequential.
Python exposes a deliberately restricted native model through
`IndependentArcVariationalSimulator`; see
[the API and ownership contract](python_independent_arcs_api.md). The generic
C++ factory is not exposed as an unrestricted Python callback.

## Configuration and ownership

```cpp
// Construct without integrating on creation.
MultiArcVariationalEquationsSolver<double, double> solver(
    bodies, multiArcSettings, parameters, false);

solver.setIndependentArcExecution(4,
    [](std::size_t arc,
       const std::shared_ptr<SingleArcPropagatorSettings<double, double>>& original,
       const Eigen::VectorXd& currentLocalParameters) {
        // Application-specific: create fresh bodies, force models, caches,
        // integrator, settings, and estimatable parameter objects matching
        // `original`. Set currentLocalParameters in the same parameter order.
        // Return a SingleArcVariationalEquationsSolver with integration OFF.
        return createIsolatedSolver(arc, original, currentLocalParameters);
    });

solver.integrateVariationalAndDynamicalEquations(
    multiArcSettings->getInitialStateList(), true);
// Existing result getters, interpolated matrix interface, and global ephemerides.
solver.clearIndependentArcExecution(); // restore the original sequential path
```

`createIsolatedSolver` is illustrative, not a new built-in Tudat function.
A working Mars implementation is `makeFactory` in
`benchmarks/cpp/nativeMultiArc.cpp`. The lower-level dynamics simulator setter
accepts a factory returning a fresh `SingleArcDynamicsSimulator`, given the arc
index and its initial state/variational matrix.

The factory is required: existing propagator settings contain already-bound
models, which cannot be generically deep-copied. All factory calls are serial,
on each integration, before launching any workers. This includes reintegration
after parameter updates. Six local state parameters and a shared Mars mu are
used in the acceptance benchmark, but the API does not hard-code Mars or six
parameters. It checks parameter dimensions, current values, and non-state
parameter names/order. Shared parameter values are supplied to each isolated
environment without duplicating their global columns.

Factories must reconstruct the same force model, state/body order, frames,
integrator, termination and output settings; they must not integrate on creation.
Known body, gravity, ephemeris, rotation, acceleration, simulator and integrator
aliases are rejected before threads start. These checks do not establish deep
ownership of every pointer or callback in arbitrary custom models. The factory
remains responsible for the complete model graph, callbacks and external library
state. It must not reuse original propagator settings/acceleration objects as a
shortcut. External calls such as SPICE must be synchronized or otherwise
isolated. The acceptance benchmark wraps its Mars input models with a shared
mutex; this is application-specific, not a blanket SPICE lock in Tudat.

## Supported first slice and sequential fallback

Only Cowell translational arcs are enabled for independent execution. Existing
raw-to-output conversions can access model state; other representations remain
sequential until those paths have been isolated. The isolated replacement must
also use Cowell and the same propagated-body and central-body order.
Consider parameters are not supported by the variational factory API yet;
configuring that combination raises an explicit error. Propagated and estimated
state dimensions must match.

Any initial arc state containing NaNs (the existing state-transfer mechanism)
causes the entire call to follow the existing sequential path, without invoking
the factory. This is latched for that solver, including subsequent calls after
transfer NaNs have been resolved; construct a fresh solver to declare a new,
fully independent configuration. Other representations likewise fall back. Hybrid-arc scheduling is
not introduced. Configuring this API explicitly asserts that arcs have no other
dependencies, such as dependence on another arc's newly integrated ephemeris;
those cannot all be inferred automatically from bound/custom callbacks.

At most `min(workers, arcs)` threads execute, with dynamic assignment of complete
arc tasks. `getLastIntegrationUsedIndependentArcs()` on the dynamics simulator
indicates successful use of the isolated path (also true with one configured
worker). Zero workers or a null factory are configuration errors.

## Native results and errors

Workers use `SingleArcDynamicsSimulator::integrateIndependentArc` to write into
the original native per-arc result objects. They finalize numerical metadata but
do not update integrated ephemerides or dependent-variable interpolators on the
worker environment. No manual global matrix history is materialized.

After every worker joins, the existing multi-arc result/status processing and
global ephemeris processing run serially. The variational solver then uses its
existing per-arc interpolator reset and shared-parameter mapping. Result objects
and the combined matrix interface retain their identity across reintegrations.

The original derivative-model objects remain configuration/metadata owners;
their per-evaluation scratch caches are not copies of the worker caches. Use
the native result objects for recorded states and function-evaluation counts.

Worker exceptions and non-successful termination are reported with arc index
after all threads have joined. In that case global postprocessing is skipped;
the call does not mark the multi-arc propagation as performed. Some local arc
results may already exist, but no successful global result is advertised. A
subsequent integration creates fresh isolated workers. The executor continues
already available tasks rather than implementing cancellation.

See [Building and validating independent arcs](independent_arc_development.md)
for environment setup, native/Python acceptance commands and ON/OFF testing.

## Acceptance checks and measurement

```bash
cmake --build build-parallel-arcs --target test_propagators_IndependentArcExecution -j 6
ctest --test-dir build-parallel-arcs -R '^test_propagators_IndependentArcExecution$' \
    --output-on-failure

# Optional separate executable, reusing the configured compiler/link flags:
python benchmarks/cpp/build_prototype.py \
    --source benchmarks/cpp/nativeMultiArc.cpp --name native_multi_arc
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    build-parallel-arcs/benchmarks/native_multi_arc
```

The acceptance test checks every stored state, STM, sensitivity and dependent
variable, plus termination reasons, global ephemerides, shared-column interface
and result/interface identity. It tests changed initial/shared parameters,
clearing the configuration, rejected aliases and stale parameters, worker failure
and recovery, linked-arc fallback on initial propagation and reintegration, and
dynamics-only output. Numerical comparison
uses the existing native sequential solver with the same physical configuration.

The benchmark measures 1, 2, 4 and 8 workers for a fixed workload. Total time
includes original native setup, factory reconstruction, integration and native
postprocessing; validation and destruction are excluded. This is a different
measurement from the earlier standalone prototype's manual merge. Measure this
native route independently rather than copying the prototype's speedup claim.

## Thread-safe Runge–Kutta coefficient tables

The acceptance checks exposed a pre-existing shared-table race. The DP87 getter
tested `higherOrder != 7`, whereas its initializer sets order 8. Consequently it
rewrote a shared static table on every lookup, potentially while another thread
was copying it. Coefficient lookup now uses function-local `static const` tables,
initialized once with C++ thread-safe publication, for all 21 methods. The
initializers and numerical table contents are unchanged. A concurrent cold-start
and repeated-copy check runs before the native propagation tests.

## Local validation and initial native measurement

The native acceptance test passed 30 consecutive executions. The existing
`RungeKuttaCoefficients` and `RungeKutta87DormandPrinceIntegrator` tests passed,
as did the executor's AddressSanitizer/UBSan checks (leak detection disabled
under ptrace). No full-library ThreadSanitizer run has been performed.
To reproduce the existing integrator tests, also build `copy_test_data` so their
MATLAB reference files are present in the build tree. The test-support target
now explicitly propagates its required Cereal header dependency, including when
serialization implementation builds are disabled.

One native scaling run on eight Mars arcs of two days, degree/order 50, step
30 s, with integrated ephemerides and dependent-variable interpolation enabled:

| Workers | Total seconds | Speedup vs 1 |
|---:|---:|---:|
| 1 | 31.070 | 1.00 |
| 2 | 15.615 | 1.99 |
| 4 | 7.977 | 3.90 |
| 8 | 4.717 | 6.59 |

At every worker count, all 46,088 saved arc epochs passed state, STM,
sensitivity, dependent-variable and metadata comparisons against the native
sequential reference. Combined matrix and interpolated dependent-variable
interfaces and global ephemerides also matched. Compilation was stopped and
BLAS/OpenMP thread counts were one. This is an initial local measurement, not a
general mission-analysis speedup guarantee. Raw output is stored in
`benchmarks/cpp/native-scaling-forward.jsonl`.
