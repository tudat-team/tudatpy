# Parallel spherical-harmonic gravity: implementation contract

Status: proposed production contract, 14 November 2026. The current branch
contains feasibility experiments, not an implementation of this API. The
[acceleration/Legendre experiment](../benchmarks/inner_loop/GRAVITY_BLOCK.md)
provides the measured motivation: approximately 1.85x dynamics-only speedup
on a one-day degree/order-100 Mars workload with four threads. Variational
propagation measured a smaller gain in its separate workload. These measurements
are conditional on the documented reused-library setup and are not acceptance
evidence for a fresh production build.

## Scope and branch boundaries

The first implementation covers the spherical-harmonic acceleration sum and
its Legendre values/derivatives, for the existing translational gravity model.
It preserves the force model, coefficients, frames, integrator, output and term
accumulation order. Hessian parallelization remains an isolated reference
experiment and is not enabled as part of this feature.

`experiment/parallel-spherical-harmonics` is a sufficiently specific branch
name for this work. The experiments and first production integration can remain
on that branch; changing its prefix to `feature/` later is optional bookkeeping.
It is based on develop and does not require the independent-arc PR to function.

Panelled aerodynamics and radiation pressure are a separate force family and
should receive a separate experimental branch, for example
`experiment/parallel-panel-forces`. Choose its base from the actual dependencies
at that time. Prefer develop when independent; only stack it on this work if it
really reuses an unmerged implementation. Do not combine unrelated force-family
changes into one gravity PR. Free-threaded Python is outside both tasks.

## User-facing configuration

Use an explicit, model-local number of threads, defaulting to one. The intended
Python shape is an additional optional keyword on the existing factory:

```python
# Proposed API, not yet implemented:
gravity = propagation_setup.acceleration.spherical_harmonic_gravity(
    maximum_degree=100,
    maximum_order=100,
    number_of_threads=4,
)
```

Existing calls, including the point-mass-removal option, retain their meaning.
The C++ settings constructor/factory receives a matching optional thread count
without reordering current parameters. Reject zero/negative counts. Do not
silently select all CPUs or promise that a requested team is always available.
The requested count is an upper bound subject to runtime resource limits;
benchmarks should record actual team size separately.

The acceleration factory transfers this policy into each newly created native
model and its own spherical-harmonics/Legendre caches. No `omp_set_num_threads`,
process-global enable switches or diagnostic counters belong in the public
implementation. Changing one model must not change another model's policy.

Thread count is a local execution policy, not a physical estimated parameter.
The first implementation should keep the existing scientific archive schema:
execution policy is not persisted in old acceleration-setting archives and
loads default to one thread. Document this explicitly on the factory/property
and cover it in roundtrip tests. Equality of the scientific settings retains
its existing meaning. Persisting execution policy would require a separately
designed archive-format/versioning change; do not silently append a binary
field or pretend old archives contain it.

## Ownership and numerical execution

Each model/cache owns its execution configuration and scratch storage. Copies
must have independent mutable scratch buffers. There is no shared mutable
scratch table between accelerations, bodies, arcs or solvers. As in the existing
API, a single mutable model/cache cannot be evaluated concurrently by unrelated
callers without external synchronization.

Retain the established dependency order:

1. Set the latitude input and compute sectoral diagonal values serially.
2. Compute independent vertical recurrences by order, retaining each column's
   degree order. Finish all columns before derivative work.
3. Compute first derivatives by independent rows and complete that pass.
4. Compute second derivatives only when required, after the first derivatives.
5. Complete the spherical cache update before acceleration-term workers read it.
6. Compute independent gradient terms into model-owned storage, then add them
   in the original degree/order sequence on the caller.

Reuse scratch allocation across evaluations; resize safely when coefficient
limits change. Preserve finite-value/shape checks and the original near-pole
and saved-per-term-output behaviour. Initial production support may retain
serial fallbacks for those branches, but must document and test them.
No SPICE, Python callbacks or numerical-integration steps run in the inner team.

Do not infer a universal workload threshold from measurements at degrees 50
and 100 alone. Start with explicit opt-in. Any later automatic threshold needs
its own measured crossover and documented semantics.

## Failure handling and parallel composition

Every worker exception must be captured inside the parallel region, all work
must finish/join, and the caller must receive the failure afterwards. Never
allow a C++ exception to escape an OpenMP worker region. Preallocate required
scratch/error storage on the caller where possible. A failed cache update must
invalidate its cached-input marker so a retry cannot accept partially computed
values as complete.

One thread follows the ordinary sequential path without creating a parallel
region. A call made inside an existing OpenMP parallel region should use a
serial inner path initially. This does not detect independent arcs implemented
with `std::thread`: callers combining both levels must explicitly budget outer
workers and inner threads. Default inner count one avoids accidental nesting;
this feature must not introduce a dependency on the independent-arc scheduler.

With an optional OpenMP backend, an explicit request for multiple threads in a
build without that backend raises a clear configuration error before
propagation. Do not pretend that the request succeeded or silently change the
scientific model.

## Build integration

Introduce a documented optional build flag, disabled by default, and discover
the compiler's OpenMP C++ target through CMake. Propagate its compile/link
requirements to targets instantiating the relevant header code and to installed
Tudat consumers. Export the matching feature information. Keep normal builds
free of the dependency when the option is OFF.

The acceptance build must compile the native libraries and actual Python
bindings consistently from the updated source. Standalone replacement objects
and reused libraries are useful for feasibility measurements, but cannot serve
as acceptance evidence once cache/model layouts change. Do not mix objects from
incompatible definitions of those classes.

## Acceptance gates

- Test validated defaults and object-local policy transfer through C++ and the
  existing Python acceleration factory. Check that separate models do not
  inherit each other's settings or share scratch buffers.
- Compare Legendre values/derivatives against the unchanged serial route for
  full/truncated order, normalized/unnormalized coefficients, cache reuse,
  degree/order changes and near-pole inputs. Test both error and retry paths.
- Compare force values and all saved states/STM/sensitivities for representative
  dynamics and variational propagations. Preserve the current tolerances;
  report observed equality without promising cross-platform bitwise identity.
- Test saved-per-term output, serial fallbacks, unsupported-backend errors,
  scientific-setting serialization defaults and nested/concurrent-model cases.
- Build OFF and ON configurations from consistent sources; run native gravity,
  Legendre and integration tests plus the full offline Python suite. Record any
  external resource blocker explicitly rather than adding ad hoc exclusions.
- Remeasure setup and complete integration with 1/2/4/8 threads, rotating order
  and reporting medians/ranges. Confirm numerical results for every run and
  assess regressions on small workloads. Do not make a fixed speedup a portable
  test assertion.

The contract and experiments are AI-assisted. Maintainer review should focus on
cache ownership, recurrence dependencies, failure recovery, archive semantics
and build/export consistency before exposing the production API.
