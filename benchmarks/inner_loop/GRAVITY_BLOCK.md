# Second candidate: acceleration summation and Legendre cache

This experiment isolates two parts of a single Mars propagation: the
spherical-harmonic acceleration sum and the Legendre values/derivatives used by
that sum. It does **not** enable the earlier Hessian experiment. Panel methods,
arc scheduling and Python/free-threading are outside this experiment.

The proposed next-stage API and engineering contract is in
[parallel_spherical_harmonic_gravity.md](../../docs/parallel_spherical_harmonic_gravity.md).
It describes intended behaviour, not an already implemented production API.

## Dependency structure and trial implementation

The Legendre recurrence cannot be split into arbitrary independent degrees.
Sectoral values depend on the preceding diagonal, and vertical values depend
on preceding degrees of the same order. The trial computes the diagonal first,
then evaluates independent order columns in parallel. After that barrier, first
derivatives can be computed by independent rows. A further completed pass is
needed before computing second derivatives. Every slot has one writer;
normalization multipliers and recurrence inputs are read only.

The acceleration trial updates its spherical-harmonic cache before starting
workers. It computes independent gradient terms in parallel, then adds them
in the original degree/order sequence. No regrouping of floating-point sums,
interpolation or model truncation is used to obtain a speedup. Per-term output
maps and the existing near-pole branch fall back to their original serial code.
SPICE and the integrator remain on the calling thread.

Each test compares the original acceleration-factory reference against a fresh
native solver using the same force model, delegates, frames, initial state,
integrator and output recipe. The measured solver uses an experimental subclass
of the existing spherical-harmonic acceleration class to compile its update
method with the trial header. The replacement Legendre object is linked ahead
of the existing libraries. Runtime counters and actual team-size checks verify
that these trial paths really run. No production library or Python binary is
overwritten.

The three variants are acceleration only, Legendre only, and both. The common
baseline uses one thread with both switches off. Hessian calculations stay
sequential in every variant. This prevents attributing the earlier Hessian
speedup to this second candidate.

## Reproduction

The helper assumes an existing compatible Release Unix Makefiles CMake build,
a compiler with OpenMP support, estimation tools and genuine standard SPICE and
JGMRO120D resources. It reuses the configured compiler/include/link settings,
but rebuilds the driver and Legendre implementation in a separate output folder.

```bash
python3 benchmarks/inner_loop/build_gravity_trial.py --build-dir build-python-arcs

# Arguments: degree/order, duration in seconds, mode, repetitions.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE \
  build-python-arcs/benchmarks/gravity_loop/gravity_trial 100 21600 variational 3
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE \
  build-python-arcs/benchmarks/gravity_loop/gravity_trial 100 86400 dynamics 3
```

Each invocation checks 60 normalized/unnormalized Legendre-cache cases against
the serial implementation: degrees 1/20/100, available orders 0/5/100 and five
polynomial arguments from -0.99999 to +0.99999. It compares values, first and
second derivatives, including truncated-order layouts. It then creates an
untimed original-factory propagation as the reference.

A fresh solver is constructed for each timed run. Ten settings (sequential,
plus 2/4/8 threads for each trial variant) are repeated with cyclic order
rotation. Setup and integration are reported separately. Integration includes
cache work, parallel-region synchronization, temporary storage, original-order
term reduction, propagation and normal result processing. Reference generation,
validation, destruction and output are excluded. No speedup assertion is used.

Both modes use a 300-km, nearly circular Mars orbit at inclination 93 degrees,
RKDP87 with fixed 30-second steps and degree/order 100. Variational mode
propagates the state, STM and a sensitivity to Mars mu. Dynamics mode constructs
a genuine dynamics-only simulator; it does not retain the unnecessary
second-derivative cache setup of a variational solver.

## Initial results, 14 November 2026

AMD Ryzen 9 7900X, 12 physical cores / 24 logical CPUs, GCC 14.4. Three
repetitions per setting, with no CPU pinning. Times are integration medians.
All speedups refer to the baseline of the **same invocation**, not to the
previous day's Hessian executable or its different build layout.

Six-hour variational propagation (721 stored epochs):

| Variant | Threads | Median seconds | Range seconds | Speedup |
|---|---:|---:|---|---:|
| Sequential | 1 | 1.407 | 1.404–1.413 | 1.00 |
| Acceleration | 2 | 1.338 | 1.338–1.348 | 1.05 |
| Acceleration | 4 | 1.270 | 1.250–1.277 | 1.11 |
| Acceleration | 8 | 1.455 | 1.425–1.492 | 0.97 |
| Legendre | 2 | 1.338 | 1.333–1.345 | 1.05 |
| Legendre | 4 | 1.247 | 1.239–1.247 | 1.13 |
| Legendre | 8 | 1.265 | 1.247–1.339 | 1.11 |
| Combined | 2 | 1.273 | 1.260–1.280 | 1.11 |
| Combined | 4 | 1.071 | 1.070–1.085 | 1.31 |
| Combined | 8 | 1.185 | 1.161–1.400 | 1.19 |

One-day dynamics-only confirmation (2,881 stored epochs):

| Variant | Threads | Median seconds | Range seconds | Speedup |
|---|---:|---:|---|---:|
| Sequential | 1 | 2.625 | 2.607–2.633 | 1.00 |
| Acceleration | 2 | 2.259 | 2.255–2.273 | 1.16 |
| Acceleration | 4 | 1.981 | 1.977–2.006 | 1.33 |
| Acceleration | 8 | 2.785 | 2.622–3.088 | 0.94 |
| Legendre | 2 | 2.385 | 2.368–2.436 | 1.10 |
| Legendre | 4 | 2.122 | 2.115–2.167 | 1.24 |
| Legendre | 8 | 2.323 | 2.306–2.473 | 1.13 |
| Combined | 2 | 2.015 | 1.996–2.033 | 1.30 |
| Combined | 4 | 1.421 | 1.402–1.431 | 1.85 |
| Combined | 8 | 2.315 | 2.240–2.445 | 1.13 |

The combined four-thread variant reduces complete integration time by about
24% in the variational workload and 46% in the longer dynamics workload. A
six-hour dynamics run also measured a reduction of about 43% (0.662 to 0.375 s).
Eight threads perform worse than four. The trials establish neither linear
scaling nor the precise causes of the eight-thread overhead.

All 90 measured runs (30 per workload) match their original-factory references
with **zero observed difference** at every stored state epoch and, where
applicable, every STM/sensitivity entry. All values must be finite and meet the
unchanged `1e-14 + 1e-11 * abs(reference)` criterion. Short dynamics and
variational smoke runs passed as well. These comparisons establish numerical
parity on these trajectories, not absolute scientific orbit accuracy.

## Scope of the evidence

- Experimental compilation requires both
  `TUDAT_EXPERIMENT_PARALLEL_ACCELERATION` and
  `TUDAT_EXPERIMENT_PARALLEL_LEGENDRE`. No regular CMake configuration enables
  them; default math and public APIs remain unchanged.
- Per-run switches/counters are process-global diagnostics for **one** solver,
  not a production control API. Concurrent outer arc execution, nested teams
  and worker exception translation are not implemented or validated.
- The existing native static libraries were previously built/validated at
  `d011f865f`. The driver and Legendre object are rebuilt against this branch's
  headers. This is a controlled feasibility comparison, not a clean-build
  acceptance test of a production patch. No full Python-suite pass is claimed.
- No full TSan/ASan run, exact-pole sweep, arbitrary force/model audit or
  cross-platform benchmark was performed. The special near-pole and per-term
  output fallbacks need dedicated regression coverage before adoption.
- This experiment does not measure a combined Hessian-plus-gradient-plus-cache
  implementation. The gains from the separate candidates must not be added or
  multiplied to claim an unmeasured speedup.

The measured gain warrants further evaluation for expensive gravity workloads,
particularly dynamics-only propagation. A production design still needs explicit
ownership, exception handling, thread/workload policies, reusable scratch memory
and validation in a fresh library build. Results/code are AI-assisted; review
recurrence dependencies, reference construction and sum order before adoption.

Local metadata, raw JSONL, summaries and build logs are under
`local-profiles/20261114/gravity-block-trial/`. Source, binary and library hashes
and the gravity resource fingerprint are recorded there. Large/raw artifacts
remain outside Git.

## Follow-up profiling

Gated native profiling on 14 November 2026 compares the same one-day dynamics
workload sequentially and with the combined four-thread variant. Only the
integration call is sampled; reference generation, construction and validation
are excluded. A separate repeat with the original uninstrumented executable
confirms medians of 2.627 s and 1.406 s, or 1.87x speedup. All 30 control runs
pass with zero observed differences.

The four-thread cycles profile assigns about 43% of **aggregate user CPU-cycle
weight** to OpenMP runtime self time, predominantly barrier-wait functions.
This is neither 43% of wall time nor an estimate of recoverable speedup.
Separate gated counters show 17.70 CPU seconds versus 7.88 sequential CPU
seconds over three integrations; work completes sooner but uses more total CPU.
Inclusive runtime callchains also contain useful numerical work and must not
be treated as pure overhead.

Generic cache-miss counts increase markedly. Additional cache-miss sampling
resolves C-library vector-store stacks through construction of the temporary
`std::vector<Eigen::Vector3d>` term buffer in the acceleration sum. Repeated
buffer initialization is therefore a concrete memory cost, alongside the
Legendre recurrence/derivative accesses. These architecture-specific events do
not by themselves prove false sharing or quantify a wall-time penalty.

The next controlled candidates are reusable model-owned term storage and fewer
parallel-region transitions, retaining required dependency barriers and original
sum order. No implementation optimization was applied during profiling. All
18 instrumented integrations passed with zero observed differences; no lost
samples were reported. Raw traces, exact commands, summaries and provenance
remain local under `local-profiles/20261114/gravity-profiling/`.

The isolated follow-up [BUFFER_REUSE.md](BUFFER_REUSE.md) measures model-owned
term storage separately; it does not change the number of parallel regions.
