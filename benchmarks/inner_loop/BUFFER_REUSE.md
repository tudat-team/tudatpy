# Isolated term-buffer reuse experiment

14 November 2026. This follow-up changes only the lifetime of the acceleration
term buffer in the [gravity/Legendre trial](GRAVITY_BLOCK.md). Parallel-region
boundaries, scheduling, recurrence formulas and final sum order stay unchanged.
Hessian parallelization remains disabled.

## Implementation and checks

The experimental derived acceleration model owns its own
`std::vector<Eigen::Vector3d>`. It passes that vector to a compile-guarded helper
in the existing acceleration implementation. The buffer is resized only when
the required size differs, and every consumed term is overwritten before the
original-order merge. There is no global or thread-local shared scratch buffer.
The existing base model receives no extra data members or virtual slots.
These experimental hooks are not a production API.

For the control path, the vector retains its original size-constructor on each
force evaluation. In an initial pilot, routing this control through generic
`resize` made it substantially slower. That confounded baseline was rejected.
The final comparison restores the original constructor; the earlier compiled
executable is also rerun separately as a historical control. The pilot data are
retained locally and are not the evidence for the reported gain.

Direct tests exercise grow/shrink sequences, repeated calls and buffers filled
with NaNs before reuse. Results match the fresh-buffer force calculation and
remain finite. The serial path does not resize the scratch. The earlier 60
Legendre-cache test cases also run before each invocation.

For a one-day propagation, counters verify 37,453 temporary vector constructions
in the non-reuse case versus one initial resize and no temporary vectors in the
reuse case. These count C++ buffer operations, not operating-system allocations.

## Reproduction

The compatible existing CMake build and OpenMP compiler requirements are the
same as for the gravity trial. Output goes into a separate build folder, leaving
the earlier executable intact.

```bash
python3 benchmarks/inner_loop/build_buffer_trial.py --build-dir build-python-arcs
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE \
  build-python-arcs/benchmarks/buffer_reuse/buffer_trial 100 86400 dynamics 3
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE \
  build-python-arcs/benchmarks/buffer_reuse/buffer_trial 100 21600 variational 3
```

Seven settings are rotated across repetitions: one-thread sequential, then
2/4/8 threads with temporary or reused buffers. Each measured propagation owns a
fresh model; the first allocation is included in its integration time. Reference
propagation, result validation and model destruction are outside integration
measurements. The underlying model, resources, integrator and validation recipe
are unchanged from the gravity trial.

## Results

AMD Ryzen 9 7900X, degree/order 100, three repetitions per setting. Times are
full integration medians; no profiler is running.

One-day dynamics-only propagation:

| Threads | Buffer | Median seconds | Range seconds | Speedup vs sequential |
|---:|---|---:|---|---:|
| 1 | Sequential | 2.620 | 2.613–2.641 | 1.00 |
| 2 | Temporary | 2.032 | 1.990–2.064 | 1.29 |
| 2 | Reused | 1.808 | 1.806–1.812 | 1.45 |
| 4 | Temporary | 1.410 | 1.409–1.480 | 1.86 |
| 4 | Reused | 1.207 | 1.180–1.222 | 2.17 |
| 8 | Temporary | 1.860 | 1.845–1.870 | 1.41 |
| 8 | Reused | 1.408 | 1.334–1.514 | 1.86 |

Six-hour variational propagation:

| Threads | Buffer | Median seconds | Range seconds | Speedup vs sequential |
|---:|---|---:|---|---:|
| 1 | Sequential | 1.399 | 1.398–1.411 | 1.00 |
| 2 | Temporary | 1.249 | 1.243–1.252 | 1.12 |
| 2 | Reused | 1.198 | 1.195–1.198 | 1.17 |
| 4 | Temporary | 1.072 | 1.069–1.077 | 1.31 |
| 4 | Reused | 1.021 | 1.020–1.042 | 1.37 |
| 8 | Temporary | 1.195 | 1.169–1.392 | 1.17 |
| 8 | Reused | 1.126 | 1.051–1.157 | 1.24 |

At four threads, reuse reduces time by a further 14.4% for dynamics and 4.8% for
variational propagation relative to the same binary's temporary-buffer variant.
The overall dynamics speedup rises from about 1.86x to 2.17x. The separately
rerun pre-change binary measured 2.614 s sequentially and 1.393 s with four
threads, consistent with the restored control. Do not multiply historical
speedup claims to infer the new result: it is measured against its own baseline.

All 42 final measured propagations passed, with zero observed differences at
every stored state epoch and, where applicable, every STM/sensitivity entry.
Shapes, finite values and the existing tolerances remain enforced. Two short
smoke runs and six pre-change control propagations also passed. Four threads
remain the best tested setting; more threads do not imply better performance.

## Limits and provenance

This is still a standalone feasibility test linked to existing native libraries,
not a fresh production-library build. The explicit scratch hook adds no base
object fields, but production adoption still requires uniform build integration,
reviewed ownership/exception handling and broader regression checks. No new
Python API, phase fusion, runtime waiting-policy adjustment, full project-suite
pass or cross-platform guarantee is claimed. The measurements are AI-assisted.

Local output is under `local-profiles/20261114/buffer-reuse-trial/`:
`constructor-dynamics.jsonl`, `constructor-variational.jsonl`,
`constructor-metadata.json` and `final-summary.json` are the final comparison.
`dynamics.jsonl` and `variational.jsonl` are the rejected generic-resize pilot.
Source snapshots before the change and for that pilot are preserved there;
the original gravity-trial executable was not overwritten. No result here
quantifies the benefit of reducing parallel-region transitions.
