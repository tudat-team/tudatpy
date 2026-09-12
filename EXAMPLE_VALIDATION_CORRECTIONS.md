# Example validation corrections

This index identifies every correction from the PR #1000 / examples PR #157 validation, its verification, and where to find its commit. Full execution and plot-review status is in [EXAMPLE_VALIDATION_PROGRESS.md](EXAMPLE_VALIDATION_PROGRESS.md). The user authorized incremental commits and pushes on 2026-09-12. Unfinished validation is explicitly identified below.

## Completed example corrections

Commit: [28e5337 — Correct validated example execution and scientific plots](https://github.com/tudat-team/tudatpy-examples/commit/28e5337114f02c01b7cbf695c9cc883c4e479f0c), branch `feature/data-refactor` (PR #157).

| File | Correction | Verification |
| --- | --- | --- |
| [estimation/mro_range_estimation.py](examples/tudatpy/estimation/mro_range_estimation.py) | Remove machine-specific kernel import path. | Complete three-mission run, eight plots reviewed; no deprecations. |
| [mission_design/hodographic_shaping_mga_optimization.py](examples/tudatpy/mission_design/hodographic_shaping_mga_optimization.py) | Guard multiprocessing entry points so imported worker modules cannot recursively start islands. | Full 40 generations, population 1,000; both plots reviewed; no deprecations. |
| [mission_design/low_thrust_earth_mars_transfer_window.py](examples/tudatpy/mission_design/low_thrust_earth_mars_transfer_window.py) | Correct meter coordinate labels and remove a second mass division from thrust acceleration. | Complete 200 × 442 grid and both trajectories; all five plots reviewed; no deprecations. |
| [pygmo/asteroid_orbit_optimization/aoo_design_space_exploration.py](examples/tudatpy/pygmo/asteroid_orbit_optimization/aoo_design_space_exploration.py) | Use each varied parameter's constraint colors, correct distance/angular units, and select inclination for the inclination/node response surface. | Full 2,000 Monte Carlo and 2,401 factorial simulations; corrected plotting blocks replayed on saved arrays with coordinate/color assertions; all three plots reviewed. |

Local evidence: `.validation/examples-pr157/pr1000-modern/<path-with-slashes-replaced-by-__>/`, including `result.json`, figures and numerical arrays. The asteroid design-space directory also contains `plot-revalidation.json`.

## Completed library plotting corrections

Files: [porkchop plotting](src/tudatpy/trajectory_design/porkchop/_plot_porkchop.py), [plotting helpers](src/tudatpy/plotting/_helpers.py).

Porkchop colorbars now identify km/s and km²/s², matching the scaled values. A raw plotting-helper docstring prevents an invalid-escape deprecation warning. Full Earth–Mars and low-thrust examples and their plots were reviewed. The combined plotting/MPC focused suite passed 31 tests with deprecations treated as errors (`.validation/examples-pr157/plotting-and-mpc-tests.log`). Commit: [fcfaa048a — Correct porkchop units and plotting docstring escapes](https://github.com/tudat-team/tudatpy/commit/fcfaa048ac6157df35bf704809408df68efb0b2a), branch `codex/pr900-review-by-category-20260912` (PR #1000).

## Implemented corrections awaiting remaining verification

| File | Correction | Remaining work |
| --- | --- | --- |
| [MPC reader](src/tudatpy/data_input/tracking_data/mpc/mpc.py), [MPC tests](tests/test_tudatpy/test_mpc_query.py) | Restore observatory names, counts, catalog filtering and unknown-code handling. | Focused tests passed; complete related example and legacy comparison. |
| [estimation/estimation_with_mpc.py](examples/tudatpy/estimation/estimation_with_mpc.py) | Explicitly enable weighting and catalog debiasing to preserve old conversion defaults. | Six figures reviewed; investigate original-reader comparison before acceptance. |
| [estimation/retrieving_mpc_observation_data.py](examples/tudatpy/estimation/retrieving_mpc_observation_data.py) | Correct wrapped right ascension, UTC time axes, legends, and unsupported satellite-conversion claims. | Four plots reviewed; final label-spacing rerun queued. |
| [estimation/mex_open_loop_residuals.py](examples/tudatpy/estimation/mex_open_loop_residuals.py) | Pass named X-band values to IFMS conversion. | Full run and both figures reviewed; compare numerical difference with original example. |
| [estimation/load_pds_files.py](examples/tudatpy/estimation/load_pds_files.py) | Preserve all cached files on a day instead of repeating the first; add archive-listing timeouts and HTTP status checks. | Focused duplicate-file check passed; heavy callers pending. |
| [pygmo/asteroid_orbit_optimization/aoo_optimization.py](examples/tudatpy/pygmo/asteroid_orbit_optimization/aoo_optimization.py) | Fix overlapping generation colorbar labels and final orbit-panel titles. | Full optimization passed and figures inspected; final layout rerun queued. |

## Exclusions and open validation issues

- `estimation/improved_estimation_with_mpc.py` is excluded at the user's request and flagged for a separate correction. Its satellite option is ignored, space observations are dropped, and the weight plot is skipped. It has not been accepted or corrected here.
- Space-Track live validation needs locally configured credentials and has not run.
- The initial original-MPC ingestion diagnostic fails when the original script indexes an empty WISE-observation selection. This is an open backwards-compatibility finding, not a passed original-example run. Evidence: `.validation/examples-pr157/compare-mpc-readers.log`.
- GRAIL/MRO full runs, original-example compatibility, and PR #905 validation are still in progress or pending; see the progress log.
