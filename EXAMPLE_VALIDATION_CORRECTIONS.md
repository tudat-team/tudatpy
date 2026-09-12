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

## Completed MPC library corrections

- [Legacy reader](src/tudatpy/data/mpc/_legacy.py), [regression test](tests/test_tudatpy/test_data_legacy_calls.py): restore the pre-refactor Note 2 filter inside the removable compatibility module, retain space records for inspection, and resolve identifiers even when filtering is disabled. Commit [44bab9d3d](https://github.com/tudat-team/tudatpy/commit/44bab9d3d1dcde47e063328e1604e442833c2a97). All 51 legacy-call and compatibility-boundary tests passed. The original basic-MPC ingestion prefix runs unchanged, including its WISE preview. Evidence: `.validation/examples-pr157/legacy-reader-regressions.log` and `compare-mpc-readers-fixed.log`. Full original-example compatibility remains pending.
- [Modern MPC reader](src/tudatpy/data_input/tracking_data/mpc/mpc.py), [regression test](tests/test_tudatpy/test_mpc_query.py): restore observatory names, counts, catalog/space filtering, empty-batch behavior and unknown-code handling. Commit [5b43b9540](https://github.com/tudat-team/tudatpy/commit/5b43b954048a81e725ac2aa0914c92f09feaff59). Included in the 31 passing focused MPC/plotting tests; related example runs and figures reviewed.

## Completed MPC estimation correction

[estimation/estimation_with_mpc.py](examples/tudatpy/estimation/estimation_with_mpc.py): explicitly enable weights and catalogue debiasing on ingestion **and** `apply_corrections=True` during collection construction. Correct stale space-observation and body-creation prose. Commit [93770dc](https://github.com/tudat-team/tudatpy-examples/commit/93770dca6a1ed2f59b519a2165519ffbef92bb5e), pushed to examples PR #157.

The complete rerun passed in 44.4 seconds on the spare single CPU; all six figures were visually reviewed and agree with the stored reference trends. Final position difference from SPICE: 35.23 km. Original/new ingestion agrees within 7.1e-14 radians and 2.1 microseconds, with identical weights and link ordering. No modern deprecation warnings; expected legacy warnings were observed. Evidence: `.validation/examples-pr157/mpc-ingestion-comparison.json`, `mpc-unapplied-corrections-diagnosis.json`, and the full `pr1000-modern/estimation__estimation_with_mpc/` run artifacts.

## Completed MPC retrieval plots

[estimation/retrieving_mpc_observation_data.py](examples/tudatpy/estimation/retrieving_mpc_observation_data.py): correct wrapped right ascension, UTC date axes, per-object legends and label spacing; clarify optional satellite environment setup and the UTC numeric filter convention. Commit [fe295ac](https://github.com/tudat-team/tudatpy-examples/commit/fe295ac20f9b319afa396d4e0439e49dd5e67f35), pushed to PR #157. Full final run passed in 70.5 seconds on the spare CPU with all four figures visually reviewed and no deprecations. The final prose-only UTC correction was checked to leave the executable AST unchanged (`source-revalidation.json`). The optional TESS block lacked its external kernel and is not claimed validated.

## Completed observation-wrapper compatibility correction

[observations_setup/__init__.py](src/tudatpy/estimation/observations_setup/__init__.py), [regression test](tests/test_tudatpy/test_data_deprecation_compatibility.py): restore lazy parent access to `observations_setup.observations_wrapper`; legacy function lookup retains its warning and user call-site attribution. Commit [e5b0de89b](https://github.com/tudat-team/tudatpy/commit/e5b0de89b76dcd41ae5ea8e779ced7a2b80b1b10), PR #1000. All 163 deprecation-compatibility tests passed (`.validation/examples-pr157/observations-wrapper-regressions.log`). The original MEX script now passes the missing-attribute failure and proceeds with the expected IFMS warning. Full residual comparison subsequently passed with identical arrays; see the MEX entry below. The original failure is preserved in `pr1000-legacy/estimation__mex_open_loop_residuals/missing-wrapper-failure.log`.

Removal instructions for the wrapper bridge are in [the compatibility README](src/tudatpy/data/README.md), added in commit [1c5430188](https://github.com/tudat-team/tudatpy/commit/1c5430188).

## Completed MEX example correction and comparison

[estimation/mex_open_loop_residuals.py](examples/tudatpy/estimation/mex_open_loop_residuals.py): use named X-band values for IFMS ingestion. Commit [efb9162](https://github.com/tudat-team/tudatpy-examples/commit/efb916297a2373ea86967681c3357de70551e21b), pushed to PR #157. Full modern run: 124.6 seconds; original-script run on the spare CPU: 289.7 seconds. Both figures visually reviewed, all 66,118 residuals and 25,969 subset residuals exactly identical between the two implementations; RMS 19.25/7.12 mHz. The stored-reference difference is not a migration regression.

Modern warnings: none. Original warnings: 28 DeprecationWarnings (including two pre-existing invalid-escape warnings) and one deprecated-interface UserWarning; native reference-point overwrite messages are retained by the old script. Evidence: `.validation/examples-pr157/mex-legacy-modern-comparison.json` and both phases' `estimation__mex_open_loop_residuals/` artifacts, including copied CSV outputs in `runtime-files/mex_output/`.

## Implemented corrections awaiting remaining verification

| File | Correction | Remaining work |
| --- | --- | --- |
| [estimation/load_pds_files.py](examples/tudatpy/estimation/load_pds_files.py) | Preserve all cached files on a day instead of repeating the first; add archive-listing timeouts and HTTP status checks. | Focused duplicate-file check passed; heavy callers pending. |
| [pygmo/asteroid_orbit_optimization/aoo_optimization.py](examples/tudatpy/pygmo/asteroid_orbit_optimization/aoo_optimization.py) | Fix overlapping generation colorbar labels and final orbit-panel titles. | Full optimization passed and figures inspected; final layout rerun queued. |

## Exclusions and open validation issues

- `estimation/improved_estimation_with_mpc.py` is excluded at the user's request and flagged for a separate correction. Its satellite option is ignored, space observations are dropped, and the weight plot is skipped. It has not been accepted or corrected here.
- Space-Track live validation needs locally configured credentials and has not run.
- The initial WISE-preview compatibility failure is corrected by the legacy-filter commit above. The original ingestion prefix now passes without modifying the example; this does not replace its pending full execution.
- GRAIL/MRO full runs, original-example compatibility, and PR #905 validation are still in progress or pending; see the progress log.
