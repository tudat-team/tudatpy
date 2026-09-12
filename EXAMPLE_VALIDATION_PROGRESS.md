# Example validation progress

Last updated: 2026-09-12T20:01:27+00:00

**The user authorized incremental commits and pushes on 2026-09-12. Publish each completed, verified part.**

Correction locations, evidence, and published commits: [EXAMPLE_VALIDATION_CORRECTIONS.md](EXAMPLE_VALIDATION_CORRECTIONS.md).

## Validation sequence

1. Build PR #1000; run and review all updated examples from examples PR #157 without deprecation warnings.
2. Commit and push verified kernel/example fixes incrementally; run the original examples against the PR #1000 kernel and verify expected deprecation warnings.
3. Merge the validated PR #1000 branch into PR #905 locally; branch the validated examples and rebuild the kernel.
4. Run the PR #1000-compatible examples unchanged against PR #905, checking compatibility and warnings.
5. Migrate examples minimally to PR #905, rerun/review all examples without deprecation warnings, and commit and push each verified part.

## Current setup

- Final MRO input audit: 328/328 published 2012 TNFs match exact archive sizes (20.94 GB); 0 mismatches.
- MRO download/repair progress: 91/91 files verified against published byte sizes; inputs are installed only after successful verification.
- Validate Python scripts and every generated plot only. Notebooks are left untouched; the two earlier notebook edits were undone.
- `improved_estimation_with_mpc.py` is excluded at user request and flagged for separate correction (ignored satellite option and skipped weight plot).
- Resource limit: one example at a time; hard affinity to eight CPUs; at most seven pool workers plus the parent; numerical libraries use one thread. The full-year MRO run uses two workers to limit memory while retaining all six intervals. All earlier processes were terminated and verified stopped before applying this limit.
- Kernel branch: `codex/pr900-review-by-category-20260912`; starting PR #1000 commit `4fca01f63`.
- Examples branch: `examples/tudatpy`, `feature/data-refactor`; starting PR #157 commit `403f590`.
- Focused MPC/plotting regressions: 31 passed with deprecations treated as errors. Fresh compilation of package and included example sources reports no warnings.
- PR #1000 kernel rebuilt successfully; imported extension verified as `build/src/tudatpy/kernel.so`.
- Original examples for compatibility: examples `origin/master` (PR #157 base). The kernel remains PR #1000 for that pass.
- Existing MRO edit preserved in the examples git stash and `.validation/examples-pr157/preexisting-examples.patch`.
- Indexed 1,409 cached mission input files; reusing the 29 GB MRO cache and existing GRAIL cache.
- MRO input audit: 328 published 2012 TNFs, 90 initially missing (5.57 GB) and one truncated January file. Sequential downloads validate published sizes before atomic replacement; details in mro-prefetch-results.json.
- Space-Track live execution awaits locally configured credentials; no credentials were present initially.

## PR #1000: updated examples

`Execution passed` means the complete configured Python run returned successfully with Python deprecations treated as errors. It is not yet full numerical/figure verification. `Validated` additionally requires result review and a native deprecation-message check. Helper modules are verified through their callers.

| File | Status | Runtime | Figures | Review / issue |
| --- | --- | ---: | ---: | --- |
| [data_retrieval/spacetrack_example.py](examples/tudatpy/data_retrieval/spacetrack_example.py) | Awaiting credentials | — | — | Live account-dependent download has not been run; SPACETRACK_USER/PASS unavailable at start. |
| [estimation/covariance_estimated_parameters.py](examples/tudatpy/estimation/covariance_estimated_parameters.py) | Validated | 7.5 s | 2 | Full run; both covariance figures reviewed against notebook, finite outputs and expected ellipsoids; no deprecations. |
| [estimation/covariance_propagation_example.py](examples/tudatpy/estimation/covariance_propagation_example.py) | Validated | 24.4 s | 2 | Full run; Cartesian/RSW uncertainty figures reviewed against notebook, expected magnitudes and evolution; no deprecations. |
| [estimation/estimation_dynamical_models.py](examples/tudatpy/estimation/estimation_dynamical_models.py) | Validated | 33.0 s | 3 | Full run; all three figures reviewed against notebook, residual and state-error trends consistent; no deprecations. |
| [estimation/estimation_with_mpc.py](examples/tudatpy/estimation/estimation_with_mpc.py) | Correction application fixed; rerun pending | 16.8 s | 6 | Found the collection factory default left requested catalogue corrections unapplied. Set apply_corrections=True; diagnostic confirmed old values were debiased and new values raw (bias difference agreement 2.2e-16 rad). Corrected comparison running; full six-figure rerun queued. |
| [estimation/full_estimation_example.py](examples/tudatpy/estimation/full_estimation_example.py) | Validated | 12.4 s | 3 | Full run; all three scientific figures reviewed, converged residuals consistent with reference; the extra old notebook figure is an empty placeholder; no deprecations. |
| [estimation/galilean_moons_state_estimation.py](examples/tudatpy/estimation/galilean_moons_state_estimation.py) | Validated | 606.5 s | 2 | Full 2031–2035 fit; both figures visually match reference, position errors 0–15 km and Laplace-angle differences about 0.001 deg. All numerical outputs finite; no deprecations. |
| [estimation/grail_examples_functions.py](examples/tudatpy/estimation/grail_examples_functions.py) | Helper — callers pending | — | — |  |
| [estimation/grail_odf_estimation.py](examples/tudatpy/estimation/grail_odf_estimation.py) | Pending | — | — |  |
| [estimation/grail_residuals.py](examples/tudatpy/estimation/grail_residuals.py) | Pending | — | — |  |
| [estimation/grail_spice_fit.py](examples/tudatpy/estimation/grail_spice_fit.py) | Running | — | — | 12/20 full-day fits completed; all five dates and four model/parameter setups retained, at most seven workers. |
| [estimation/improved_estimation_with_mpc.py](examples/tudatpy/estimation/improved_estimation_with_mpc.py) | Excluded — correction required | — | — | Excluded at user request. Satellite-data option is ignored, space observations are dropped, and the weight plot is skipped. Requires a separate correction; not validated. |
| [estimation/kosmos482_reentry.py](examples/tudatpy/estimation/kosmos482_reentry.py) | Validated | 12.6 s | 2 | Full supplied-TLE rerun; both ground-track and altitude figures match reference. Correct PROJ data path eliminates native environment error; no deprecations. |
| [estimation/load_pds_files.py](examples/tudatpy/estimation/load_pds_files.py) | Helper fix verified; callers pending | — | — | Fixed duplicate-first-file bug for cached days with multiple TNFs (present on 2012 day 152). Focused check retained both files exactly once without network. Added timeouts/status checks for archive listings; heavy callers pending. |
| [estimation/mex_open_loop_residuals.py](examples/tudatpy/estimation/mex_open_loop_residuals.py) | Plots reviewed; legacy comparison pending | 124.6 s | 2 | Full run passed after X-band argument fix; both residual/histogram figures reviewed, all arrays finite, no deprecations. RMS 19.25 mHz overall / 7.12 mHz subset differs from stored reference (20.31 / 9.93); compare original script against this kernel before final acceptance. |
| [estimation/mission_data_downloader.py](examples/tudatpy/estimation/mission_data_downloader.py) | Validated | 57.4 s | 0 | Full mission downloader run passed using existing archives and required downloads, including archive cleanup. No figures by design and no deprecations. |
| [estimation/mro_range_estimation.py](examples/tudatpy/estimation/mro_range_estimation.py) | Validated | 70.5 s | 8 | All three missions complete; final RMS 5.94/3.20/5.32 m; all eight figures reviewed; no deprecations. Removed hard-coded kernel path. Python script only; notebooks untouched. |
| [estimation/mro_tnf_estimation.py](examples/tudatpy/estimation/mro_tnf_estimation.py) | Pending | — | — |  |
| [estimation/mro_tnf_residuals_analysis.py](examples/tudatpy/estimation/mro_tnf_residuals_analysis.py) | Pending | — | — |  |
| [estimation/mro_utils.py](examples/tudatpy/estimation/mro_utils.py) | Helper — callers pending | — | — |  |
| [estimation/retrieving_mpc_observation_data.py](examples/tudatpy/estimation/retrieving_mpc_observation_data.py) | Plot layout adjustment; rerun pending | 8.7 s | 4 | All four corrected plots reviewed: complete sky coverage, per-object legends, correct wrapped RA and UTC time panels. Increased sky-axis label spacing to avoid tick overlap; final rerun pending. Optional TESS environment block unavailable and not claimed validated; satellite conversion explicitly deferred. |
| [estimation/tudat_azimuth_elevation_example.py](examples/tudatpy/estimation/tudat_azimuth_elevation_example.py) | Validated | 10.7 s | 1 | Full run; all four station panels reviewed. Tudat and Horizons azimuth/elevation points agree visually, with sensible horizon crossings and complete labels; no deprecations. |
| [mission_design/cassini1_mga_optimization.py](examples/tudatpy/mission_design/cassini1_mga_optimization.py) | Validated | 2.3 s | 2 | Full 800 generations; convergence and trajectory figures match reference; finite outputs, no deprecations. |
| [mission_design/earth_mars_transfer_window.py](examples/tudatpy/mission_design/earth_mars_transfer_window.py) | Validated | 17.3 s | 5 | Full grid rerun; all five figures reviewed. Transfer windows and contour minima match reference; delta-v and C3 colorbars now correctly labeled km/s and km²/s². No deprecations. |
| [mission_design/hodographic_shaping_mga_optimization.py](examples/tudatpy/mission_design/hodographic_shaping_mga_optimization.py) | Validated | 572.8 s | 2 | Full 40 generations / 1,000 population; both figures match stored reference, final delta-v 43.546 km/s. Fixed script multiprocessing main guards; no deprecations. |
| [mission_design/low_thrust_earth_mars_transfer_window.py](examples/tudatpy/mission_design/low_thrust_earth_mars_transfer_window.py) | Validated | 817.7 s | 5 | Full 200×442 grid rerun and both complete trajectories; all five plots reviewed with corrected kilometer-based colorbars, meter coordinates, and acceleration scaling. All grid values finite; no deprecations. |
| [mission_design/mga_trajectories.py](examples/tudatpy/mission_design/mga_trajectories.py) | Validated | 1.1 s | 3 | Full run, all three plots reviewed against reference: both transfer geometries and TNW thrust histories agree, flyby markers/units consistent; no deprecations. |
| [propagation/coupled_translational_rotational_dynamics.py](examples/tudatpy/propagation/coupled_translational_rotational_dynamics.py) | Validated | 1223.8 s | 14 | Full 1,224 s execution including all damping iterations; all 14 plots reviewed. Orbital histories match reference, free libration mode is suppressed while orbital harmonics remain, units and labels consistent; finite arrays and no deprecations. |
| [propagation/impact_manifolds_lpo_cr3bp.py](examples/tudatpy/propagation/impact_manifolds_lpo_cr3bp.py) | Validated | 54.5 s | 1 | Full run; both planar projections reviewed, libration orbit and inward/outward manifold branches are visible relative to Phobos, with correct kilometer axes. Numerical arrays finite; no deprecations. |
| [propagation/juice_flybys.py](examples/tudatpy/propagation/juice_flybys.py) | Validated | 17.2 s | 2 | Full run; all 2 figures reviewed against notebook; expected trends and finite numerical outputs; no Python/native deprecations. |
| [propagation/keplerian_satellite_orbit.py](examples/tudatpy/propagation/keplerian_satellite_orbit.py) | Validated | 0.9 s | 1 | Full 86,400 s / 8,641 states; energy variation 3.83e-10, angular momentum 1.90e-10; orbit figure checked; no Python/native deprecations. |
| [propagation/linear_sensitivity_analysis.py](examples/tudatpy/propagation/linear_sensitivity_analysis.py) | Validated | 2.3 s | 2 | Full run; all 2 figures reviewed against notebook; expected trends and finite numerical outputs; no Python/native deprecations. |
| [propagation/panelled_radiation_target.py](examples/tudatpy/propagation/panelled_radiation_target.py) | Validated | 1.3 s | 2 | Full run; all 2 figures reviewed against notebook; expected trends and finite numerical outputs; no Python/native deprecations. |
| [propagation/perturbed_satellite_orbit.py](examples/tudatpy/propagation/perturbed_satellite_orbit.py) | Validated | 2.3 s | 4 | Full run; all 4 figures reviewed against notebook; expected trends and finite numerical outputs; no Python/native deprecations. |
| [propagation/reentry_trajectory.py](examples/tudatpy/propagation/reentry_trajectory.py) | Validated | 1.7 s | 7 | Full run; all 7 figures reviewed against notebook; expected trends and finite numerical outputs; no Python/native deprecations. |
| [propagation/separation_satellites_diff_drag.py](examples/tudatpy/propagation/separation_satellites_diff_drag.py) | Validated | 11.3 s | 4 | Full run; all 4 figures reviewed against notebook; expected trends and finite numerical outputs; no Python/native deprecations. |
| [propagation/solar_system_propagation.py](examples/tudatpy/propagation/solar_system_propagation.py) | Validated | 2.2 s | 2 | Full run; all 2 figures reviewed against notebook; expected trends and finite numerical outputs; no Python/native deprecations. |
| [propagation/thrust_between_Earth_Moon.py](examples/tudatpy/propagation/thrust_between_Earth_Moon.py) | Validated | 1.9 s | 3 | Full run; all 3 figures reviewed against notebook; expected trends and finite numerical outputs; no Python/native deprecations. |
| [propagation/thrust_satellite_engine.py](examples/tudatpy/propagation/thrust_satellite_engine.py) | Validated | 101.5 s | 3 | Full rerun; all three figures rechecked with complete titles/labels. Inclination response, mass consumption, and two thrust windows match reference; no deprecations. |
| [propagation/two_stage_rocket_ascent.py](examples/tudatpy/propagation/two_stage_rocket_ascent.py) | Validated | 1.7 s | 1 | Full run; all 1 figures reviewed against notebook; expected trends and finite numerical outputs; no Python/native deprecations. |
| [propagation/walker_constellation.py](examples/tudatpy/propagation/walker_constellation.py) | Validated | 2.5 s | 3 | Full run; all 3 figures reviewed against notebook; expected trends and finite state outputs; intentional NaN longitude breaks checked; no Python/native deprecations. |
| [pygmo/asteroid_orbit_optimization/aoo_custom_environment.py](examples/tudatpy/pygmo/asteroid_orbit_optimization/aoo_custom_environment.py) | Validated | 1.3 s | 1 | Full run; orbit/asteroid geometry matches reference, all eight numeric arrays finite, meter axes complete; no deprecations. |
| [pygmo/asteroid_orbit_optimization/aoo_design_space_exploration.py](examples/tudatpy/pygmo/asteroid_orbit_optimization/aoo_design_space_exploration.py) | Validated | 248.2 s | 3 | Full 2,000 Monte Carlo and 2,401 factorial-design simulations passed. All three figures reviewed after plotting-only replay of saved arrays: corrected distance/angular units, per-variable constraint colours, and inclination/node coordinates; no deprecations. Plot/data assertions passed. |
| [pygmo/asteroid_orbit_optimization/aoo_optimization.py](examples/tudatpy/pygmo/asteroid_orbit_optimization/aoo_optimization.py) | Plot layout fixes; rerun pending | 114.0 s | 5 | Full 25-generation / 48-member optimisation passed with five inspected figures, finite arrays and no deprecations. Pareto front and near-polar final orbits sensible; adjusted overlapping first-generation colorbar labels and final orbit-panel titles. Rerun pending. |
| [pygmo/himmelblau_minimization.py](examples/tudatpy/pygmo/himmelblau_minimization.py) | Validated | 8.3 s | 3 | Full run; all three plots reviewed. Optimiser convergence reaches about 3e-9, contour minima and zoom around (3, 2) agree with the objective; finite arrays and no deprecations. |

Current counts: Awaiting credentials: 1; Correction application fixed; rerun pending: 1; Excluded — correction required: 1; Helper fix verified; callers pending: 1; Helper — callers pending: 2; Pending: 4; Plot layout adjustment; rerun pending: 1; Plot layout fixes; rerun pending: 1; Plots reviewed; legacy comparison pending: 1; Running: 1; Validated: 31.

## Later passes

Original examples / PR #1000: pending completion of the updated-example validation.

PR #905 compatibility and migration: pending the completed PR #1000 passes.

## Evidence

Per-example logs, figures, arrays, warnings and source hashes: `.validation/examples-pr157/pr1000-modern/`.
Input inventory: `.validation/examples-pr157/input-cache-inventory.json`.
Previously extracted reference figures: `.validation/examples-pr157/notebook-source-comparison.json` and `notebook-reference/`.
