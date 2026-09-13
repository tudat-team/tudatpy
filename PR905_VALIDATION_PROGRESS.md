# PR 905 validation progress

This log records validation of PR 905 after integrating the PR 1000
compatibility review. PR 1000's completed example-validation record remains in
`EXAMPLE_VALIDATION_PROGRESS.md`.

## Integration

- PR 905 base: `origin/feature/observation_redesign_refactor` at
  `861e77ec1885206e52c97f61213f921dd65c8557`.
- PR 1000 implementation merge: `b0f0f0bdb`.
- Examples revision: `0113d35864519867da66acd0a3ca3a3f2f2e5d63`.
- Merge resolution retains PR 905's `ObservationDataset` representation for
  collection creation, set mutation, and splitting. The PR 1000 Python parent
  import compatibility bridge and TNF link-delay regression coverage were
  migrated to the PR 905 API.

## Validation status

- Kernel configuration: passed with the `tudatpy-dev` Conda environment.
- Kernel build: passed; only the `kernel` target was built.
- Observation vector interface: renamed consistently from
  `FlattenedObservationData` to `ObservationVectorData`; Python dataset access
  is now `observation_vector_data()`. The C++ filenames, factory/accessor names,
  row lookup, validation names, Python bindings, tests, docs, and MPC example
  use the same terminology. The final kernel-only rebuild passed with six jobs.
- Final focused observation and tracking-residual suite: 80 passed in 29.27 s
  against the rebuilt kernel. Deprecation and future warnings were treated as
  errors; the old Python class and accessor are absent.
- Tracking residual tests: 5 passed in 3.73 s with deprecations treated as
  errors. The MRO TNF test now uses a 5.1 MB repository fixture containing the
  required earlier ramp records and bounded Doppler data. Its runtime fell from
  more than 125 s with the full daily file to 1.42 s in the focused run, while
  preserving the 59-sample, link-delay, mean, and RMS checks. Downloads made by
  this test module now use a temporary file and become visible only after a
  complete transfer.
- [propagation/juice_flybys.py](examples/tudatpy/propagation/juice_flybys.py):
  passed in 16.9 s with zero deprecation diagnostics.
- [estimation/load_pds_files.py](examples/tudatpy/estimation/load_pds_files.py):
  passed with zero deprecation diagnostics.
- [estimation/mro_range_estimation.py](examples/tudatpy/estimation/mro_range_estimation.py):
  passed in 51.5 s with zero deprecation diagnostics after conversion to
  `ObservationDataset`. The eight generated figures were visually checked:
  the geometric residual series and final spectral-density plot are complete
  and physically consistent.
- [estimation/retrieving_mpc_observation_data.py](examples/tudatpy/estimation/retrieving_mpc_observation_data.py):
  passed in 12.0 s with zero deprecation diagnostics after conversion to
  `ObservationDataset`. The four figures were visually checked: the RA/DEC
  residual series and object tracks have the expected dates, scales, and
  grouping. The focused MPC dataset test also passed in 8.0 s.
- [estimation/estimation_with_mpc.py](examples/tudatpy/estimation/estimation_with_mpc.py):
  the final converted run completed in 11.9 s with zero deprecation diagnostics and
  six plausible figures. Its scalar residual plotting now uses
  `observation_vector_data()` for matching scalar-aligned times and set IDs;
  this resolved the interface issue without example-side expansion.
- [estimation/covariance_estimated_parameters.py](examples/tudatpy/estimation/covariance_estimated_parameters.py),
  [estimation/covariance_propagation_example.py](examples/tudatpy/estimation/covariance_propagation_example.py),
  [estimation/full_estimation_example.py](examples/tudatpy/estimation/full_estimation_example.py),
  and [estimation/estimation_dynamical_models.py](examples/tudatpy/estimation/estimation_dynamical_models.py):
  passed in 7.1 s, 21.4 s, 9.3 s, and 26.1 s respectively. All ran without
  warnings and their numerical arrays match the PR 1000 example runs exactly.
  The generated figures were visually checked.
- [estimation/galilean_moons_state_estimation.py](examples/tudatpy/estimation/galilean_moons_state_estimation.py):
  passed in 591.2 s without warnings. Its seven captured numerical arrays
  match the PR 1000 run exactly, and both generated figures were visually
  checked.
- [estimation/mex_open_loop_residuals.py](examples/tudatpy/estimation/mex_open_loop_residuals.py):
  passed in 34.3 s without warnings. Its numerical arrays match the PR 1000
  run exactly; the overall and selected-interval residual RMS values remain
  19.25 mHz and 7.12 mHz. Both figures were visually checked.
- [estimation/grail_residuals.py](examples/tudatpy/estimation/grail_residuals.py):
  passed in 43.4 s without Python warnings or deprecation diagnostics. After
  alignment by observation time, all four model variants have exactly the
  same observation times and residual values as the PR 1000 run. Per-set RMS
  and mean values also match exactly after alignment by set time bounds. Both
  figures were visually checked and match the reference plots.
- [estimation/grail_odf_estimation.py](examples/tudatpy/estimation/grail_odf_estimation.py):
  passed in 343.5 s. Across all seven arcs, all 3,513 observation times and
  station names, residuals with respect to SPICE, pre-fit residuals, post-fit
  residuals, and every RSW orbit-difference sample match the PR 1000 run
  exactly. The seven generated figures are byte-for-byte identical to the
  previously validated figures and were visually checked. Existing messages
  for unsupported ODF type 11 and three high least-squares condition numbers
  have the same text and counts as the reference run; no deprecation warning
  was emitted.
- [propagation/two_stage_rocket_ascent.py](examples/tudatpy/propagation/two_stage_rocket_ascent.py),
  [propagation/walker_constellation.py](examples/tudatpy/propagation/walker_constellation.py),
  [propagation/thrust_between_Earth_Moon.py](examples/tudatpy/propagation/thrust_between_Earth_Moon.py),
  [propagation/solar_system_propagation.py](examples/tudatpy/propagation/solar_system_propagation.py),
  [estimation/mission_data_downloader.py](examples/tudatpy/estimation/mission_data_downloader.py),
  and [propagation/thrust_satellite_engine.py](examples/tudatpy/propagation/thrust_satellite_engine.py):
  passed against the final PR 905 kernel in 3.94 s, 5.10 s, 2.60 s,
  3.65 s, 57.74 s, and 50.65 s. All ran without warnings, deprecation
  diagnostics, or native errors. Across the five numerical examples, 52 shared
  arrays match the PR 1000 references exactly and one Walker array agrees
  within numerical tolerance. All 12 figures were visually checked; the
  solar-system, Earth-Moon thrust, and JUICE-thrust images are pixel-identical
  to their references, while the two-stage and Walker figures contain the same
  results with different canvas dimensions.
- [mission_design/low_thrust_earth_mars_transfer_window.py](examples/tudatpy/mission_design/low_thrust_earth_mars_transfer_window.py):
  passed against the final PR 905 kernel in 275.1 s without warnings,
  deprecation diagnostics, or native errors. Its departure epochs, arrival
  epochs, and complete delta-v grid are exactly equal to the PR 1000 results.
  All five figures were visually checked; the first three are pixel-identical,
  and the final two show the same curves and trajectories with a one-to-two
  pixel difference in canvas width.
- [propagation/coupled_translational_rotational_dynamics.py](examples/tudatpy/propagation/coupled_translational_rotational_dynamics.py):
  passed against the final PR 905 kernel in 936.8 s without warnings,
  deprecation diagnostics, or native errors. All 13 saved orbital, rotational,
  libration, and spectral arrays match the PR 1000 results exactly. All 14
  generated figures are pixel-identical to the reference figures and were
  visually checked.
- [estimation/grail_spice_fit.py](examples/tudatpy/estimation/grail_spice_fit.py):
  passed against the final PR 905 kernel in 4780.9 s using the existing local
  GRAIL archive, without warnings, deprecation diagnostics, or native errors.
  The complete 2,881-by-7 RSW difference history matches the PR 1000 result
  exactly. All five generated figures are pixel-identical to the PR 1000
  reference and were visually checked.
- The modern dataset API now supplies the concise operations these examples
  need directly: precise `observation_time_bounds`, observable-specific link
  definitions, selected dependent-variable addition, reference-point switch
  histories, scalar-aligned link-definition IDs, and a dataset overload of
  `compute_residuals_and_dependent_variables`. These avoid manual array
  expansion and keep the example code direct.
- The remaining PyGMO examples passed against the final kernel without
  warnings: Himmelblau minimization, the three asteroid-orbit optimization
  scripts, Cassini MGA optimization, and hodographic-shaping MGA optimization.
  Their saved numerical arrays match the PR 1000 results exactly, and all
  generated figures were visually checked.
- The requested PR 905 local Python suite and runnable `.py` examples are
  complete, subject to the explicit exclusions below. The SpaceTrack example
  remains blocked by unavailable account credentials and cached response data;
  this external requirement does not indicate a kernel or example failure.
- [estimation/improved_estimation_with_mpc.py](examples/tudatpy/estimation/improved_estimation_with_mpc.py):
  excluded from this validation and migration pass, as requested; it is already
  flagged for correction.
- [estimation/mro_tnf_estimation.py](examples/tudatpy/estimation/mro_tnf_estimation.py)
  and [estimation/mro_tnf_residuals_analysis.py](examples/tudatpy/estimation/mro_tnf_residuals_analysis.py):
  excluded from the current analysis and all remaining validation steps by the
  user's latest instruction. A partial MRO estimation run was terminated and
  is not evidence for any conclusion in this report.

## Selective example cleanup

- Nonfunctional formatting and prose edits were removed from all 13 migrated
  examples. Their module prose and every unchanged comment now match the
  examples branch exactly. Comments removed with deleted legacy operations are
  the only comment differences.
- The executable Python syntax of all 13 cleaned files was compared with the
  exact saved sources used for the completed PR 905 runs. All comparisons are
  identical after excluding standalone Markdown strings, so the numerical and
  plot evidence from those runs remains valid without repeating the runs.
- The final rebuilt kernel is
  `build-pr905/src/tudatpy/kernel.so`, SHA-256
  `eda937bd9f051acdcd8a4cfe9b4199697d922cf0db973e8ee5a519fb7d39814c`.
- The complete 42-file Python suite passed against this kernel: 846 passed,
  14 skipped, and 1 expected-pass in 75.22 s. The 88 reported warnings come
  from tests that intentionally invoke compatibility interfaces and include
  the expected deprecation guidance. The focused modern observation suite
  separately passed 80 tests with deprecation and future warnings treated as
  errors.

## Build isolation

- This branch now uses the primary `tudatpy` worktree and its `build-pr905`
  directory.
- The earlier temporary build was terminated before this primary-worktree build
  was started.

## PR 1000 example compatibility with the PR 905 kernel

- The initial 27-file group completed with 26 direct passes. The untouched
  `mex_open_loop_residuals.py` exposed a legacy collection bug: calling
  `remove_empty_observation_sets()` with no empty sets detached the collection
  from its live `ObservationDataset`, so later residual updates were lost.
- `ObservationCollection::removeSingleObservationSets` now leaves the live
  dataset connection unchanged when no set is removed. Fixed and time-varying
  reference-point updates likewise preserve that connection, and link-end
  changes are propagated to observations saved by an earlier filter.
- Two focused compatibility tests cover the live dataset and saved-filtered-set
  behavior. Both pass against kernel SHA-256
  `aec0e713d42a9166181026a5edbeb1a3173251e7242ec7a05f6370b7c67a7b0c`.
- The untouched PR 1000 MEX example then passed in 60.1 s. Its 66,118 complete
  residuals and 4,643 FDETS residuals are exactly equal to the validated PR 1000
  reference arrays, including RMS values of 0.019248211606 Hz and
  0.008262467295 Hz. Both generated figures are pixel-identical. It emitted
  263 expected `DeprecationWarning` instances from its repeated use of legacy
  observation interfaces, with no native errors or unrelated warnings.
- The five shorter PyGMO examples and the longer hodographic-shaping example
  all passed against the PR 905 kernel without warnings or native errors. All
  49 captured arrays match their PR 1000 references exactly. Thirteen plots are
  pixel-identical; the two Cassini plots and two hodographic plots use the older
  scripts' canvas dimensions but show the same results and were visually
  checked.

## Current-master example compatibility with the PR 905 kernel

- The 27 shorter runnable examples were checked at current-master commit
  `723a2c9a4a2987f8b55d637dd54224a3ca582e49`. Twenty-four run unchanged. Three
  have pre-existing standalone-script defects unrelated to the observation
  interfaces: `mex_open_loop_residuals.py` contains an absolute developer path,
  while `mission_data_downloader.py` and `earth_mars_transfer_window.py` omit
  their `matplotlib.pyplot` import. Temporary validation copies containing only
  those functional corrections passed; the tracked current-master files were
  not edited.
- The corrected MEX run exposed and verified fixes for grouped legacy
  collections. Fixed and ephemeris reference-point changes now retain or adopt
  a live dataset, and collection membership changes rebuild that dataset while
  preserving existing observation-set wrappers and saved filtered sets. Its
  five captured arrays and both figures exactly match the PR 1000 results.
- The five shorter PyGMO examples pass unchanged. All 45 captured arrays and all
  14 figures exactly match the PR 1000 compatibility results.
- The shorter-example figures that differed at the file level were inspected
  side by side. The covariance plots have the same envelopes and observation
  windows with labels changed between the example revisions; the full
  estimation, MPC retrieval, reentry, and two-stage ascent plots show the same
  numerical curves. The KOSMOS map differs only in presentation added by the
  newer example. All figures are complete and visually consistent with their
  PR 1000 counterparts.
- The shorter current-master runs emitted the expected compatibility guidance
  for old observation collection, simulation, weighting, data, Horizons, and
  MPC imports. Each observation warning identifies a supported replacement;
  no tested old observation call failed silently. The old sources also produce
  unrelated Python invalid-escape warnings and 34,567 pending warnings from
  repeated NumPy matrix use in the covariance example. Those originate in the
  untouched example code and are not observation-refactor diagnostics.
- The current-master hodographic-shaping script has a pre-existing process
  startup defect: unlike the PR 1000 version, it creates and evolves its PyGMO
  island at module scope without a `__main__` guard. Worker processes therefore
  re-enter the script and fail while starting further workers. The failed run
  was terminated after this was confirmed. PR 157 already contains the required
  guard and its corresponding PR 1000-compatible script passed in 649 s with
  exact numerical results. This is recorded as an old-example defect and does
  not delay the remaining compatibility checks.
- The untouched current-master Galilean moon estimation passed in 797.0 s. All
  seven captured arrays and both figures are exactly equal to the PR 1000
  compatibility result. It emitted the expected old-import warning and three
  observation-interface deprecation warnings, each directing users to the
  replacement API; no native errors were found.
- The untouched current-master low-thrust Earth-Mars transfer window passed in
  424.8 s. Its departure grid, arrival grid, and complete delta-v surface are
  exactly equal to the PR 1000 result, as are all five figures. It has no old
  observation calls; its warnings are pre-existing string-escape, plotting,
  file-handle, and old SPICE-import diagnostics. No native errors were found.
- The untouched current-master coupled translational-rotational dynamics
  example passed in 935.4 s. All 13 captured orbital, rotational, libration,
  and spectral arrays and all 14 figures are exactly equal to the PR 1000
  result. Its warnings are limited to old string escapes and the old SPICE
  import; no native errors were found.
- The untouched current-master GRAIL residual analysis passed in 51.9 s using
  the reused 63-file cache. All six captured arrays, including the full 11,481
  residual samples, and both figures are exactly equal to the PR 1000 result.
  Its main and four worker logs contain 110 expected compatibility diagnostics;
  the observation warnings identify the replacement dataset calls and migration
  guide. No native errors were found.
- The untouched current-master GRAIL ODF script completed all seven arcs and
  generated seven figures in 692.9 s, but the validation is flagged: arc 2
  reaches 2012-04-09 00:02:09 TDB before the available GRAIL-A CK coverage and
  SPICE reports `NOFRAMECONNECT` twice during its two propagations. The script
  continues and returns success itself, but the saved worker log correctly
  causes the stricter result check to fail. This is an old example/data-coverage
  issue, not an observation-interface exception. It is recorded without a
  source change so it does not delay the remaining work. Despite those logged
  SPICE diagnostics, all six captured result arrays, including the complete
  1,887-by-7 RSW history, and all seven figures are exactly equal to the PR 1000
  result.
- The untouched current-master GRAIL SPICE fit passed in 7,051.7 s using the
  reused cache. Its complete 2,881-by-7 RSW history and all five figures are
  exactly equal to the PR 1000 result. Its three warnings identify the current
  GRAIL mass reader and SPICE import; no native errors were found.

## Final Python validation

- The complete Python suite passed after the final test cleanup against kernel SHA-256
  `eda937bd9f051acdcd8a4cfe9b4199697d922cf0db973e8ee5a519fb7d39814c`:
  850 passed, 14 skipped, and 1 expected-pass in 48.77 s. Only three warnings
  remain: the expected old SPICE import warning and two compatibility warnings
  in the separate inter-arc constraint test, each with migration guidance.
- The 15 PR 1000 test files were then run with `DeprecationWarning` and
  `PendingDeprecationWarning` promoted to errors. Incidental data and optical
  tests now use `ObservationDataset` directly. Legacy-call and legacy
  serialization tests explicitly assert their public warning boundaries and
  suppress only nested expected warnings while checking returned values. The
  strict result is 474 passed and 8 skipped in 28.63 s, with no deprecation
  warning escaping a test.
- The catalogue-sign experiment continues to use its two frozen Eros series,
  high-accuracy Earth setup, Horizons ephemerides, and full numerical
  assertions. Only its observation container changed to the modern dataset;
  the test copies the read-only residual view before applying its local angular
  wrapping transformation.
