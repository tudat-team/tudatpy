# PR 905 session output index

This index identifies the outputs created during the current PR 905 validation.

## Progress and findings

- Validation state and completed-file list:
  [`PR905_VALIDATION_PROGRESS.md`](PR905_VALIDATION_PROGRESS.md)
- Findings that may also apply to PR 1000:
  [`PR905_POTENTIAL_PR1000_ISSUES.md`](PR905_POTENTIAL_PR1000_ISSUES.md)

## Build output

- Initial configuration: `.validation/pr905-configure.log`
- Initial kernel build: `.validation/pr905-kernel-build.log`
- Metadata conversion rebuild: `.validation/pr905-kernel-rebuild-metadata.log`
- Empty-ancillary-settings rebuild:
  `.validation/pr905-kernel-rebuild-empty-ancillary.log`
- Final current kernel rebuild, including the dataset operations used by the
  migrated examples: `.validation/pr905-kernel-rebuild-modern-example-api.log`
- Final post-cleanup kernel build:
  `.validation/pr905-final-kernel-build.log`

## Tests and examples

- Focused MPC dataset tests:
  `.validation/pr905-mpc-dataset-test.log` and
  `.validation/pr905-mpc-dataset-test-v2.log`
- Observation-vector API tests:
  `.validation/pr905-observation-vector-data-tests.log` and
  `.validation/pr905-observation-vector-data-legacy-tests.log`
- Combined modern observation API and tracking-residual tests:
  `.validation/pr905-observation-modern-api-and-residual-tests.log`
- Final five-file focused suite, including modern APIs, legacy warnings, MPC,
  and tracking residuals: `.validation/pr905-final-focused-python.log`
- Complete final 42-file Python suite: 846 passed, 14 skipped, and 1
  expected-pass in `.validation/pr905-final-full-python-normal-warnings.log`.
  The companion `.validation/pr905-final-full-python.log` records why a global
  warning-as-error policy is unsuitable for tests that deliberately invoke
  deprecated compatibility interfaces.
- Final post-cleanup Python suite: 850 passed, 14 skipped, and 1 expected-pass
  in `.validation/pr905-final-after-test-cleanup.log`. The earlier
  `.validation/pr905-final-eda937bd-python-suite.log` records the same passing
  suite before incidental PR 1000 test warnings were removed.
- Strict PR 1000 test group, with deprecation and pending-deprecation warnings
  treated as errors: 474 passed and 8 skipped in
  `.validation/pr1000-python-tests-no-deprecations-final.log`.
- Executable-syntax comparison for all 13 selectively cleaned examples:
  `.validation/pr905-example-executable-equivalence.txt`
- Final tracking-residual test file and timing:
  `.validation/pr905-tracking-residuals-final.log` and
  `.validation/pr905-tracking-residuals-final-timing.txt`
- Focused bounded-TNF result and timing:
  `.validation/pr905-tracking-residuals-tnf-fixture.log` and
  `.validation/pr905-tracking-residuals-tnf-fixture-timing.txt`
- Example terminal output: `.validation/pr905-examples-*.log`
- Per-example results, numerical records, and figures:
  `.validation/pr905-examples/`
- The MRO range run is in
  `.validation/pr905-examples/estimation__mro_range_estimation_v3/`.
- The MPC data example run is in
  `.validation/pr905-examples/estimation__retrieving_mpc_observation_data_v2/`.
- The MPC estimation example run is in
  `.validation/pr905-examples/estimation__estimation_with_mpc_final/`. Its
  interface finding is recorded in `PR905_POTENTIAL_PR1000_ISSUES.md`.
- Final short estimation-example results use directories ending in `_final`
  beneath `.validation/pr905-examples/`.
- The validated Galilean-moons run is in
  `.validation/pr905-examples/estimation__galilean_moons_state_estimation_pr905/`.
- The validated MEX open-loop run is in
  `.validation/pr905-examples/estimation__mex_open_loop_residuals_final_v2/`.
- The validated GRAIL residual run is in
  `.validation/pr905-examples/estimation__grail_residuals_final/`.
- The validated seven-arc GRAIL estimation run is in
  `.validation/pr905-examples/estimation__grail_odf_estimation_final/`; terminal
  output is in `.validation/pr905-examples-grail-odf-estimation-final.log`.
- Final shorter-example batch output is in
  `.validation/pr905-pending-fast/`; its sequence log is
  `.validation/pr905-pending-fast-batch.log` and its manifest is
  `.validation/pr905-pending-fast.json`.
- The final low-thrust Earth-Mars transfer-window result is in
  `.validation/pr905-long/mission_design__low_thrust_earth_mars_transfer_window/`;
  its PR 1000 numerical and figure comparison is
  `pr1000-comparison.json`, and terminal output is in
  `.validation/pr905-low-thrust-batch.log`.
- The final coupled translational-rotational dynamics result is in
  `.validation/pr905-long/propagation__coupled_translational_rotational_dynamics/`;
  its PR 1000 numerical and figure comparison is
  `pr1000-comparison.json`, and terminal output is in
  `.validation/pr905-coupled-dynamics-batch.log`.
- The final GRAIL SPICE-fit result is in
  `.validation/pr905-long/estimation__grail_spice_fit/`; its PR 1000 numerical
  and figure comparison is `pr1000-comparison.json`, and terminal output is in
  `.validation/pr905-grail-spice-fit-batch.log`.
- Current-master compatibility outputs are stored under
  `.validation/pr905-compat-master-fast-8aaf675a/`,
  `.validation/pr905-compat-master-pygmo-eda937bd/`,
  `.validation/pr905-compat-master-long-eda937bd/`, and
  `.validation/pr905-compat-master-grail-final-eda937bd/`. Each completed
  scientific example contains `arrays.npz`, generated figures,
  `pr1000-comparison.json`, and its complete run log where applicable.
- The durable index for all PR 1000 and PR 905 records is
  [`OBSERVATION_REFACTOR_LOG.MD`](OBSERVATION_REFACTOR_LOG.MD).

New output produced during this validation will be added here as it is created.

The partial MRO TNF estimation output is retained only as an interrupted-run
record. Both MRO TNF examples are excluded from the current analysis and next
steps and must not be counted as validated.
