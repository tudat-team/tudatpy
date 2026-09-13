# Observation refactor record index

This file indexes the persistent review, implementation, and validation records
for PR 1000 and PR 905. All paths are relative to the current worktree:
`/home/dominic/Tudat/tudat-monorepo/tudatpy`.

## Current authoritative records

### PR 1000

- `EXAMPLE_VALIDATION_PROGRESS.md` — main execution record for the PR 1000
  kernel and examples. It lists the validation sequence, completed files,
  warning checks, plot checks, exclusions, and remaining work at the end of that
  stage.
- `EXAMPLE_VALIDATION_CORRECTIONS.md` — correction ledger for PR 1000. It lists
  changed files, the reason for each correction, validation evidence, branch and
  commit locations, and publication state.
- `examples/tudatpy/EXAMPLE_MIGRATION_STATUS.md` — detailed per-example record
  created during the PR 1000-compatible examples migration. It includes
  runtimes, numerical observations, figure comparisons, input-data locations,
  and blocked examples. References to notebook synchronization describe the
  earlier work; the current PR 905 work changes and validates `.py` files only.
- `claude_review_followup_report.md` — follow-up report for the AI review of the
  data-refactor branch. It records the reviewed revision, issue groups,
  implemented changes, and supporting diff.

### PR 905

- `PR905_VALIDATION_PROGRESS.md` — current authoritative PR 905 progress log.
  It records the PR 1000 merge, kernel build, test results, validated `.py`
  examples, numerical and visual comparisons, exclusions, and the proof used to
  preserve completed run results during the selective source cleanup.
- `PR905_SESSION_OUTPUT_INDEX.md` — index of detailed PR 905 build logs, test
  logs, example results, saved arrays, and generated figures under
  `.validation/`.
- `PR905_POTENTIAL_PR1000_ISSUES.md` — findings made during PR 905 work that may
  also apply to PR 1000. These entries are recorded for later review and are not
  acted on during the PR 905 stage unless separately required.
- `PR905_REVIEW.md` — earlier high-level PR 905 review notes, reviewed revision,
  issue assessment, and CI status.

## Supporting and historical records

### Observation refactor and PR 1000

- `codex_handoff_observation_dataset_2026-07-01.md` — early implementation
  handoff for the observation-dataset work. It is useful for design history but
  predates the current names and validation results.
- `python_script_migration_handoff_2026-07-01.md` — early handoff for Python
  script migration. Later validation logs supersede its status information.
- `observation_dataset_redesign_pr_text.md` — draft PR text describing the
  observation-dataset redesign.
- `created_branches.md` — branch inventory created during the earlier review and
  migration work.
- `pr_messages.md` — saved PR titles and message drafts from the earlier branch
  reorganization.

### PR 905 architecture review

- `.validation/pr905-architecture-review-20260907/review-report.md` — detailed
  architecture review, issue findings, implemented corrections, and CI state.
- `.validation/pr905-architecture-review-20260907/verification.md` — commands
  and results for the focused C++ and Python checks used during that review.
- `.validation/pr905-architecture-review-20260907/architecture-notes.md` —
  design notes gathered while reviewing ownership, ordering, mutation, and
  compatibility behavior.
- `.validation/pr905-architecture-review-20260907/workflow-evidence.md` —
  evidence for end-to-end estimation and covariance workflows.
- `.validation/pr905-architecture-review-20260907/continuation-state.md` —
  recovery record for the earlier long-running review session. Its status is
  historical; `PR905_VALIDATION_PROGRESS.md` contains the current state.
- `.validation/examples-pr157/pr905-api-preparation.md` — initial PR 905 API
  inventory used before migrating examples. Later source and validation logs
  supersede preliminary conclusions in this file.

## Detailed machine-readable and visual evidence

- `.validation/pr905-final-kernel-build.log` — final PR 905 kernel-only build
  output using six compile jobs.
- `.validation/pr905-final-focused-python.log` — final focused Python result:
  80 tests passed in 29.27 seconds with deprecation and future warnings treated
  as errors.
- `.validation/pr905-examples/` — per-example PR 905 result files, numerical
  arrays, figures, and complete run logs for the already validated estimation
  and mission examples.
- `.validation/pr905-pending-fast/` — results and figures for the final group of
  shorter PR 905 examples.
- `.validation/examples-pr157/` — PR 1000 modern and legacy example results,
  comparison data, figures, and scripts used to record the validation.
- `.validation/pr905-compat-master-*/` — current-master example runs against
  the final PR 905 kernel, including saved arrays, figure comparisons, warning
  records, and the complete long GRAIL outputs.
- `.validation/pr905-final-after-test-cleanup.log` — final complete Python suite
  result: 850 passed, 14 skipped, and 1 expected-pass, with three expected
  warnings outside the PR 1000 test group.
- `.validation/pr1000-python-tests-no-deprecations-final.log` — strict PR 1000
  test result: 474 passed and 8 skipped with deprecation warnings treated as
  errors.
- `.validation/pr905-examples-before-scope-cleanup.diff` — exact saved diff from
  before nonfunctional example formatting and prose edits were removed.
- `.validation/pr905-tested-examples-path.txt` — location of the saved exact
  example sources used to prove executable equivalence after cleanup.
- `tests/test_tudatpy/fixtures/mro_short_arc_2012_077.tnf` and
  `tests/test_tudatpy/fixtures/README_mro_short_arc.md` — bounded MRO TNF test
  input and its provenance, used to keep the tracking-residual test well below
  the two-minute target.

## Current execution policy and exclusions

- Run at most six threads in total; the final kernel build uses at most six.
- Change and validate `.py` examples only. Do not modify notebooks.
- Preserve simple example code and avoid manual scalar expansion such as
  `np.repeat`.
- Commit completed PR 905 source and examples work once per repository, then
  perform the PR 1000 and older-example backward-compatibility runs.
- Record issues more difficult than moderately easy in
  `PR905_VALIDATION_PROGRESS.md` and continue without delaying the overall run.
- `improved_estimation_with_mpc.py`, `mro_tnf_estimation.py`, and
  `mro_tnf_residuals_analysis.py` are excluded from the remaining analysis as
  requested.
