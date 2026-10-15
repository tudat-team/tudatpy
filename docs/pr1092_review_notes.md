> **15 October: numerical verification passed.** All 29 native executables, 810 executed full-suite Python cases,
> 14 focused/architecture checks and 209 autosummary names passed. The full offline suite retains 101 documented skips.
> The seven required resource files and their CI installation step are included and checked, without recompilation.
> Dominic explicitly authorized committing and pushing this follow-up; the earlier local-only restriction is superseded.

# PR 1092 review follow-up — implementation status (handover document)

The earlier review stages finished on 15 October 2026. Their commits are listed below. The subsequent five local handover
commits and the current build/test follow-up are in `review/pr1092-external-simplify`. Work in that review worktree;
the older coordinator checkout is not the source for this validation. Dominic subsequently authorized a push after the
repository resource files and CI installation step are added.

## Pushed commits (on top of 226d87de1)
1. `69dff7426` Simplify PR 1092 time handling and restore base compatibility (proposed diff items, S1, compatibility fixes C1/C3 etc.)
2. `b1bf14869` Use the station-history file reader in the binding test; add PR 1092 review and notes
3. `6954ab6e1` Read Desai-Sibois and official ITRF tables from tudat resources (S3; resource files now included below)
4. `1969a1c49` Use one fourth-order central-difference helper for station and Earth orientation rates (S4)
5. `12bddd0c6` Remove per-query station-loading overlap bookkeeping from the generic deformation aggregator (S7 part)
6. `92b9c7b16` Update PR 1092 review notes with final status

Net effect of 1–5 against 226d87de1: native code +575 / −2,046 lines (net −1,471), Python bindings +68 / −179 (net −111), tests +136 / −376 (net −240); the review documents under docs/ are added.

Keep this file current. Another agent must be able to continue from it.

## Locations
- Review worktree: `/home/dominic/Tudat/tudat-monorepo/tudatpy-pr1092-external-review`, local branch `review/pr1092-external-simplify`.
  Review start: `226d87de1`; pre-handover revision: `92b9c7b16`; five-commit handover revision: `70aabd7c6`.
- Most complete build directory: `/home/dominic/Tudat/pr1092-review-build`.
  It uses the review worktree and the original Release/GNU 15.2 options, with serialization now enabled for the full Python suite.
  Use at most eight compiler jobs in total. The 29 selected native targets are recorded in
  `/home/dominic/Tudat/pr1092-review-build/verification-allpass-20261015/selected-native-targets.txt`.
- Review: `/home/dominic/Tudat/tudat-monorepo/REVIEW-DONE.md`. The first proposed diff is `REVIEW-PROPOSED-CHANGES.diff`.
  The list of riskier redesigns is `REVIEW-RISKIER-REDESIGNS.md`.
- Resource additions are now included under `resources/` in the review worktree. Their source copies are in
  `/home/dominic/Tudat/tudat-monorepo/tudat-resources-additions/`; identical files were used in the completed validation.

## Repository resources and CI installation
Dominic requested that the new mandatory tables be included in this repository so CI does not depend on a separate
resource publication. The following files are copied unchanged under `resources/`:
- `earth_orientation/desaiSibois2016PolarMotionOceanTidesAmplitudes.txt`
- `earth_orientation/desaiSibois2016PolarMotionOceanTidesFundamentalArgumentMultipliers.txt`
- `earth_orientation/desaiSibois2016Ut1OceanTidesAmplitudes.txt`
- `earth_orientation/desaiSibois2016Ut1OceanTidesFundamentalArgumentMultipliers.txt`
- `earth_orientation/README_desaiSibois2016.txt`
- `station_locations/Transfo-ITRF2014_ITRFs.txt` (SHA-256 6cc1c068…0543, the official IGN file)
- `station_locations/Transfo-ITRF2020_TRFs.txt` (SHA-256 0d8dc872…1a68, the official IGN file)

CI copies this directory into the native resource directory returned by `tudatpy.data_input.resource_paths.get_resource_path`
after package installation and before Python/native tests. It merges these additions with the installed standard resources.
The existing C++ readers and default paths stay unchanged. For a source checkout outside CI, run the same operation from
the repository root with the built/installed tudatpy package on Python's path:

```python
from pathlib import Path
import shutil
from tudatpy.data_input.resource_paths import get_resource_path

shutil.copytree(Path("resources"), Path(get_resource_path()), dirs_exist_ok=True)
```

## Decisions from Dominic (14 October 2026)
1. Nothing may break backwards compatibility with `feature/data-refactor` (`3070b20e8`). Existing code must compile and behave the
   same. Deliberate accuracy improvements are the purpose of the PR; list each one.
2. Deprecating existing interfaces (keeping a thin deprecated wrapper) is allowed when it brings a large simplification.
3. C1: restore the base behaviour, i.e. `J2000`/`ECLIPJ2000` apply the SOFA frame bias. Keep `ICRF` (no bias) and `EME2000` as
   additions.
4. C2: HARDISP (342 constituents) stays the default; the 11-constituent synthesis must stay available.
5. C3: `short_term_eop_interpolation_settings` must keep working (accepted, no error), since it exists in the base.
6. The original plan put mandatory input files in tudat-resources. On 15 October, Dominic instead requested checked-in
   copies and CI installation from this repository; see the resource section above.
7. This follow-up initially required local commits only. On 15 October, after numerical validation, Dominic explicitly
   authorized committing and pushing after the resource update. Include these notes as `docs/pr1092_review_notes.md`
   and then write a reviewer handover.
8. On 15 October, Dominic removed the backwards-compatibility requirement for serialization. The follow-up uses this
   exception to fix existing value comparisons directly; it does not change propagation or Earth-model behavior.

## Status per item
| Item | Status |
| --- | --- |
| Proposed diff items (REVIEW-DONE §6): text bindings, default converter reuse, matrix-constructor argument, nominal-CIP rate, ocean-pole file constructor, relativistic settings collapse, station-creation overload merge, EOP evaluator internals | Done |
| S1: rounded-double boundary machinery removed; 1 µs end-point tolerance (`EarthOrientationDataEvaluator::coverageToleranceSeconds`, station-loading grid) | Done. 26 native tests built and passing (after test updates: sub-µs "outside" checks moved to 1 ms) |
| Compatibility with `feature/data-refactor` | Done and tested: C1 reverted (`J2000`/`ECLIPJ2000` apply the frame bias; `ICRF` identity and `EME2000` = `J2000` kept as additions); C3 short-term interpolation settings accepted again (no effect on the precise model; test added); old spin-only `calculateRotationRateFromItrsToGcrs` templates restored and used again for derivatives of custom-interpolator Earth models (base behaviour); base `BodyCentricToBarycentricGroundStationMotionSettings( centralBody, useGR )` constructor restored, default compatibility = unscaled (Eq. 11.18) so existing code no longer throws; `BodyDeformationModel::calculateDisplacement( time, position )` restored; `OceanTideBodyDeformationSettings::getBlqFiles` restored; `LinearStationMotionModel` default reference epoch restored to the base value |
| Test fix (pre-existing at PR head, test not in the coordinator's selection): `unitTestItrsToGcrsRotationModel` `test_RotationBetweenIntermediateFrames` built its expected matrices from a rounded double UT1 and a TDB TIO locator; now uses the model's precise TT/UT1 | Done |
| S3: Desai–Sibois via resource files | Done, tested and pushed (stage 2): the embedded table (`desaiSibois2016EopData.inc`), the IVS reader (`readHighFrequencyEopData`), `highFrequencyEarthOrientationModel` and the matrix constructor of `ShortPeriodEarthOrientationCorrectionCalculator` are removed. `OceanTideEopModel` (now in `shortPeriodEarthOrientationCorrectionCalculator.h`) selects the `desaiSibois2016*` files instead of the conventional ocean files in the factory. Python bindings for the raw table and `read_high_frequency_eop_data` are removed (new in this PR, not base API) |
| Official ITRF2014/2020 transformation tables moved from `pinnedItrfTransformationTables.h` to resource files | Done, tested and pushed (stage 2): `readOfficialItrf2014/2020TransformationData` read `station_locations/Transfo-ITRF2014_ITRFs.txt` and `Transfo-ITRF2020_TRFs.txt` (byte-identical, same SHA-256) |
| S4: one `numerical_derivatives::computeFourthOrderCentralDifference` (Time or double epoch) replaces the guarded stencils in `bodyDeformationModel.cpp`, the relativistic station motion, `getGeocentricTtRateWithRespectToTdb` and `calculateShortPeriodCorrectionRate` | Done, tested and pushed (stage 3). Same step, pairing and Time epochs, so results for valid inputs are unchanged; the relativistic default step is 16 s, which is the step whose result the PR returned (it checked it against 32 s). Removed: convergence/resolution checks, representability checks and the "restore callbacks to the central epoch" logic (tests pinning these removed) |
| S7: custom-Doodson-callback constructors (ocean tide, IERS 2010) | **Not done**: they exist in `feature/data-refactor` (compatibility rule) |
| S7: station-loading overlap/inventory bookkeeping (`validateStationLoadingComposition`, `allowDeclaredOverlap`, `allowUnknownInventory`, the `hasFittedSeasonalTerms` chain through `StationMotionModel`/`GroundStationState`/combined/history motion) | Done, tested and pushed (stage 4). No per-query loading checks in the generic aggregator any more. Avoiding double-counted signals is now a documented configuration responsibility; the explicit frame/origin/basis checks at construction remain. If a check is wanted, a single check at station creation is the suggested alternative |
| S7: atmospheric endpoint/origin policies, ocean-pole latitude policy, reference-point rate policy | Implemented in `023a5323f`, following Dominic's decisions recorded in the build/test handover; native validation passes |
| S10 (time converter in deformation models) | Not started |
| Stage 1 (all of the above) | Committed as `69dff7426` plus the Python test/notes commit, pushed to `origin/fix/earth_deformation_rotation` |
| `DONE-WITH-EDITING.md` | When all editing stops |

## Not implemented (kept as recommendations; reasons)
- S5 (station-history time model to plain `Time`) and S6 (remove trajectory wrapper classes): API redesign is permitted,
  but these changes are outside this build/test follow-up.
- Python exposure redesign (settings factories instead of mutable classes, one factory per deformation model, `ign_station_history`
  returning settings): permitted, but outside this build/test follow-up.
- S2 complete (remove the custom-interpolator EOP route): would remove a base interface; only possible with deprecation, not yet done.
- S8 (shared reading helpers, per-record provenance), S9 (single grid interpolator), smaller S11 items, the broad
  comment/format pass, and test-generator reorganisation: lower value; not started.
- Behaviour changes already in the PR that are deliberate accuracy improvements and are kept: HARDISP default (C2), the exact
  `a = 1/(1+Z)` in `calculateRotationFromCirsToGcrs` (C7), the corrected IERS 11.19 relativistic station formula, and the
  continuous-TAI EOP interpolation.

## Compatibility notes for maintainers (flag, not changed)
- C4 and C5 were resolved in handover commits `9024aa95b` and `ceade07f5`: factories use the registered rotation directly,
  and frame definitions are processed after station/camera creation as on the base branch. The revised tests cover both.
- Base `LinearStationMotionModel` default reference epoch is `JULIAN_DAY_ON_J2000` used as seconds (looks like a unit error). It is kept for
  compatibility; the environment factory always passes the epoch explicitly.
- The PR replaces the base relativistic station formula (sign/factor errors with respect to IERS 2010 Eq. 11.18/11.19) by the IERS
  formula. This is a deliberate correction; results of existing code using this model change.
- The PR changes the C++ virtual interfaces `StationMotionModel::getBodyFixedStationMotion` (now `Time`) and `BodyDeformationModel`
  (non-virtual `calculateDisplacement`, protected virtual `evaluateDisplacement`). Callers compile unchanged; user-written C++
  subclasses would need updating.

## Verification log
- 14 Oct: proposed diff — `-fsyntax-only` on all changed files: OK.
- 14 Oct: S1 — `-fsyntax-only` on 124 files: OK after test updates. Native build of the 26 test targets started in `/home/dominic/Tudat/pr1092-review-build`.
- 14 Oct: S1 — 26 native tests built in `/home/dominic/Tudat/pr1092-review-build` and run serially: all pass after test updates.
- 14 Oct: compatibility fixes — 26 native tests rebuilt and run: all pass (tests adapted: frame-bias fixtures for J2000/ECLIPJ2000,
  SPICE comparison back to the base GCRS form, station round-trip tolerance 8e-9 m).
- 14 Oct: Python kernel built in the review build; every tudatpy module imports; the PR's focused Python tests (9) pass;
  `tests/architecture` (4) pass. The full `tests/test_tudatpy` suite has 214 failures/19 errors, all environmental: serialization is
  disabled in this build configuration (`TUDAT_BUILD_WITH_SERIALIZATION=OFF`, as in the coordinator's build) and tracking-data tests need
  data files that are not in this build tree. None involves the changed code.
- Note: the SPICE comparison test passes only with the GCRS base (SPICE `J2000` has ICRF axes). With the base behaviour restored,
  tudat's `J2000` Earth-rotation base frame differs from SPICE `J2000` by the frame bias (~23 mas), exactly as in `feature/data-refactor`.
- 14 Oct: stage 2 (S3 + ITRF tables to resources) — 28 native tests (26 + ReferenceFrameTransformations + TerrestrialFrameTransformationData
  reader) pass; kernel rebuilt; all modules import; focused Python + architecture tests (12) pass.
- 14 Oct: stage 3 (S4) — 26 native tests pass after removing assertions that pinned the removed guards.
- 15 Oct: stage 4 (S7 loading bookkeeping) — 26 native tests pass; kernel rebuilt; imports ok; 12 focused/architecture Python tests pass.

### 15 October: handover build and test follow-up

- Validated source: `review/pr1092-external-simplify` at `70aabd7c6`, containing the five local handover commits, plus
  the test corrections described below. This first pass changed only tests; numerical tolerances remain unchanged.
- Most complete build directory: `/home/dominic/Tudat/pr1092-review-build`, configured from
  `/home/dominic/Tudat/tudat-monorepo/tudatpy-pr1092-external-review`. The existing Release/GNU 15.2 configuration
  (`-O1`, serialization OFF) was retained for this first pass. All compilation uses one job pool with at most eight jobs. Native tests run
  serially with numerical-library thread counts set to one.
- All 29 targets listed in handover section 3 compiled without errors. All **209 Boost cases and 14,718 assertions pass**,
  including the 23 cases in `test_simulation_EnvironmentModelSetup` that exercise environment creation after the frame-order
  restoration. Successful earlier runs were retained when only a later test's expectations changed.
- Three test files needed updates for the approved behavior:
  - `unitTestAtmosphericTideDeformation.cpp`: replaced the old rejection of default RP03 conventions with displacement
    and state comparisons against the explicitly configured model, at the same three epochs and original tolerances.
  - `unitTestCompleteEarthStationVelocity.cpp`: physical leap/coverage checks use the required TDB environment configuration.
    Interior factory checks retain all five declared rotation scales and compare against the registered model directly;
    both finite-difference velocity checks and their tolerances are unchanged.
  - `unitTestEarthDeformationTimeRouting.cpp`: solid-tide references use the registered rotation at the supplied environment
    epoch. The basic solid-tide factory is checked against the existing degree-two displacement equation for both TT- and
    TDB-declared rotations. The original position, callback-time and velocity checks retain their tolerances.
- These failures were obsolete RP03-default and hidden-TDB-copy expectations. The revised tests preserve the approved
  behavior from handover section 2; no numerical failure was worked around.
- Before rebuilding the kernel, its pre-handover package was preserved and its provenance checked against the stage-4
  build log. A separate checkout at `92b9c7b16` reproduced the full Python baseline: **214 failed, 568 passed, 109 skipped,
  1 xfailed, 19 errors**. Complete node, phase and traceback records were saved for all 233 failed/error cases.
  The causes are 207 failures with serialization disabled (111 file-IO and 96 pickle), three ATDF failures and 19 ATDF
  setup errors from missing fixtures, and four tracking-residual failures from missing fixtures/kernels.
- Kernel preparation left the CMake cache byte-identical and did not rewrite compiler flags. Recompiled libraries depend
  on refreshed Git/date metadata in `config.hpp`; mathematical switches and source are unchanged. This does not warrant
  repeating the passing native tests.
- Evidence is in `/home/dominic/Tudat/pr1092-review-build/verification-20261015/`: `selected-native-targets.txt`,
  `native-summary.json`, `native-commands.txt`, the per-executable logs and XML reports, `baseline-provenance.json`,
  `baseline.log`, `baseline.xml`, `baseline-reports.json`, and `native-config-flags-audit.txt`. Original failures are retained.
- The first kernel build succeeded. All 14 focused Python/architecture tests passed, and all 209 autosummary names on
  the 11 edited RST pages resolved. The full suite reproduced all baseline totals and all 233 failed/error node, phase
  and normalized exception causes exactly, with no added, removed or changed failures. The remaining failures were
  accounted for by disabled serialization, 23 missing-fixture reports and three missing-SPICE-kernel reports.
- Dominic then required **all tests to pass**. The baseline comparison is therefore not the final acceptance result.
  `TUDAT_BUILD_WITH_SERIALIZATION` was changed to ON, with other settings retained, and real test fixtures were
  recovered locally. The same 29 native executables and the Python kernel were rebuilt with one shared eight-job pool.
  New evidence was saved separately under `/home/dominic/Tudat/pr1092-review-build/verification-allpass-20261015/`.
  No final commit or push had been made at that point.

### 15 October: serialization-enabled verification

- Enabled `TUDAT_BUILD_WITH_SERIALIZATION=ON` in the existing review build. The only other CMake cache change is the
  generated makefile count; compiler, optimization and numerical settings are unchanged. The kernel and the same 29 native
  targets share one `make -j8` graph. Exact commands and before/after cache/compiler flags are saved in
  `verification-allpass-20261015/` under the build directory.
- All 29 native executables rebuilt successfully and passed a fresh serial run: **209 cases, 14,718 assertions,
  zero failures**. This rerun was necessary because serialization changes compiler definitions.
- Restored ten required/dependent fixture files selectively from the repository's `tests/test_tudat/data.zip`: ATDF
  measurements/ramps, MRO troposphere/ionosphere corrections, and the MEX/Voyager/JUICE kernel/measurement pairs.
  Archive and destination hashes agree with existing local caches. No downloads, fabricated data or source changes
  were needed; existing differing files were not overwritten. Both fixture-restoration JSON manifests record provenance.
- Removed eight invalid `pytest.importorskip("spice")` checks from the two serialization test files. Their actual
  `tudatpy.interface.spice` imports remain, so these tests will now execute. The prior 109 skips comprised these eight,
  99 exclusions explicitly selected by `--no-remote-data`, and two existing empty observation-serialization placeholders.
- The existing formatting hooks passed for Python and only joined one C++ declaration across a line break. That
  whitespace-only edit does not affect the verified native binaries. No assertion, tolerance or production behavior changed.
- The serialization-enabled kernel compiled successfully. All 14 focused Python/architecture tests and all 209 autosummary
  names across 11 pages/modules passed. The first full offline run had **801 passed, eight failed, 101 skipped and one XPASS**.
  All 233 original failed/error cases now pass; the eight failures were setup errors in the newly enabled tests. These full
  results are preserved under `verification-allpass-20261015/attempt-01/`. No final commit or push has been made.

### 15 October: serialization equality correction

- Dominic explicitly relaxed serialization backwards compatibility. The correction below is limited to serialization
  equality; propagation, Earth rotation/deformation behavior and numerical tolerances remain unchanged.
- The newly enabled tests also needed their existing setup calls corrected: `SystemOfBodies.create_empty_body`,
  `save_to_binary`/`load_from_binary`, the concrete result loader, and the existing variational solver/parameter factories.
  Hybrid-termination cases now supply the missing RK4 integrator, with the same 30-second step and 600-second arc as their
  companion cases. The obsolete mode-coupled Love-number xfail was removed after its original assertions passed.
- These fixes exposed two result-equality failures. Diagnostics proved that both pickle round trips preserve all public
  state histories, exact `Time` keys, state-transition and sensitivity matrices. Original/restored JSON files are byte
  identical. Every differing binary byte lies in unused long-double padding; no numerical data differ.
  The full evidence is in `verification-allpass-20261015/pickle-equality-diagnostic/diagnostic-summary.md`.
- Corrected the existing equality methods in three headers: `SingleArcSimulationResults` compares its owned processing
  settings and dependent-variable interface by value, and the save/print interval comparisons recognize two NaNs as the
  same unset interval. Comparisons of numerical results remain exact. No new helper layer or comparison tolerance was added;
  the original Python `assert_roundtrip` assertions remain intact.
- Only 49 compiled objects depend on these header changes (26 library objects and 23 kernel objects). The kernel and
  `EnvironmentModelSetup` rebuilt successfully in one shared eight-job graph, without reconfiguration; the CMake cache
  remained byte-identical. The rebuilt environment test passed all 23 cases and 1,114 assertions. Evidence is in
  `verification-allpass-20261015/equality-fix-build/`. The other 28 native tests exercise unchanged code; their passing
  results are retained. Final Python verification is recorded below.

### 15 October: final validation

- Final kernel: `/home/dominic/Tudat/pr1092-review-build/src/tudatpy/kernel.so`, SHA-256
  `3f29c7f14a75591f230b72edcbb8bbf65369934186823a91185697966591e604`. It contains the five handover commits plus the
  reviewed serialization equality correction. Every final Python run asserted this exact kernel and the authoritative
  package/source paths; the recorded test-source hashes identify the corrected tests.
- Nine targeted serialization cases passed, including the original single-arc and variational pickle equality assertions.
  The Love-number case also passed after its obsolete expected-failure marker was removed.
- All **14 focused Python/architecture tests passed**. All **209 autosummary names on 11 pages in 11 public modules resolved**.
- Full `tests/test_tudatpy --no-remote-data`: **810 passed, 101 skipped, zero failures/errors/xfails/xpasses**, in 49.96 seconds.
  The skips are exactly 99 explicit offline exclusions and two unchanged, empty observation-serialization placeholders.
  The only warning is the existing `tudatpy.interface.spice` deprecation. No new skip or expected-failure marker was added.
- Comparison of individual test nodes and phases confirms that all **233 original baseline failures/errors now pass**;
  none is missing or skipped. Exact numerical comparisons and all original tolerances were retained.
- The passing native result remains **29 executables, 209 cases and 14,718 assertions**. Following the equality correction,
  the affected environment executable was rebuilt and passed its 23 cases/1,114 assertions again. Native tests were serial,
  numerical-library runtime threads were one, and compilation never used more than eight jobs in total.
- Final evidence is under `/home/dominic/Tudat/pr1092-review-build/verification-allpass-20261015/`:
  `python-final-summary.md`, `python-final-commands.txt`, `final-baseline-comparison.json`,
  `serialization-final.*`, `focused-final.*`, `full-final.*`, their `*-reports.json` files, `autosummary-final*`, and
  `equality-fix-build/kernel-provenance.json`. Original baseline and intermediate failures remain available separately.
- The subsequent resource/CI update changes data placement and workflow setup only. Dominic explicitly requested no
  recompilation for it; no compiled source or numerical behavior is changed by that update.

### 15 October: resource files included for CI

- Included all seven required files under the repository's `resources/` directory: four Desai–Sibois tables, their existing
  provenance README and the two official ITRF tables. Their 26,899 bytes match the source additions and the files used in
  the completed numerical tests. A scoped `.gitattributes` rule prevents Git from changing their newlines on Windows.
- Added one setup step to `.github/workflows/build_and_test.yml`, after package installation and before both test steps.
  It uses the existing native resource-path binding and `shutil.copytree(..., dirs_exist_ok=True)` to merge the repository
  additions into that directory. Standard resources remain available alongside them; no loader, installer, CMake or
  numerical source change was needed.
- Verified the workflow YAML and embedded Python syntax, step ordering, all seven resource hashes, and the copy operation
  into an empty temporary directory. The final kernel resolves the same native resource location used in the numerical
  validation. The Git text attribute is unset for all seven files, preserving their exact contents across checkouts.
- Evidence: `verification-allpass-20261015/vendored-resource-validation.json`, `vendored-resource-ci.patch` and
  `vendored-resource-untracked-files.txt`. No recompilation or broad test rerun followed this data/workflow-only change,
  as explicitly requested. A separate tudat-resources publication is no longer needed for this CI setup.

### 15 October: whole-repository CI formatting

- The first push passed its local commit hooks, but CI's `pre-commit --all-files` found formatting left in nine earlier
  PR files. Applied only whitespace/wrapping changes and one namespace-closing comment. The executable tokens, strings,
  numerical values, assertions and resource bytes are unchanged.
- The complete cached `pre-commit --all-files` check now passes for both black and clang-format. The reviewed patch and
  both formatter runs are in `verification-allpass-20261015/format-ci*`. No recompilation or numerical rerun is needed
  for these formatting-only changes.

### 15 October: review-comment patch and CI without unity builds

- Applied the supplied review-comment patch unchanged as `da9fbf194`, after the local test-data archive commit
  `3613e38d3`. It documents settings, restores the historical Delta-T default, combines the angles/rates records,
  moves short-period rates into their calculator, templates the polar-motion entry point, replaces the selected
  lambdas and adds the requested braces. The remaining design decisions are unchanged.
- Linux CI job `116421227351` failed when linking `test_ephemerides_CartesianStateExtractor`: a unity-generated
  basic-astrodynamics object required unlinked polyhedron mathematics symbols. Commit `05292ca2b` explicitly sets
  `CMAKE_UNITY_BUILD=OFF` in both CI configure steps. No other CI option changed; the workflow YAML passes validation.
- The existing `/home/dominic/Tudat/pr1092-review-build` was reused with unity explicitly disabled and all other
  cache options retained. Its source checkout is clean and detached at `05292ca2b`; the previous review branch reference
  is preserved. One Sol worker used the entire eight-job allowance. The actual CMake dependency graph rebuilt callers
  of the enlarged angles record, avoiding stale objects with its previous layout.
- Fourteen native targets are selected: the eight Earth-orientation tests; ITRS/GCRS rotation and the CI-failing
  Cartesian-state extractor; complete station velocity, deformation time routing, solid pole displacement and ocean pole
  displacement. Static review found no changed control flow in the other brace-only edits, so their previous validation
  is retained. Python bindings changed only in docstrings/formatting, so no kernel build or Python run was included
  in this follow-up.
- All fourteen selected executables built successfully and passed: **100 cases, 6,508 assertions, zero failures**.
  In particular, `CartesianStateExtractor` links with unity disabled and passes seven cases/twelve assertions;
  rotation derivatives pass thirteen cases/292 assertions. The final ITRS/GCRS rotation test passes six cases and
  1,403 assertions. Each executable ran once, serially, with numerical-library threads limited to one.
- No source fix or numerical/tolerance change was needed. The only CMake cache change was explicitly setting unity
  to OFF, and Earth-orientation compiler flags were byte-identical. No unity compilation unit was built.
- Commands, source identity, cache comparison and per-executable hashes, timestamps and test results are recorded in
  `/home/dominic/Tudat/pr1092-review-build/verification-ci-review-20261015/native-summary.json`, with build/runtime
  logs and XML reports alongside it. The tested source is `05292ca2b`; this later notes update changes no compiled file.
  All changes are local: Dominic explicitly requested no push, and GitHub comments remain untouched. Remote CI has
  not been rerun with these local commits.

### 15 October: second review patch and CI follow-up

- Applied `PR1092-REVIEW-COMMENT-FIXES-2.patch` unchanged as `8cfd4532d`. The C++ calculations now share one EOP-based
  route; the input-data and evaluated-angle structures have their requested names. Python signatures are retained.
  The removed C++ entry points and constructor arguments were checked for stale callers. Unity remains disabled.
- EOP node recognition now allows `2 * epsilon(long double) * 3600` seconds of fractional-hour rounding before selecting
  an interval. The existing outer coverage margin is unchanged. Both-sided node/gap and right-interval rate checks were
  added to existing cases; requests beyond the rounding margin still have to fail across an unavailable stencil.
- Only the five precision-sensitive test files were adjusted. Offset bounds use the fractional-hour resolution; clock
  derivative bounds propagate it through the existing 32-second stencil; differenced-arrival bounds account for both
  receiver clocks. Existing tolerance floors and exact node checks remain on wider-long-double platforms. Boost requires
  matching residual/tolerance types, so bounds are evaluated in long double and converted to double for double residuals.
- The deformation nonfinite-input case uses NaN, which `Time` can represent, instead of an infinity that triggers undefined
  integer conversion before model validation. Scalar step infinity checks remain. No generic `Time` change was made.
- Tracking TXT and IFMS expected values now use the same default EOP-backed clock as collection creation, with their
  original tolerances. The live MPC case selects the established 2022–2023 interval and explicitly requires nonempty data.
  Its existing remote-data marker remains; the observed failure was EOP coverage, not a connection failure.
- Dominic authorized increasing only the older MGS first-block upper residual bound from 0.12 to 0.15 Hz. In CI run
  `38752589589`, both operating-system logs execute all three DSN cases and attribute their only two failures to that
  MGS bound. The MRO case's 59-observation, 2 mHz mean and 3 mHz RMS checks passed there and remain unchanged.
- The rate documentation now describes the single conversion route, the available interpolated TT/TDB clock rate and
  the fixed 32-second nominal SOFA derivative. Stale descriptions of removed APIs and interpolation derivatives were removed.
- Source changes through `4cfe1decf` were first committed without executing code, as instructed. Dominic subsequently
  authorized affected builds and tests. All 22 selected native executables and the Python kernel built successfully in
  `/home/dominic/Tudat/pr1092-review-build`, using a separate detached source checkout and at most eight compiler jobs
  total. Release, serialization ON and unity OFF were retained. Tests ran serially with numerical-library threads set to one.
- Compilation required two explicit observation-factory includes in the DSN test/helper (`711f47037`). The routing suite
  needed three old exception expectations replaced with checks of supported finite clock rates/deformation states
  (`39aff926a`). Both corrections are test-only. Missing local fixtures were extracted from the existing `data.zip`
  without replacing existing files. Only the affected/blocked executables were rerun; initial logs were preserved.
- **All 22 native executables passed: 170 cases and 24,872 assertions.** MRO residual mean was
  **0.0012364529657288123 Hz**, RMS **0.001906654846884383 Hz**, within its unchanged 0.002/0.003 Hz limits.
- **All 64 focused Python checks passed**, with no failures/errors/skips or reported pytest warnings: the live MPC
  integration case with remote data enabled, plus 63 time/EOP/rotation/converter-binding and architecture checks.
  MPC retained its original numerical assertions and passed the nonempty-data check; no connection fallback was needed.
- Tested clean source: `d867b1619238563bf85eac763bc7dbbc4a607085`. Its production source is identical to authoritative
  `39aff926a`; three formatting-only test differences do not change the tested expressions. Python verified the interpreter,
  build package and kernel paths. Kernel SHA-256: `2988970ddfe1fee2fb9ad20d85d07167ffc6639a2c98e861303b38e42ab902bf`.
- Evidence is in `/home/dominic/Tudat/pr1092-review-build/verification-review2-20261015/`: `native-final-summary.json`,
  `python-review2-summary.json`, exact commands, per-attempt logs/XML and source/kernel provenance. Black, clang-format
  and `git diff --check` pass. The full unrelated suite was not rerun; macOS/ARM64 execution remains a CI check.
- No changes have been pushed or GitHub comments modified. The separate reviewer report is
  `/home/dominic/Tudat/tudat-monorepo/PR1092-REVIEW-FIXES-2-REPORT.md`.

### 15 October: shorter tests and complete archived fixtures

- Moved the remaining 39 loose C++ fixture assets into `tests/test_tudat/data.zip`; 37 were missing
  from the archive and two licences were already present. All 111 PR fixture/provenance files now
  have manifest hashes. Published data, generated reference values and existing archive payloads
  retain their exact bytes. The source-tree copies are removed.
- C++ and Python tests use the existing build-tree test-data path. Test targets consuming archived
  headers or files depend on `unpack_test_data`, and their include path contains the extracted data.
- Removed both architecture tests, their standalone CMake companions and five separate reader suites.
  Their useful published-format checks run in the retained model tests. SINEX loading now has three
  focused cases, IMLS loading one and the three GRAIL mass formats share one case.
- The reviewed test sources contain 290 cases instead of 382, and 35 exception checks instead of 321.
  Their source length falls from 16,027 to 13,223 lines, excluding the relocated fixture files.
  Removed cases are repetitive validation, endpoint/leap-second or interface checks. Independent
  IERS/SOFA/HARDISP/Desai/IGN/ITRF references, numerical tolerances, state/velocity checks and the
  actual Mac near-node precision regression remain. No production calculation changed.
- Added purpose and assertion-group comments to the new or substantially changed tests.
  `KeepEmptyLinesAtTheStartOfBlocks` makes clang-format preserve the requested separation before comments.
- Static verification checked archive CRCs, unique names, original bytes, manifest hashes, extracted
  fixture paths, test registration, Python syntax and formatting. No local builds or tests were run,
  following the current instruction. CI is the compile/runtime validation for this cleanup.
- Dominic explicitly authorized committing and pushing this cleanup, followed by CI checks every
  ten minutes. CI outcomes are recorded in the local `REVIEW-IMPLEMENTATION-STATUS.md`.
