> **15 October, later:** four further commits (policies replaced by single behaviours, C4, C5, documentation pass) are
> described in `docs/pr1092_handover_build_test.md`, which also lists the build and test steps still to do.

# PR 1092 review follow-up — implementation status (handover document)

**Editing finished on 15 October 2026.** See `DONE-WITH-EDITING.md`. Final pushed head of `fix/earth_deformation_rotation`: see the
commit list below. The coordinator's checkout `tudatpy-earth-station-models` is behind the remote and must be updated (`git pull`) before
further work; any uncommitted work there must be rebased onto these commits.

## Pushed commits (on top of 226d87de1)
1. `69dff7426` Simplify PR 1092 time handling and restore base compatibility (proposed diff items, S1, compatibility fixes C1/C3 etc.)
2. `b1bf14869` Use the station-history file reader in the binding test; add PR 1092 review and notes
3. `6954ab6e1` Read Desai-Sibois and official ITRF tables from tudat resources (S3; **needs tudat-resources files**)
4. `1969a1c49` Use one fourth-order central-difference helper for station and Earth orientation rates (S4)
5. `12bddd0c6` Remove per-query station-loading overlap bookkeeping from the generic deformation aggregator (S7 part)
6. (this notes update)

Net effect of 1–5 against 226d87de1: native code +575 / −2,046 lines (net −1,471), Python bindings +68 / −179 (net −111), tests +136 / −376 (net −240); the review documents under docs/ are added.

Keep this file current. Another agent must be able to continue from it.

## Locations
- Review worktree: `/home/dominic/Tudat/tudat-monorepo/tudatpy-pr1092-external-review`, local branch `review/pr1092-external-simplify`.
  Base: `226d87de16d2e962cbff7326a830efbc6f89fc3f` (= `origin/fix/earth_deformation_rotation` when fetched, unchanged).
- Review build directory (my own; do not use the coordinator's): `/home/dominic/Tudat/pr1092-review-build`.
  Configured with the same options as `/home/dominic/Tudat/earth-station-python-build`. Build with `make -j8 <targets>` (never more than 8 jobs).
  Selected test targets are listed in `/tmp/claude-1000/-home-dominic/da490dc0-d594-4ff0-aa23-705088333291/scratchpad/selected.txt`
  (26 earth_orientation / ground_stations / ItrsToGcrs / TerrestrialFrame tests).
- Review: `/home/dominic/Tudat/tudat-monorepo/REVIEW-DONE.md`. The first proposed diff is `REVIEW-PROPOSED-CHANGES.diff`.
  The list of riskier redesigns is `REVIEW-RISKIER-REDESIGNS.md`.
- Resource additions for tudat-resources (not yet a PR):
  `/home/dominic/Tudat/tudat-monorepo/tudat-resources-additions/earth_orientation/` (Desai–Sibois 2016 tables). These are already copied
  into the local `~/.tudat/resource/earth_orientation/`.

## RESOURCES CHANGE REQUIRED (flag)
tudat-resources needs these new files. They are installed locally in `~/.tudat/resource` and staged for a PR in
`/home/dominic/Tudat/tudat-monorepo/tudat-resources-additions/`:
- `earth_orientation/desaiSibois2016PolarMotionOceanTidesAmplitudes.txt`
- `earth_orientation/desaiSibois2016PolarMotionOceanTidesFundamentalArgumentMultipliers.txt`
- `earth_orientation/desaiSibois2016Ut1OceanTidesAmplitudes.txt`
- `earth_orientation/desaiSibois2016Ut1OceanTidesFundamentalArgumentMultipliers.txt`
- `earth_orientation/README_desaiSibois2016.txt`
- `station_locations/Transfo-ITRF2014_ITRFs.txt` (SHA-256 6cc1c068…0543, the official IGN file)
- `station_locations/Transfo-ITRF2020_TRFs.txt` (SHA-256 0d8dc872…1a68, the official IGN file)
Without them, selecting Desai–Sibois and the legacy ITRF conversion functions (`convertGroundStationState*Itrf*`) fail at run time.

## Decisions from Dominic (14 October 2026)
1. Nothing may break backwards compatibility with `feature/data-refactor` (`3070b20e8`). Existing code must compile and behave the
   same. Deliberate accuracy improvements are the purpose of the PR; list each one.
2. Deprecating existing interfaces (keeping a thin deprecated wrapper) is allowed when it brings a large simplification.
3. C1: restore the base behaviour, i.e. `J2000`/`ECLIPJ2000` apply the SOFA frame bias. Keep `ICRF` (no bias) and `EME2000` as
   additions.
4. C2: HARDISP (342 constituents) stays the default; the 11-constituent synthesis must stay available.
5. C3: `short_term_eop_interpolation_settings` must keep working (accepted, no error), since it exists in the base.
6. Mandatory input files belong in tudat-resources. Add them locally under `~/.tudat/resource` and flag the resources change.
7. Pushing this work to `fix/earth_deformation_rotation` is authorized by Dominic. Rebase onto the latest remote; never force-push.
   Include this notes file in the repository commit (as `docs/pr1092_review_notes.md`).

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
| S7: atmospheric endpoint/origin policies, ocean-pole latitude policy, reference-point rate policy | Not done (each changes which behaviour users get; left as review recommendations) |
| S10 (time converter in deformation models) | Not started |
| Stage 1 (all of the above) | Committed as `69dff7426` plus the Python test/notes commit, pushed to `origin/fix/earth_deformation_rotation` |
| `DONE-WITH-EDITING.md` | When all editing stops |

## Not implemented (kept as recommendations; reasons)
- S5 (station-history time model to plain `Time`) and S6 (remove trajectory wrapper classes): large public API redesign of new PR
  features; needs Dominic's go-ahead.
- Python exposure redesign (settings factories instead of mutable classes, one factory per deformation model, `ign_station_history`
  returning settings): large user-facing change; needs Dominic's go-ahead.
- S2 complete (remove the custom-interpolator EOP route): would remove a base interface; only possible with deprecation, not yet done.
- C4 (hidden TDB re-wrapping of the Earth rotation) and C5 (frame-definition ordering in `createSystemOfBodies`): need a design
  decision; flagged only.
- S8 (shared reading helpers, per-record provenance), S9 (single grid interpolator), smaller S11 items, the broad
  comment/format pass, and test-generator reorganisation: lower value; not started.
- Behaviour changes already in the PR that are deliberate accuracy improvements and are kept: HARDISP default (C2), the exact
  `a = 1/(1+Z)` in `calculateRotationFromCirsToGcrs` (C7), the corrected IERS 11.19 relativistic station formula, and the
  continuous-TAI EOP interpolation.

## Compatibility notes for maintainers (flag, not changed)
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
