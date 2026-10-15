# PR 1092 — independent review: simplification, code paths, exposure and style

Reviewed on 14 October 2026, independently of the coordinator's own static review.
`MODIFICATIONS-DONE.md` was read as background information only; where my
recommendation differs from the restrictions recorded there, I say so explicitly
and give the numerical argument.

## 0. Review basis and deliverables

| Item | Value |
| --- | --- |
| PR | https://github.com/tudat-team/tudatpy/pull/1092/ |
| Reviewed (starting) commit | `226d87de16d2e962cbff7326a830efbc6f89fc3f` (branch `fix/earth_deformation_rotation`) |
| Comparison base | `3070b20e83e18a749767fe0332a5280464f59b16` (`origin/feature/data-refactor`) |
| Proposed changes (diff) | `/home/dominic/Tudat/tudat-monorepo/REVIEW-PROPOSED-CHANGES.diff` — a diff **against `226d87de1`** |
| Worktree holding the same changes (uncommitted, local branch only) | `/home/dominic/Tudat/tudat-monorepo/tudatpy-pr1092-external-review`, branch `review/pr1092-external-simplify` (never pushed) |

Apply with: `git -C /home/dominic/Tudat/tudat-monorepo/tudatpy-earth-station-models apply --index /home/dominic/Tudat/tudat-monorepo/REVIEW-PROPOSED-CHANGES.diff`
(on a clean checkout of `226d87de1`).

**Follow-up implementation (in progress after this review was written).** Further
changes from Section 2 are being implemented in the same review worktree
(`/home/dominic/Tudat/tudat-monorepo/tudatpy-pr1092-external-review`). It is local only: no
commits, nothing pushed, and the authoritative checkout is untouched. The cumulative diff against
`226d87de1` will be at `/home/dominic/Tudat/tudat-monorepo/REVIEW-IMPLEMENTATION.diff`. Progress
and verification status are in `/home/dominic/Tudat/tudat-monorepo/REVIEW-IMPLEMENTATION-STATUS.md`.
Use the diff above for the items in Section 6; check the status file before starting work on
S1/S4/S7 so the same items are not implemented twice.

**Verification of the diff.** Nothing was built, linked or run. Each changed C++ file in the diff,
plus several files that include the changed headers, passed a compiler syntax check
(`-fsyntax-only`, using the existing compile flags from `/home/dominic/Tudat/earth-station-python-build/compile_commands.json`
but pointed at the review worktree). This covers the native sources, the pybind11 sources and the
changed unit tests. Nothing was written into the build directory. The unit-test edits are
straightforward, but the changed tests were **not** executed.

**Scope of the PR (base → head).** Native code +13,108 / −1,450 lines, Python bindings +1,583,
tests +18,932, docs +1,880. My main conclusion is that the native part can lose several thousand
lines without numerical loss. The largest contributors to avoidable complexity are:

1. A system-wide mechanism to recognize "the rounded `double` image of an exact EOP/loading/civil
   boundary". It is the reason for most double/`Time` duplicate overloads (Section 2, S1).
2. Two complete Earth-orientation data routes: the new continuous EOP evaluator and the old
   custom daily interpolators. They are kept side by side in every class, although after this PR
   no production code constructs the old route (S2).
3. Six hand-written, heavily guarded finite-difference stencils (S4).
4. A station-history time model with split source coordinates, convention strings, optional rates
   and user time callbacks. Its numerical benefit is below 1e-7 m (S5).
5. Many explicit "policy" enums, booleans and alternative syntheses, most with a single sensible
   choice (S7).
6. Provider data embedded in source (Desai–Sibois) with its own reader. The existing tudat mechanism
   (resource files plus `EopCorrectionSettings`) already does the same job (S3).

---

## 1. Required corrections and decisions (address before merge)

These are behaviour changes or inconsistencies that need a maintainer decision or a fix. They are
not optional cleanup.

**C1. Earth rotation `J2000` / `ECLIPJ2000` base frames changed meaning for every user.**
`include/tudat/astro/ephemerides/itrsToGcrsRotationModel.h:55-75`. Before this PR, `J2000` and
`ECLIPJ2000` applied the SOFA frame bias (≈23 mas). They no longer do, and the old behaviour is
renamed `EME2000`. I think the new behaviour is physically right: SPICE `J2000` vectors are ICRF
axes, so the old bias put the Earth rotation and the SPICE ephemerides on different axes. However:

- this silently changes every existing Earth-rotation simulation by ≈23 mas, which is ≈0.7 m at the
  Earth's surface and tens of metres at lunar distance;
- it is unrelated to station deformation.

Required: an explicit maintainer decision, a changelog/"breaking change" note, and preferably a
separate PR. `EME2000` is a new frame name used nowhere else in tudat. If it is kept, it should
also be added to the frame documentation.

**C2. The default ocean-tide synthesis changes** from 11 constituents to the 342-constituent
HARDISP expansion (`expose_shape_deformation.cpp`, `ocean_tidal`). This is a millimetre-level
change in every existing ocean-loading simulation and needs a changelog note. I recommend also
deleting the legacy synthesis (S7).

**C3. `short_term_eop_interpolation_settings` now throws.**
`createRotationModelEarthOrientation.cpp:30-35`. The argument is still exposed in
`environment_setup.rotation_model.gcrs_to_itrs`, so existing scripts that pass it now fail at
environment creation. Either remove the parameter from the C++/Python signature (and say so in the
changelog), or keep it working. Keeping it is easy: interpolate the short-period corrections as
before; the UT1/TT arguments are the same either way. The factory also still passes
`getShortTermInterpolatorSettings( )` (always null by then) to the correction calculators in the
`else` branch (`:93-106`), which is dead code.

**C4. Hidden re-wrapping of the Earth rotation model in TDB.**
`createGroundStations.h:233-238` (`createReferencePointEphemeris`), `createBodyDeformationModel.cpp:31-38`
(`deformationRotationInTdb`) and `createGroundStations.cpp:~286` (relativistic motion) each
`dynamic_pointer_cast` the body's rotation to `GcrsToItrsRotationModel` and build a second model
with `tdb_scale`. A user who sets a non-TDB `inputTimeScale` therefore gets different rotations
in different parts of the same environment. Decide one rule and enforce it in one place. Inside an
environment every rotation is called with TDB, so the simplest rule is that
`createRotationModel` rejects a non-TDB `inputTimeScale` for a body rotation. Then all three
re-wrappings disappear.

**C5. The global ordering of `processBodyFrameDefinitions` changed in `createSystemOfBodies`.**
`createBodiesFactory.h:291-293`. Frame processing now runs before ground-station creation. This
affects every environment, not only stations. Make sure a test covers an environment with ground
stations and a non-SSB global origin (ephemeris origin conversion). Otherwise, restrict the
dependency to the relativistic motion factory by building its callbacks lazily.

**C6. `createDefaultTimeConverter` decided reuse by comparing every settings field.**
`terrestrialTimeScaleConverter.cpp:282-288`. Adding a field to `EarthOrientationEvaluationSettings`
later would silently share a converter with the wrong settings. **Fixed in the diff:** the settings
argument is now `std::optional`, and only the "settings not given" case is shared.

**C7. `calculateRotationFromCirsToGcrs` now uses the exact `a = 1/(1+Z)`** instead of the series
`0.5 + (X²+Y²)/8`. This is a small numerical improvement (≲1e-17 in the matrix today) and is fine,
but it changes a pre-existing function and should be mentioned in the PR description.

**C8. A dead pre-existing API was modified.** `earthOrientationCalculator.h:87-135`: the templated,
spin-only `calculateRotationRateFromItrsToGcrs` overloads have no remaining caller. The PR changed
their numerics (dropped `ut1 - 1.42E-3`). **Fixed in the diff:** the three dead overloads are
removed, leaving only the complete product derivative.

---

## 2. Major simplifications (recommended; ranked by complexity removed)

### S1. Remove the "rounded double boundary image" mechanism (largest single simplification)

**What it is.** About 35 sites recognize a `double` epoch whose value equals
`static_cast<double>(exact boundary)` and replace it by that exact `Time` boundary:

- `EarthOrientationDataEvaluator::convertTtToTai(double)` (`earthOrientationDataEvaluator.cpp:60-75`);
- the `if constexpr( double )` block and `candidateWindow` in `TerrestrialTimeScaleConverter::getTaiTime`
  (`terrestrialTimeScaleConverter.h:528-580`);
- `invertContinuousTime(..., roundedDoubleInput)`;
- the `convertTimeScale` double branch (`terrestrialTimeScaleConverter.h:815-822`);
- the templated `getRotationAnglesAndTimesFromItrsToGcrs<TimeType>`;
- `PolarMotionCalculator::getPositionOfCipInItrsAtTt<TimeType>`;
- `basic_astrodynamics::normalizeScalarTdb`, with its ambiguity check
  (`earthDeformationTime.h/.cpp`, which then only keeps `getGeocentricTimeConverter`);
- the virtual `BodyDeformationModel::normalizeScalarTdb`, overridden by
  `StationLoadingDeformation::normalizeScalarTdb` and its second (TDB) knot grid
  (`stationLoadingDeformation.cpp:115-150, 305-320`);
- `getIgnNominalUtcStationTimeFromTdb(double)` with civil-midnight image snapping
  (`stationHistoryTime.cpp:121-155`).

Above all, it requires a separate **double overload of every `StationMotionModel` entry point**
(`groundStationState.h:40-70`, overridden in `BodyDeformationStationMotionModel`,
`BodyCentricToBarycentricRelativisticStationMotion`, `CombinedStationMotionModel`
(`evaluateMotion<Epoch>`, `evaluatePosition<Epoch>`) and `StationHistoryMotionModel`
(three templated helpers)).

**Why it can go without numerical loss.** `Time(double)` represents any `double` exactly, so for
any epoch strictly inside coverage the double and `Time` paths already give the same result. The
mechanism only matters for an epoch that is the rounded image of the **first or last** EOP (or
loading) sample, which `Time(double)` places at most half a double ULP (≈60 ns in 2026) outside
coverage. Replace all of the above by **one** tolerance in the coverage check:
`EarthOrientationDataEvaluator::evaluate` (and the station-loading grid) accepts
`|tai − boundary| ≤ 1 µs` and evaluates the boundary interval there (extrapolating at most 1 µs).
The largest possible effect is:

- UT1−TAI: |dUT1/dt| ≲ 3e-8 s/s, so ≤ 3e-14 s;
- polar motion: ≲ 1e-13 rad/s, so ≤ 1e-19 rad;
- loading series: ≲ 1e-8 m/s, so ≤ 1e-14 m.

These are far below double resolution of the quantities themselves. The `MODIFICATIONS-DONE.md`
instruction "do not add broad coverage epsilons" is understandable, but a 1 µs edge tolerance is
not broad: no input that differs from a valid one by more than the numerical noise is accepted.
In exchange, the double/`Time` overload pairs reduce to **one non-virtual double→`Time`
forwarding function in the `StationMotionModel` base class**. All overrides become `Time`-only.

**What else disappears.** The "ambiguous double endpoint" errors and their tests:
`unitTestEarthDeformationTimeRouting.cpp` (most of its 461 lines), the boundary-image cases in
`unitTestEarthOrientationDataEvaluator.cpp`, and the endpoint cases in
`unitTestCompleteEarthStationVelocity.cpp`.

### S2. Keep a single Earth-orientation data route

After this PR, every factory (`createGcrsToItrsRotationModel`, `createStandardEarthOrientationCalculator`,
`createDefaultTimeConverter`) builds the continuous `EarthOrientationDataEvaluator`. No production
code constructs `PolarMotionCalculator`, `PrecessionNutationCalculator` or
`TerrestrialTimeScaleConverter` with custom daily interpolators any more; only the unit-test
"cookbook" calculator does. Yet each class keeps both routes:

- `if( evaluator_ ) … else …` in `PolarMotionCalculator::getPositionOfCipInItrs` (both overloads),
  `getPositionOfCipInItrsAtTt`, `PrecessionNutationCalculator::getPositionOfCipInGcrs` (both
  overloads), `EarthOrientationAnglesCalculator::getRotationAnglesFromItrsToGcrs` and
  `getRotationAnglesAndTimesFromItrsToGcrs`, and in the converter `getCurrentTime`,
  `getCurrentTimeDifference`, `updateTimes`, `getUt1MinusUtc`…;
- getters that throw in the other mode (`getDailyIersValueInterpolator`, `getDailyCorrectionInterpolator`,
  `getDailyUtcUt1CorrectionInterpolator`);
- `usesContinuousEop`.

Recommendation: make the evaluator the only daily-EOP source in these three classes, and remove the
interpolator members, constructor arguments and branches. Keep the **EOP-free** default
`TerrestrialTimeScaleConverter( )`, which is used widely (ODF/PSF processing, troposphere
corrections) for atomic/UTC conversions. That remains a valid single route: UT1 requested without
EOP data simply throws. Port the cookbook test fixture to a small `EarthOrientationData` built in
code. This is a C++-only API change (none of these constructors is exposed to Python).

*Done in the diff (first step):* the legacy custom-interpolator **rate** path
(`calculateRates` 5-point stencil, `requireLegacyStencilInSingleInterval`, about 60 lines) is replaced
by an early "rates require a model constructed from EOP data" error. The cookbook-based
derivative-mapping test now uses `createStandardEarthOrientationCalculator()`.

### S3. Desai–Sibois (2016) as resource files, not as embedded data plus a reader

The PR embeds the 159-term table in `src/tudat/io/desaiSibois2016EopData.inc` (166 lines) and adds:

- a reader (`readHighFrequencyEopData.h/.cpp`, 150 lines) and its structs;
- `highFrequencyEarthOrientationModel.h/.cpp` (108 lines);
- a new matrix constructor in `ShortPeriodEarthOrientationCorrectionCalculator`;
- the `OceanTideEopModel` enum and settings accessors;
- a factory branch that compares default `EopCorrectionSettings` to detect custom files
  (`createRotationModelEarthOrientation.cpp:66-82`);
- three Python bindings and an architecture test.

tudat already handles this with resource files read by `readAmplitudesAndFundamentalArgumentMultipliers`.
The IERS libration and ocean tables are just lists of files in `EopCorrectionSettings`
(`createRotationModel.h:304-322`). Converting the Desai–Sibois table into the same two-file format
(ocean amplitudes and multipliers, for polar motion and for UT1) in tudat-resources and selecting
those file names would remove all of the above. The model enum can stay as the user-facing switch
in `gcrs_to_itrs`, where it just selects file lists.

Numerically this is identical: the same coefficients, the same `[GMST+π, l, l', F, D, Ω]`
argument convention, the same unit factors (µas → rad, µs → s), the same group order (libration,
then ocean), and the same RSS amplitude cut-off criterion. This needs a tudat-resources addition,
which is the established place for such tables.

### S4. One numerical-derivative helper for `Time` epochs

Before the diff there were six separate four-point stencils with their own guards:

- `bodyDeformationModel.cpp:15-95` (`differentiateDisplacementAtEpoch`, including a second
  "resolution step" stencil and a "restore callback" catch block);
- `groundStationState.cpp` `guardedStationRate` (coarse/fine convergence, a resolution check and a
  restore callback);
- `PrecessionNutationCalculator::getNominalPositionAndRate` (**simplified in the diff**, see below);
- `TerrestrialTimeScaleConverter::getGeocentricTtRateWithRespectToTdb`;
- `calculateShortPeriodCorrectionRate` (`earthOrientationRates.h`);
- the legacy `calculateRates` (**removed in the diff**).

Recommendation: one small function in `tudat/math/basic/numericalDerivative.h`. It takes a `Time`
centre and a `double` step, forms the four epochs in `Time`, and returns
`(8(f₊₁−f₋₁) − (f₊₂−f₋₂)) / (12h)` (paired differences, as now, so there is no
cancellation penalty). Use it at all remaining sites. Drop:

- the representability checks: offsets of ±32 s or ±1 s on a `Time` are always exact;
- the "resolution step" double check, which guards against callbacks that round to `double`.
  Production callbacks no longer do this after the PR, and those that do are a callback bug;
- the coarse/fine convergence tests. A fixed step chosen for the known signal bandwidth is what
  tudat does elsewhere;
- the "restore callbacks to the central epoch" `catch(...)` blocks. Callbacks are functions of
  time; any caching inside them is the callback's own concern.

The existing `computeCentralDifferenceFromFunction` sums weighted samples rather than paired
differences. For the quantities involved here (≤1 m displacements, unit matrices, ms-level clock
offsets) the extra rounding is <1e-17 per second of step, so even that helper would be
acceptable. A paired-difference `Time` overload is cleaner, though.

*Done in the diff:* `getNominalPositionAndRate` now always differentiates the **direct SOFA model**
with one 32 s stencil. This replaces about 120 lines that reconstructed derivatives of arbitrary
linear/Lagrange/spline interpolators from interior samples. It also removes
`setNominalDerivativeTimeStep`, `getNominalDerivativeTimeStep`, `getNominalCipInterpolator` and the
catch-and-restore block. With the default settings (no CIO interpolation) the numbers are
**bit-identical** to the PR: same stencil, same `Time` epoch arithmetic, same pairing. With CIO
interpolation enabled, the rate is now that of the true model rather than the interpolant. The
difference is the derivative of the interpolation error, which is ≲1e-16 rad/s for the default
interpolation step, i.e. <1e-9 m/s at the Earth's surface.

### S5. Station-history time: evaluate with `Time`, not with `StationHistoryTime`

`StationHistoryTime` (`stationHistoryTrajectory.h:22-34`) carries TAI, an integer MJD, an MJD
fraction, a calendar year, a year fraction, two optional rates and a convention string. Users must
build it (it is exposed with eight mutable fields), and it comes with:

- public IGN-specific converters (`getIgnNominalUtcStationTimeFromTdb` ×2, `…FromTai`,
  `getIgnNominalUtcStationTime`, `convertSinexCivilUtcToTai`), also exposed in Python;
- a custom `std::function< StationHistoryTime( const Time&, bool ) >` callback threaded through the
  settings, the motion model and `createStationReferencePoint.h`;
- convention-string checks in `validateTime`.

Numerically, the split UTC civil coordinates versus plain TAI seconds differ only by the leap
seconds accumulated between the reference epoch and the evaluation epoch (≤ about 10 s over the
SINEX era):

- linear term: 5 cm/yr × 10 s / 3.16e7 s ≈ 1.6e-8 m;
- seasonal phase: 5 mm × 2π × 10 s / 3.16e7 s ≈ 1e-8 m;
- logarithmic post-seismic term right after an event: amplitude/τ × 10 s ≈ 1e-7 m for 0.1 m / 0.1 yr.

Recommendation:

- `StationHistoryTrajectory::getPosition/getState( const Time& tai )`, with elapsed time
  `(tai − tai_ref)` in `Time` divided by `JULIAN_YEAR`;
- the calendar-year phase computed internally from the TAI epoch;
- SINEX civil epochs converted to TAI once at construction with the existing
  `sofa_interface::convertUTCtoTAI`;
- the rate with respect to TDB applied once in `StationHistoryMotionModel` (dTAI/dTDB);
- remove the struct, the custom time callback and the IGN helpers from the public C++/Python API.

If strict IGN civil-day semantics are wanted for traceability, keep the civil-day arithmetic as a
private detail of the trajectory. This still removes the public struct, the optional rates, the
callbacks and the double overload.

### S6. Replace the "sometimes-a-proxy" `StationHistoryTrajectory` base

`TransformedStationHistoryTrajectory` and `StationReferencePointTrajectory` derive from
`StationHistoryTrajectory`, which therefore has a nullable `source_`. Five getters branch on it
(`source_ ? source_->x() : x_`), plus there is a protected proxy constructor and virtual
`getPosition/getState/getMetadata` overrides (`stationHistoryTrajectory.h:105-160`).
`getStationReferencePointIdentity` then downcasts through the chain. Simpler design:

- `StationHistoryTrajectory` is a plain, final data/evaluation class (no `source_`, no virtuals);
- `StationHistoryMotionModel` holds the trajectory plus an optional `TerrestrialFrameTransformationModel`
  (with target identity) and an optional `StationReferencePointOffsetHistory`, applied in that
  fixed order: frame transformation first, then the offset. The required ordering is "transform
  the marker, then add the eccentricity expressed in the target frame". Arbitrary interleavings are
  not a physical requirement;
- the metadata of the evaluated point is computed once in the motion model.

This removes both wrapper classes, their factories, the five branching getters and two Python classes.

### S7. Reduce policies, alternative syntheses and constructors

Each of these doubles a code path for an option without a demonstrated user:

| Location | Recommendation |
| --- | --- |
| `OceanTideSynthesis::legacyElevenConstituents` and 6 `OceanTideEarthDeformation` constructors (`oceanTideEarthDeformation.h:56-85`) | Keep HARDISP only; one data constructor plus the existing file-list convenience. Remove the custom Doodson-argument constructors and `usesNativeUtcPhase_`. |
| `Iers2010EarthDeformation` custom `doodsonArgumentFunction` with `usesNativeUtcPhase_` | Keep the native arguments only (no remaining production caller passes a callback). |
| `AtmosphericTideEndpointPolicy::requireConsistency` (default) | The default **rejects the official RP03 grid** because of its seam/pole rounding. Always average and report the discrepancy; drop the enum. |
| `AtmosphericTideEvaluationSettings` lookup-latitude/basis default `unspecified` (must be set explicitly), `AtmosphericLoadingOriginPolicy` | Choose the provider's documented convention as the default; keep at most the origin conversion as a bool. |
| `OceanPoleTideLatitudePolicy` | The official grid ends at ±89.75°, and no station lies beyond that. Keep one behaviour (clamp) and drop the enum and its Python bindings. |
| `StationReferencePointRatePolicy` | Keep one behaviour (right-hand rate at a switch, as everywhere else in the PR). |
| `StationLoadingEvaluationSettings::allowUnknownInventory/allowDeclaredOverlap`, `StationLoadingSignalPresence` metadata (diurnal/semidiurnal/annual/semiannual/mean/trend), `validateStationLoadingComposition` | Remove the overlap bookkeeping from the generic deformation aggregator (`bodyDeformationModel.h:166-174` calls a loading-specific validator on every query, and it is declared in two headers). Double counting is a configuration issue: document it, or check it once in `createGroundStationMotionModels`. |
| `OceanPoleTideDeformation( std::string gridFile, … )` | **Removed in the diff** (unused; settings load the grid eagerly). |
| `BodyCentricToBarycentricRelativisticStationMotion` legacy constructor plus `legacyExternalPosition_`/`legacyExternalGm_` and their branches | **Removed in the diff** (no caller). |
| Relativistic settings: 2 constructors, `centralBodyName_`, `useGeneralRelativisticCorrection_`, `requireEarthHost_`, `earthRelativisticStationMotionSettings` and Python `earth_relativistic_station_motion` plus `relativistic_station_motion` | **Collapsed in the diff** to one settings constructor, one factory `relativisticStationMotionSettings( externalBodies = {"Sun"}, compatibility = tt_to_tdb )` and one Python function `relativistic_station_motion`. The pre-existing `bodycentricToBarycentricStationMotionSettings` is kept as a thin forwarder. |
| `createGroundStation`/`createGroundStationState` duplicate "with bodies" overloads | **Merged in the diff** into default arguments. Python `add_ground_station` keeps its three call forms. |
| `ShortPeriodEarthOrientationCorrectionCalculator` matrix constructor `interpolation` argument | **Removed in the diff** (never passed). Remove the constructor entirely with S3. |

### S8. Shared reading helpers in `tudat/io`

Each of the eight new readers defines its own `trim`, `number`/`integer` (strict full-token parse,
Fortran `D` exponents), `require`/`fail`/`error`, tokenizer and `leapYear`:

- `readEarthOrientationData.cpp:19-80`;
- `readAtmosphericTideData.cpp:17-55`;
- `readOceanLoadingData.cpp:18-70`;
- `readOceanPoleTideData.cpp:18-36`;
- `readSinexStationHistory.cpp:20-70`;
- `readStationLoadingData.cpp:20-90`;
- `readTerrestrialFrameTransformationData.cpp:17-45`;
- `readHighFrequencyEopData.cpp:18-24`.

Put one set in a non-installed `src/tudat/io/parsingUtilities.h`, and use
`basic_astrodynamics::isLeapYear` (`timeConversions.h:450`). Also drop the **per-record provenance**
(`sourceLine`/`sourceRecord` on every sample, site and parameter, exposed in Python). Keep the
file-level `source` and header strings. Error messages can report the line number while parsing,
without storing it in every record.

### S9. One latitude/longitude grid interpolator

`AtmosphericTideGridInterpolator::interpolate` (`atmosphericTideLoading.cpp:150-196`) and
`interpolateValidatedOceanPoleTideCoefficients` (`oceanPoleTideDeformation.cpp:40-105`) both do
periodic-longitude bilinear interpolation, with different grid storage orders (latitude-fast versus
longitude-fast) and different grid descriptions (axis vectors versus first/step/count). Use one
grid description and one private bilinear function, templated on the coefficient matrix type
(3×4 or 3×2). Keep the atmospheric pole special case (Cartesian averaging) as a pre-processing
step on the grid data, as now.

### S10. Time-converter plumbing into deformation models

Every deformation model now takes a `TerrestrialTimeScaleConverter`. The factories retrieve it with
the same `dynamic_pointer_cast< GcrsToItrsRotationModel >(...)->getAnglesCalculator()->getTerrestrialTimeScaleConverter()`
eight times (`createBodyDeformationModel.cpp`, `createGroundStations.cpp`). Apart from S1 (scalar
boundaries), the models use it only for TDB→TT/TAI, which is EOP-free and geocentric, and for its
rate. The **only** model that needs EOP-based time is atmospheric loading (UT1 phase). After S1:

- remove the converter from the `BodyDeformationModel` base class and from the solid, pole,
  ocean-pole, ocean-tide, basic-tidal and station-loading models (use `getGeocentricTimeConverter()`
  or `sofa_interface` directly);
- pass it explicitly only to `AtmosphericTideDeformation`;
- one helper `getEarthTimeScaleConverter( body )` in `createRotationModel` replaces the repeated casts.

### S11. Smaller code-path and clarity items

- `EarthOrientationDataEvaluator`. **Done in the diff:** the interpolation scratch `std::vector`s
  became `std::array`, and the one-element-Eigen adaptation of scalar quantities was replaced by
  one generic difference-from-anchor lambda that works for `double` and `Eigen::Vector2d`. The
  arithmetic order is unchanged, so the results are identical. Optional further step:
  `validateCubicUt1Monotonicity` builds its own polynomial coefficients. It could reuse the weight
  computation by evaluating dUT1/dTAI at u = 0, ½ and 1 (the derivative of a cubic is quadratic).
- `BodyDeformationStationMotionModel::evaluateMotion` evaluates every model twice for an optional
  state request (a position preflight, then the state). Evaluate the state directly and treat an
  `std::out_of_range` from it as "not available". A position that is available while its rate is
  not is then only an error when `throwExceptionWhenNotAvailable` is set, which is the semantics the
  aggregator documents anyway.
- `TerrestrialTimeScaleConverter::getCurrentTime` contains an explicit identity `switch` that
  re-validates the scale (`terrestrialTimeScaleConverter.h:147-158`). `getCurrentTimeDifference`
  calls `getCurrentTime` "to validate the scale only". Both can be dropped; unknown scales fail
  later anyway.
- `EOPReader( const EarthOrientationData& )` together with `getRawEarthOrientationData` turns the
  legacy reader into a wrapper around the new data. `createStandardEarthOrientationCalculator(
  EOPReader )` and `createDefaultTimeConverter( EOPReader )` only use `getRawEarthOrientationData`.
  Change their parameter to `EarthOrientationData` and keep `EOPReader` only for the legacy map
  getters (or deprecate it).
- `GcrsToItrsRotationModelSettings::setEopData` plus the `eop_data` factory argument plus Python
  `set_eop_data`: keep only the factory argument.
- `createStationReferencePoint.h` contains one inline function. Fold it into `createGroundStations.h`.
- `earthDeformationTime.h` (namespace `basic_astrodynamics`, in `ground_stations`) keeps only
  `getGeocentricTimeConverter` and `validateEarthDeformationUtcInterval` after S1.
  - `getGeocentricTimeConverter` belongs in `terrestrialTimeScaleConverter.h`.
  - `validateEarthDeformationUtcInterval` exists only because the ocean-loading and solid-tide
    **rates** are finite differences of UTC-phase models across leap seconds. For ocean loading,
    the displacement is a harmonic sum, so its rate can be computed **analytically**
    (Σ −Aω sin(ωt+φ)·dUTC/dTDB). That removes the stencil and the leap-second interval check for
    that model. Solid tides need the check only for their step-2 (UTC-phase) part.

---

## 3. Exposure review (location and style)

**C++ locations** are mostly consistent with tudat:

- readers in `tudat/io` (`input_output` namespace);
- runtime models in `astro/ground_stations` (deformation models in `basic_astrodynamics`, as for
  the existing tidal models);
- Earth orientation in `astro/earth_orientation`;
- frame transformations in `astro/reference_frames`;
- settings and factories in `simulation/environment_setup`.

The exceptions are `earthDeformationTime.h` (S1/S11), `createStationReferencePoint.h` (S11), and
**provider data embedded in source** (`desaiSibois2016EopData.inc`, `pinnedItrfTransformationTables.h`).
tudat keeps such tables in tudat-resources (S3). The two ITRF tables are small and replace
previously hard-coded tables, so keeping them in code is acceptable.

**Python locations.**

- `data_input.environment_data.{earth_orientation, ocean_loading, atmospheric_loading, station_loading, station_history, terrestrial_frames}`:
  matches the existing pattern (`ilrs`, `space_weather`). Remove the `read_*_text` functions:
  they are not used by any Python test or document and exist only to parse strings
  (**removed in the diff**, together with their docs entries).
- `dynamics.environment`: the runtime station classes are fine there. The provider-specific clock
  helpers (`ign_nominal_utc_station_time*`, `sinex_civil_utc_to_tai`) and the eight-field
  `StationHistoryTime` do not belong in the public environment API (S5).
- `dynamics.environment_setup.ground_station.ign_station_history` (and the reference-point
  constructors) return **runtime objects** from a setup module. The tudat pattern is that setup
  functions return *settings*, and `create_system_of_bodies` builds the models. Recommend
  `station_history_motion( solutions, discontinuities, site, post_seismic=None, seasonal=None, ... )`
  returning a `GroundStationMotionSettings` that holds the preloaded data.
- `dynamics.environment_setup.shape_deformation`: eight new factory functions (`ocean_tidal_from_data`,
  `ocean_pole_tidal`, `ocean_pole_tidal_from_data`, `atmospheric_tide_from_grid`,
  `…_from_stations`, `…_from_grid_file`, `…_from_station_file`, `station_loading`). The usual tudat
  style is one factory per model. Recommend one function per model taking the data object, e.g.
  `ocean_pole_tidal( grid )`, `atmospheric_tidal( grid_or_station_data, … )`. Users read files
  with the `data_input` readers, which is the PR's own "preloaded data" design. Keep the existing
  `ocean_tidal( blq_files )` for compatibility. Use the existing naming (`…_tidal`).
- Settings exposed as mutable classes with `py::init<>` plus `def_readwrite`
  (`EarthOrientationEvaluationSettings`, `StationLoadingEvaluationSettings`,
  `StationLoadingOriginTranslation`, `AtmosphericTideEvaluationSettings`, `StationHistoryTime`,
  `StationReferencePointOffset`, …) are not the tudatpy pattern. That pattern is factory functions
  with keyword arguments, or keyword arguments on the model factory, e.g.
  `gcrs_to_itrs( …, eop_data=None, eop_interpolation=…, maximum_eop_gap=…, … )`. Raw data records
  from readers may stay inspectable.
- `astro.frame_conversion`: adds two classes (`TerrestrialFrameAffineTransformation` with mutable
  fields, `inverse`, `compose_after`; `TerrestrialFrameTransformationModel`). The rest of that
  module is free functions. A single free function
  `transform_terrestrial_state( state, data, source_frame, target_frame, decimal_year, decimal_year_rate )`
  covers the user need. The affine algebra can stay C++-only.

**Docstring style.** Most new bindings use one-line plain strings or dense one-paragraph strings
rather than the `R"doc(` numpy style with `Parameters`/`Returns` used across tudatpy.
`GcrsToItrsRotationModelSettings`, `IAURotationModelSettings` and `add_ground_station` still say
"No documentation found.". In the diff, the new `relativistic_station_motion` follows the standard
style as an example.

---

## 4. Tests

- `tests/architecture/*.py` (3 files, 236 lines) **are** collected by the bare `pytest` call in
  `.github/workflows/build_and_test.yml`, since there is no `testpaths` setting. They check source
  hashes and include/import rules, which is unusual for tudatpy. Keeping them is a maintainer choice;
  the Desai–Sibois hash test becomes moot with S3. They are kept unchanged in the diff.
- About a third of the new native test code exists to pin the boundary-image semantics (S1), the
  legacy-interpolator rates (S2), the policy enums (S7) and the derivative-guard behaviour (S4).
  These tests go away with those simplifications. Keep the independent IERS/SOFA/HARDISP fixtures:
  they are the valuable part.
- Fixture generators (Python, Fortran and C drivers, `sources.json`, licence copies) under
  `tests/test_tudat/src/...` are fine for traceability. Consider keeping generators in one
  `tests/test_tudat/reference_files/<topic>/` directory per topic instead of mixing them into the
  `src` test tree.

---

## 5. Code style (tudat conventions and Dominic's requests)

- Many one-line `if( … ) throw …;` and `for( … ) if( … ) …;` without braces, e.g. in
  `earthOrientationCalculator.cpp`, `earthOrientationDataEvaluator.cpp`, `groundStationState.cpp`,
  `stationHistoryTrajectory.cpp` and `oceanTideEarthDeformation.cpp`. Dominic asks for braces on
  separate lines, and the rest of tudat uses them.
- New files use a one-line header `/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */`
  or none at all, instead of the standard tudat licence block.
- Several functions have no header or implementation comment (e.g. `EarthOrientationDataEvaluator::stencil`,
  `validStencil`, most of the `stationHistoryTime.cpp` helpers). Many comments describe review
  history ("Retain…", "Preserve the existing…", "legacy…") rather than what the code computes.
- Lambda-heavy, dense expressions (e.g. `calculateRotationRateFromItrsToGcrs( const EarthOrientationAnglesAndRates& )`,
  `getTaiTime`) would benefit from named intermediate quantities with units in the comments, as
  elsewhere in tudat.

---

## 6. Contents of `REVIEW-PROPOSED-CHANGES.diff` (32 files, +169 / −644)

1. Python: remove all `read_*_text` bindings and their docs entries (6 binding files, 5 `.rst` files).
2. `createDefaultTimeConverter`: optional settings, so the shared converter is used only when
   settings are omitted (C6).
3. `ShortPeriodEarthOrientationCorrectionCalculator` matrix constructor: remove the unused
   interpolation argument.
4. `EarthOrientationAnglesCalculator::calculateRates`: remove the legacy custom-interpolator
   finite-difference path and `requireLegacyStencilInSingleInterval`; rates now require EOP data
   (S2, first step).
5. `PrecessionNutationCalculator::getNominalPositionAndRate`: one 32 s stencil of the direct SOFA
   model; remove the interpolant-branch derivative reconstruction, the step setter/getter, the
   interpolator getter and the catch/restore block (S4). Default results are bit-identical.
6. Remove the three dead templated `calculateRotationRateFromItrsToGcrs` overloads (C8).
7. Remove the unused file constructor of `OceanPoleTideDeformation`.
8. Relativistic station motion: remove the legacy constructor and members; one settings
   constructor and factory; one Python function (S7).
9. Merge the "with bodies" overloads of `createGroundStation` into default arguments, with the
    Python `add_ground_station` overloads kept (S7).
10. `EarthOrientationDataEvaluator::evaluate`: `std::array` scratch storage and a generic
    scalar/vector difference-from-anchor sum; arithmetic unchanged (S11).
11. Tests adapted to items 4, 5 and 8:
    - `unitTestEarthRotationDerivatives.cpp`: the nominal-CIP test now checks against independent
      SOFA, and that interpolated positions keep the direct rate; the legacy test checks that the
      rotation is available but the rate throws;
    - `unitTestItrsToGcrsRotationModel.cpp`: the derivative mapping test uses the standard EOP
      calculator;
    - `unitTestRelativisticStationPosition.cpp` and `unitTestGroundStationState.cpp`: new settings
      API; the legacy-constructor lines are removed.

**Not in the diff (described above):** S1, the rest of S2, S3, the full S4 helper, S5, S6, the
remaining S7 rows, S8–S10, most of S11, C1–C5 and the exposure redesign. These are larger API or
structural changes that need the decision of whoever integrates the PR.

## 7. Agreement and disagreement with `MODIFICATIONS-DONE.md`

- **Agree:** keep `Time` precision through clock conversion; interpolate UT1−TAI on continuous TAI
  and differentiate the chosen interpolant; position availability is independent of rate
  availability; keep independent reference fixtures; keep frame/origin/tide-system declarations
  explicit; readers stay free of dynamics dependencies.
- **Disagree (with argument):**
  - "do not add coverage epsilons" (S1): a 1 µs edge tolerance cannot change any result by more
    than ~1e-14 of the quantity and removes the largest source of duplicate code paths;
  - "the generic numerical derivative is not a safe drop-in" (S4): true for the existing
    weighted-sum helper only in principle; the rounding involved is <1e-17 per second of step for
    these quantities, and a paired-difference `Time` overload removes even that concern;
  - keeping split source coordinates and custom time callbacks for station histories (S5):
    the numerical effect is ≤1e-7 m, far below SINEX coordinate uncertainty;
  - keeping arbitrary ordering of frame and reference-point wrappers (S6): the physically
    meaningful order is fixed.
- **New findings not in the coordinator's list:** C1 (J2000 meaning change), C3 (short-term
  interpolation argument now throws), C4 (three hidden TDB re-wrappings), C6/C8, S3
  (Desai–Sibois via resource files), S9 (duplicate grid interpolators), S10, the analytic
  ocean-loading rate, the default atmospheric endpoint policy rejecting the official grid (S7),
  and the per-record provenance (S8).
