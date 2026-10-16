# High-accuracy satellite and lunar laser ranging

Updated: 2026-10-16 13:05 UTC.

## Current authorized task

Latest instruction: arrange the three plotted observing series in separate
vertical columns, with the five reflectors as rows. The columns are APOLLO,
Grasse 532.1 nm and Grasse 1064.2 nm: two physical stations and three optical
series. Complete: Sol HIGH worker `filtered_plots_sol` updated the plotting
helper and regenerated all six PNG/PDF pairs with unchanged residuals,
exclusions, and common time and residual axes. Column counts are 311, 92 and
1,174 respectively. Panel counts, date ranges and full/detail limits were
checked, and representative figures were visually inspected. Integrated as
local commit `ed2a720b2`; no push, model solve or native build. The links below
now open the three-column figures at the same saved paths.

Previous instruction: make residual figures for the outlier-excluded annual
selection below. Complete: Sol HIGH worker `filtered_plots_sol` generated
unconnected five-reflector plots for all three models. Six PNG/PDF pairs are
saved under `year_2017/outlier_filtered/<case>/`, beside the filtered table.
All use the same date span and retain 311 APOLLO plus 1,266 Grasse observations.
The common +/-250 mm axes include every retained residual and supplied error
bar; separate +/-100 mm details label all points outside their display window
(291 DE440, 122 INPOP with glass, 2 INPOP without separate glass).
The files were visually checked, and the README records exact exclusions,
counts, input hashes and a replay snippet. Integrated as local commit `635a5a635`.
No model rerun, native build, further exclusion or push was performed. Original
figures and raw residual files remain unchanged.

Previous instruction: remove outliers from the annual statistics table. Complete:
the same two isolated Grasse L1 observations are omitted from all three cases,
`luna17_201706.npt:58` (June 13) and `:154` (June 22). Their residuals exceed
0.3 m under every model. This is a derived-table selection; all original residual
files and figures remain unchanged, and no systematic model pattern is clipped.
The earlier March exclusion is already absent and is not counted again.

| Current filtered model | APOLLO N | Mean / RMS / P95 absolute (mm) | Grasse N | Mean / RMS / P95 absolute (mm) |
| --- | ---: | ---: | ---: | ---: |
| DE440 | 311 | +15.626 / 72.395 / 146.130 | 1266 | -10.382 / 75.992 / 150.192 |
| INPOP17a, separate glass | 311 | -51.842 / 56.255 / 89.357 | 1266 | -62.619 / 67.296 / 107.297 |
| INPOP17a, no separate glass | 311 | -8.220 / 22.986 / 46.205 | 1266 | -13.543 / 26.286 / 50.078 |

The table, exact exclusions, input hashes and reproduction snippet are saved at
`tudatpy-laser-ranging/analysis/laser_ranging/results/year_2017/outlier_filtered/`.
Sol HIGH worker `slr_sol` performed the CSV-only calculation; integrated in local
commit `3b9362cfb`. No model reruns or builds were needed. This small table update
has not been pushed. The previous raw results below retain the two observations.

Previous instruction: remake the LLR figures using unconnected points, give a
results table, and extend coverage to one year. The declared annual interval is
1 January 2017 inclusive through 1 January 2018 exclusive, for the same APOLLO
and Grasse sources, all five reflectors, and all three ephemeris/optics cases.
Freeze and audit annual source coverage before residual evaluation. No dynamics
or nuisance parameters are estimated. The following user-requested exclusion
is retained. Two additional source records with impossible lunar flight times
are now separately classified as invalid, as described below.

Dominic explicitly requests removal of the single Grasse
Apollo 15 outlier at 2017-03-09T18:34:15.297622200200 UTC,
`apollo15_201703.npt` line 171. Exclude this source identity only, record the
decision, and regenerate current statistics/figures from existing computed rows.
Current March-April results contain 305 points: 69 APOLLO and 236 Grasse. Original
306-point packages remain historical evidence; their statements below about the
retained outlier describe those earlier runs. No observation rerun is needed.
The physical cause remains unestablished; this is an explicit user selection.
The current CSV rows are copied unchanged from the original calculations.
Marker-only PNG/PDF figures and the table are committed in `fb4958915`, under
`analysis/laser_ranging/results/current_305_20170301_20170430/`.

| Model | APOLLO mean / RMS (mm), 69 | Grasse mean / RMS (mm), 236 |
| --- | ---: | ---: |
| DE440 | +1.830 / 58.488 | -12.959 / 70.101 |
| INPOP17a, separate glass | -47.439 / 49.408 | -60.826 / 63.920 |
| INPOP17a, no separate glass | -3.363 / 12.341 | -12.807 / 22.095 |

The optical reference convention of the published INPOP reflector coordinates
remains unresolved. Both predeclared choices are retained; smaller residuals
alone do not establish the correct convention.

Annual source audit: all 60 official monthly files retrieved; March-April files
are byte-identical to the frozen originals. There are 311 APOLLO points on 42
nights and 1,271 Grasse source points on 114 nights (1,268 eligible after the
requested exclusion and two invalid flight times). The same products contain
33 Matera points outside the predeclared station scope. Native reader checks
confirm Grasse correction flags, transmit epochs and recorded wavelengths.
Three June 22 Grasse points have stale measured weather: endpoint ages about
1,071, 1,889 and 6,601 s. Official daily files repeat the same last measurement.
They remain in the annual results with explicit source-identity exceptions and
weather-limited flags, using the last actual measured tuple. This is a declared
analysis assumption with unquantified meteorological error, not evidence that
weather remained constant. The ordinary 1,000 s guard remains for all other
points. No later weather is invented and no extra observation is rejected.
The source also repeats one complete Grasse NP record, at August 28 19:50:52,
in `luna21_201708.npt` lines 102/105. Both source occurrences are retained and
documented; they are not independent measurements. The loader now distinguishes
source lines even when the printed epoch and round-trip time coincide.
The source audit subsequently identified two nonlunar flight times:
`apollo15_201702.npt` lines 168/211 report 24.976/17.271 ns, or 3.74/2.59 m
half-ranges. Daily files repeat these exact records. Their cause is unknown;
neither is a usable lunar return. They are excluded by exact source identity,
with separate invalid-source reasons, independent of residuals. All raw provider
bytes and 1,582 source rows remain in the audit. Final eligibility is 1,579.
The annual manifest is now SHA-256
`4140847ee51866fed8d3d2a92d6ead4d03113a7e67885adcbfc52958d9174446`.
All three first/last controls pass: four distinct expected boundary observations
per model, with valid clock conversions, table coverage and finite predictions.
The first full runs stopped at the unchanged atmospheric elevation guard at one
of the invalid source records; their logs are preserved. No final CSV was
produced. The complete source geometry check finds the two invalid records
below -18 degrees; every other point is above 12.5607 degrees at both terminals.
The Moon-centre screening approximation has ample margin above the 3-degree
atmospheric domain, including reflector offset and light-time effects.
All three annual cases now finish successfully: 1,579 rows each. Independent
verification confirms exact source identities, exclusions, duplicate records,
weather flags and coverage. All 915 overlapping computed ranges (305 per case)
are bit-identical to the current March-April controls. Executed source blobs
agree with root `5e4934584`; binary and model formulas are unchanged. The three
parallel runs took 236.6 s wall time in total; no native build was needed.

| Annual model | APOLLO mean / RMS / P95 absolute (mm), 311 | Grasse mean / RMS / P95 absolute (mm), 1268 |
| --- | ---: | ---: |
| DE440 | +15.626 / 72.395 / 146.130 | +77.172 / 3102.942 / 151.711 |
| INPOP17a, separate glass | -51.842 / 56.255 / 89.357 | +24.734 / 3097.276 / 107.421 |
| INPOP17a, no separate glass | -8.220 / 22.986 / 46.205 | +73.833 / 3098.441 / 50.259 |

The annual extension reveals another retained discrepancy: Grasse L1 on June 22
at 09:47:06.763983741400 UTC, `luna17_201706.npt:154`, is +110.328 m in the
no-separate-glass case. Its cause is not established. It is one of the three
weather-limited records, but this does not attribute the discrepancy to weather.
No additional residual cut is applied; the full-range figures retain it and
separate +/-200 mm figures expose the central structure. Actual eligible
coverage is APOLLO January 6-December 31 (42 dates), Grasse January 6-December 21
(113 dates). The raw 114 Grasse dates include the invalid-only February 20.
Annual figures, descriptive temporal comparison and packaging are complete.
The result package is
`/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging/analysis/laser_ranging/results/year_2017/`.
Its README links all three raw CSVs, the statistics table, and unconnected
full-span and +/-200 mm figures in PNG/PDF. The adjacent `year_2017_temporal/`
package contains monthly medians, dated nightly extrema and station comparisons.
The no-separate-glass A11 medians show a common June-July trough: APOLLO
-42.8/-38.4 mm, Grasse -48.5/-48.6 mm. Shared centimetre structure remains.
The 64 nearby same-reflector night pairs are 2.13-11.46 hours apart, so they
are not simultaneous station-bias measurements. No trend was fitted or removed.
All annual and temporal artifact/source hashes have been checked after integration.
The annual milestone, scripts and report are committed and pushed as
`c5890681f` on `origin/feature/high-accuracy-laser-ranging`. The remote head is
verified and the root working tree is clean. No native build or additional
observation rerun followed the successful annual verification.

Rerun the two frozen LAGEOS-1 weeks with the highest-fidelity compatible Earth
rotation and station deformation available. Keep the external final postfit
orbits, observations, published calibration corrections and selection cuts fixed.
No propagation, estimation, fitted offset, demeaning or residual clipping.
Unrelated CI work remains paused. Dominic additionally requested a deeper LLR
investigation: Hannover station solutions, kernel loading order and representation
accuracy, consistent ephemeris comparisons, Shapiro epochs, and relativistic
frame/time conventions. Preserve the frozen eleven APOLLO and two Matera points.
If these controls leave the reflector-dependent discrepancy unresolved, extend
to a fixed one-to-two-month interval and compare stations targeting the same
reflector. The interval is frozen to March–April 2017, with 69 APOLLO and 237
Grasse normal points selected from source coverage before computing residuals.
Dominic additionally requires saved LLR residual-versus-time figures and a trend
assessment. The long-arc output will include raw per-reflector, cross-station
plots, exportable images and the plotting script; no trend is removed from data.
Dominic has now authorized further error investigations and safe, measured
runtime improvements. Preserve the frozen observations and numerical accuracy.

- Coordinator: completed annual source/model/output review, integrated all
  worker commits, updated the report and pushed the annual milestone. Complete.
- slr_sol (Sol HIGH), `tudatpy-laser-data`: completed the full 2017 source audit;
  independent verification passes for output identities, the three weather flags
  and exact agreement with every overlapping March-April computed range. Idle.
- llr_sol (Sol HIGH), `tudatpy-laser-llr`: completed the exact single-point
  exclusion, unconnected March-April figures and annual temporal/cross-station
  summary without fitting an adjustment. Idle.
- compile_sol (Sol HIGH), `tudatpy-laser-signature`: completed the annual driver,
  all three fixed-ephemeris runs and annual figure/table package. Idle.
  Sole compiler owner; zero compiler jobs were needed for this extension.

Current extension checklist:

- [x] Use geocentric TT only for linear station drift. All 11,934 CSV fields
  remain byte-identical; solver time falls 19.631 to 12.778 s (34.9%), and total
  driver time 27.066 to 18.994 s (29.8%). Maximum coordinate effect is below 2e-15 m.
- [x] Check absolute launch time conversion, including source epoch precision.
  Three events agree with independent ERFA/SOFA within 2.5 ps; no missing or
  duplicated TDB-minus-TT conversion is found.
- [x] Check actual lunar radius, Love numbers, degree normalization and once-only
  deformation. Independent/native vectors agree within 0.145 nm; scaled radius
  and shape controls agree. No residual correction follows from these checks.
- [x] Check representative NP epochs and apparatus correction conventions.
  Sources support the current launch tags and already-applied system delays.
  The daily EDC archive independently repeats the same retained +5.43 m record.
- [x] Evaluate a finite moving-body Shapiro expression at saved actual endpoints.
  Exact first-PM constant-velocity terms change individual leg lengths by at most
  2.282 micrometres; the opposite uplink/downlink terms nearly cancel. This is a
  fixed-endpoint control, not a newly iterated observable. A diagnostic originally
  omitted Earth's terminal spatial map; this is corrected in `991289322`.
  Native reconstruction now agrees within 0.0271 nm (previously 3.23 nm), and
  the moving-body differences change by less than 3.7e-14 m. Observation solves
  were unaffected; independent quadrature and zero-velocity checks close.
- [x] Compare SOFA against a separately integrated time ephemeris. Three launch
  epochs agree with the INPOP17a time SPK within 1.393 ns; manual Chebyshev and
  SPICE evaluations agree within 4.6e-18 s. No clock swap is warranted.
- [x] Run a coherent INPOP17a control after resolving documented table misprints
  through the cited INPOP13c bases and the published 17a coordinate corrections.
  Source review is complete and implementation is active. The fitted reflector
  optical reference point remains unspecified in the retrieved INPOP sources:
  two explicitly declared glass conventions will be shown without selecting one
  by residual size. Actual reflector radii and published degree-2 tides are used;
  the fixed 1738 km radius differs by several millimetres in this reduction.
  The first 12 predeclared points show materially different temporal and reflector
  patterns. Published kernel constants strongly identify the G+IR dynamics, though
  the archive README does not explicitly name that fit. Both full-arc runs complete.
- [ ] Record which remaining reduction terms are unavailable or only bounded,
  without claiming the unexplained centimetre residuals have been resolved.

New INPOP17a sensitivity results, all original points retained and no fitted
parameters or offset correction:

| Explicit optical convention | APOLLO mean / RMS (mm), 69 points | Grasse mean / RMS (mm), 237 points |
| --- | ---: | ---: |
| JPL n*h per leg | -47.439 / 49.408 | -37.637 / 358.749 |
| No separate glass term | -3.363 / 12.341 | +10.348 / 356.316 |

The second line is a predeclared sensitivity, not a convention selected by
residual size. INPOP's optical reference point remains unspecified in the
retrieved sources. Its Grasse outlier remains +5.475 m; no point is rejected.
Raw median absolute / 95th-percentile absolute residuals in the no-separate-glass
case are APOLLO 8.212 / 22.850 mm and Grasse 12.529 / 37.777 mm, using every point.
The March 22 APOLLO spread between five reflector nightly means falls from
139.25 mm (DE440) to 17.06 mm in this INPOP sensitivity, but a shared curved
temporal pattern remains. Intrinsic APOLLO measurement-noise residuals are not
achieved. INPOP Table 10 also lists a Grasse fitted bias covering this arc,
equivalent to 12.3 mm one-way; it is documented, never applied to these results.

Latest native change is `a366f436c`, the recorded-wavelength guard fix. The build is
`/home/dominic/Tudat/laser-ranging-build`, kernel SHA-256
`cab49e738cfcf0ea5f251cb140ab8919d5372dbc3e5ba09314d0d5be4eedfcef`.
The affected atmosphere suite passes all 15 cases / 9880 assertions; Python
accepts the recorded 532.1/1064.2 nm values unchanged. Earlier precision/clock
changes passed the native N-way suite (five cases / 1023 assertions) and
all three Python environment tests. The new binding test initially made
an unsupported 20 ns derivative comparison: independent SOFA evaluation showed
inherited ERA quantization equivalent to 0.12 micrometres of station motion.
That new assertion was removed; the distinct-Time, float and inverse-method
checks remain. Existing numerical test tolerances and production behavior are
unchanged. The numerical floor and original failure are preserved in the logs.

The confirmed LLR frame mismatch is corrected locally in the analysis: Tudat's
`J2000` Earth rotation applies the dynamical-frame bias, whereas modern SPICE
products label ICRF axes `J2000`. Explicit ICRF labels and an identity SPICE
alias keep the whole environment consistent. No global frame semantics changed.

Verified SLR/short-arc LLR milestone was pushed at `cd36e63a4`. The completed
long-arc diagnostic, wavelength fix, source/numerical reviews, raw figures,
runtime proposals and documentation are now pushed to
`origin/feature/high-accuracy-laser-ranging` at `3dc1574d9`. Local diff checks
pass; all 306 rows and plotting source/input hashes agree. The newly authorized
drift-clock optimization is committed locally at `596deb098`; its output is
unchanged. LLR residual accuracy remains unresolved.
Short-arc raw time figures are saved under
`analysis/laser_ranging/results/short_de440_time_plots/` in PNG and PDF. The
shared plotting script is `analysis/laser_ranging/plot_llr_residuals.py`.

All 69 APOLLO normal points over nine nights are now saved under
`analysis/laser_ranging/results/apollo_long_arc_20170301_20170430/`, including
raw CSV, geometry, provenance and PNG/PDF time figures. Mean is +1.830 mm,
RMS about zero 58.488 mm, and scatter 58.460 mm. The March 9 control is
reproduced exactly. Runtime is 5.332 seconds, including 0.633 seconds setup
and 4.244 seconds solving. The variation is nonmonotonic: several reflectors
rise toward March 22, fall toward early April and then recover. Apollo 15 has
159.364 mm raw peak-to-peak variation despite a near-zero descriptive linear
slope. Sparse nights do not establish a period. No trend or mean is removed.

The first combined run stopped before solving any Grasse points because its
recorded 1064.2 nm laser exceeds the native guard's nominal 1064.0 nm ceiling.
The fix accepts a half-nanometre reporting interval around that nominal line,
without changing the supplied wavelength or optical formula. This is an
explicit acceptance convention, not a newly published validation endpoint.
A numeric control from the unmodified official IERS FCUL routine was added
at the recorded wavelength, keeping the existing 1e-9 m test bound. The full
affected atmosphere suite passes all 15 cases / 9880 assertions. The Python
kernel rebuild and actual-wavelength checks pass; the complete run then succeeds.

All 306 frozen normal points now complete in 26.6935 seconds: 5.8127 seconds
setup and 19.3589 seconds solving, with no runtime optimization. Every APOLLO
row reproduces the earlier control exactly. Grasse's 237 points have raw mean
+10.007 mm and RMS 359.583 mm. One retained +5.429946840 m A15 point accounts
for 96.2% of its squared residuals. Its printed source fields and interpretation
are verified, but its cause is not established. The other 13 points in the same
session are −71.3 to −36.7 mm. Re-solving this point and its two neighbours
with 15-second instead of 30-second Time tables changes computed half-ranges
by at most 10.133 micrometres, excluding that interpolation as the metre-scale
cause. No point is rejected or reweighted.

The complete raw package is
`analysis/laser_ranging/results/paired_long_arc_20170301_20170430/`; combined
PNG/PDF figures are in `analysis/laser_ranging/results/paired_long_arc_time_plots/`.
Both full range and a clearly labelled ±200 mm detail retain the original data
and statistics. The detail explicitly labels its one outside-window point.
Seventeen nearby-night pairs, with actual 4.675–11.458 hour gaps, show similar
curved changes and reflector ordering at both stations. Their raw station
differences span −31.684 to +19.386 mm; they are not simultaneous observations.
There is no demonstrated monthly period or unique physical attribution.

A final original-ASCII-to-native lunar orientation control at two frozen
epochs agrees within 0.75 micrometres at the surface. It rules out the tested
fixed lunar-frame composition error. The unresolved physical reduction/fit
parameters and the isolated Grasse discrepancy remain documented; intrinsic
measurement-noise residuals have not been achieved.

Runtime proposals are saved in
`analysis/laser_ranging/results/llr_runtime_audit/README.rst`. The small profile
attributes about 28% of elapsed time to repeated topocentric clock conversions;
simulator creation is below 0.2%. Earth rotation derivatives and repeated station
evaluations are other candidates. A simpler first proposal is geocentric TT
only for slow linear station drift: the measured coordinate bound on this arc
is below 2e-15 m, while observation clocks and simultaneity remain exact. Full
rotation/displacement tabulation needs precise-Time and rate support before
it can safely replace the direct calculation. The drift-only proposal has now
been implemented and verified; other proposals remain unimplemented.

Completed SLR evidence is integrated through `5cf270065`. Both weeks use kernel
SHA-256 `10f8307ba799847f11d9614d2c3f3da3db23962c46a739d73f536f9dc31ae5d5`,
from native source `b1aa77250` in `/home/dominic/Tudat/laser-ranging-build`.
The cubic arithmetic fix passes 11 focused cases / 4052 assertions and both
full-week runs converge with the original light-time tolerance.

| Week | Retained / raw points | No-loading mean / RMS (mm) | Geodetic HARPOS mean / RMS (mm) |
| --- | ---: | ---: | ---: |
| 31 March–6 April 2024 | 23 / 29 | −0.214 / 8.769 | −0.057 / 7.221 |
| 7–13 April 2024 | 131 / 134 | +5.258 / 18.752 | +7.659 / 17.435 |

The requested highest compatible tidal profile now passes both frozen weeks on
the new kernel. Every new run uses ICRF and precise-Time rotation. The three
profiles retain identical observations, orbit/calibration inputs and cuts:

| Week | Standard EOP + solid/pole/GOT mean / RMS (mm) | Desai–Sibois/four-point mean / RMS (mm) | Plus ocean-pole/RP03 atmospheric tides mean / RMS (mm) |
| --- | ---: | ---: | ---: |
| 31 March–6 April | −0.060 / 7.221 | −0.060 / 7.222 | −0.405 / 7.142 |
| 7–13 April | +7.659 / 17.435 | +7.659 / 17.435 | +7.536 / 17.456 |

The two additional tides change individual residuals by at most 1.010 mm.
No model was chosen by residual size. Non-tidal histories remain excluded until
their basis and removed mean are demonstrably compatible with the station datum.

LLR remains diagnostic, with fixed DE440 ephemerides and no fitted parameters:

- APOLLO, 9 March 2017, 01:19:15–02:10:30 UTC: 51 min 15 s,
  eleven normal points over four reflectors. Published INPOP station coordinates
  with DE440 and corrected ICRF axes give mean −25.500 mm / RMS 49.900 mm
  without loading (previous mismatched-frame RMS 89.291 mm). Reflector means
  remain A11 −85.690, A14 +14.505, A15 −51.002 and L1 +25.689 mm. The station
  solution is not matched to DE440.
- Matera, 21 August 2022, 02:46:45.004–03:48:26.604 UTC: 1 h 1 min 41.600 s,
  two Apollo 15 normal points: corrected mean −65.200 mm / RMS 65.228 mm,
  versus the mismatched-frame RMS 152.683 mm.
- The 69-point APOLLO arc now spans March 2–April 21 within the frozen two-month
  interval. Independent lunar-tide, optical and clock checks have not explained
  its remaining centimetre-scale discrepancy. Intrinsic measurement-noise
  accuracy has not been reached; the completed Grasse comparison shows a
  similar curved signature and reflector ordering, with one isolated large point.
- Original DE440 ASCII versus binary coefficient checks show no centimetre
  representation loss: lunar libration differences are below one micrometre at
  the reflector, and translation coefficient differences are tens of micrometres.
  Direct-SPICE double epoch rounding can reach 1.75 mm in absolute body position;
  the existing precise-Time tables avoid that final-event truncation.
- Actual-event Shapiro evaluation agrees with an independent calculation within
  3.26 nm. Earth is evaluated at Earth transmit/receive, Moon at reflection, Sun
  at each leg midpoint. Alternative body epochs change one-way range by at most
  49.6 micrometres here; other planets add up to 1.82 mm. No centimetre Shapiro
  epoch error has been identified.
- All six controlled station/ephemeris runs pass on the new kernel. The DE440
  controls reproduce the previous corrected-frame CSVs exactly. Published Hofmann
  station coordinates give APOLLO mean −3.119 mm / RMS 42.083 mm and Matera
  mean −24.635 mm / RMS 25.318 mm. They are an independent station sensitivity,
  not a recovered fit consistent with DE440; no offsets are fitted here.
- Coherent DE421 SPK/PCK/MER and source constants, with the original stations,
  give APOLLO RMS 282.428 mm and Matera RMS 306.031 mm. The published MER
  reflector coordinates are frame-compatible but not recovered original DE421
  fit parameters. DE440 remains the predefined primary case.
- The APOLLO reflector spread persists under the station change, justifying the
  conditional longer arc. March–April 2017 has 69 frozen APOLLO points on nine
  nights and 237 Grasse points in 70 CRD sessions, covering all five reflectors.
  The correct Grasse telescope is CDP7845 / DOMES10002S002. Its INPOP17a G+IR
  station row and calibration/weather handling are verified. The complete run
  keeps all frozen data, with no wavelength rounding or residual-based selection.

The dated sections below retain the previous checkpoints; the current task and
results above supersede their pending-build and failed-convergence statements.

## Scope and working rules

Use fixed published postfit satellite orbits and lunar ephemerides with actual normal points.
Improve the physical observation models where necessary. Do not remove fitted offsets or trends
from residuals, tune tolerances to results, or propagate/estimate an orbit unless essential.
Keep changes local to laser ranging; reuse the existing station, time, rotation, atmosphere and
light-time machinery. SLR-specific interfaces may be refactored without compatibility wrappers.

Sol HIGH workers implement, compile and handle small debugging tasks. Coordinator designs, reviews
and handles major numerical debugging. Workers use separate worktrees once implementation starts.
At most eight compilation jobs in total, focused native tests, Python builds for the exposed feature
and requested analysis script. No execution of Orekit/GROOPS unless a last-resort stepwise comparison
is needed. Source inspection is allowed and does not replace independent physical references.

Commit coherent changes regularly. Dominic explicitly authorizes pushing substantive milestones,
including the Python analysis script. Only the coordinator pushes; no AI attribution in commits.

## Branch and agents

- New branch: `feature/high-accuracy-laser-ranging`.
- New worktree: `/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging`.
- Starting gravity branch: `fix/gravity_field_variations`, `322800bd9d6703ee3501f9c9c7341dd3b1302c7c`.
- Latest rotation/deformation branch: `fix/earth_deformation_rotation`,
  `638167cd05a8dab2a4d9b684a266bf4cb8f0f48e`.
- compile_sol (Sol HIGH): sole owner of the focused rebuild at `b1aa77250`, at most eight compilation jobs. Cubic and N-way suites pass (11 cases, 4052 assertions); Python kernel still compiling. Independent LLR clock audit complete.
- slr_sol (Sol HIGH): cubic arithmetic fix and regression integrated; prepares frozen SLR reruns and pins a second week (7–13 April 2024) before residual calculation.
- llr_sol (Sol HIGH): reflector link/cache review found no wrong station reuse; prepares minimal CI test fixes in an isolated rotation worktree, then resumes the four fresh-kernel LLR runs.
- Working research/data: `/home/dominic/Tudat/laser-ranging-work/`.

## Checklist

- [x] Integrate latest gravity and rotation/deformation work. [DONE — merge `ebe3115ec`]
- [x] Select published SLR cases and distinguish published precision metrics from exact residual references. [DONE]
- [x] Select published LLR cases and distinguish observation uncertainty from ephemeris error. [DONE]
- [x] Record immutable data URLs, hashes, units, frames, time scales and correction flags. [DONE — manifests beside both scripts]
- [x] Audit existing Tudat, Orekit and GROOPS laser-ranging model definitions. [STATIC AUDIT RECORDED — no external solver execution]
- [x] Define observation epochs, elapsed station-clock time and one-way residual normalization. [IMPLEMENTED — explicit station simultaneity and one native TAI interval conversion]
- [x] Preserve required normal-point/configuration/calibration metadata in the SLR reader. [NATIVE AND PYTHON CHECKS PASS]
- [x] Load fixed postfit satellite orbits and check interpolation accuracy. [DONE — native reader passes; independent node holdouts RMS 1.4–1.5 mm for millimetre-printed SP3 positions]
- [x] Construct station histories and deformation in the frame used by each orbit product. [SLRF2020 XYZ/velocity and solid/pole tides active; explicit ocean-loading experiment under numerical review]
- [x] Check optical atmosphere, relativistic delay and station clock/frame transformations. [NATIVE REFERENCES PASS — precise-time refinement continues]
- [x] Add a local target-signature model with observer-dependent reflection geometry and explicit
  optical weighting; prevent duplicate centre-of-mass corrections. [DONE — four native cases and Python surface checks pass]
- [x] Obtain LAGEOS raw residuals with complete correction accounting and published calibration only. [23 POINTS — RMS 8.843 mm, mean −0.474 mm; six predefined exclusions, no fitting/clipping; ocean loading omitted]
- [x] Assemble matched lunar translation/rotation, reflector coordinates and lunar deformation. [ALL ELEVEN FROZEN POINTS RUN]
- [ ] Obtain fixed-ephemeris LLR residuals and quantify any model floor above measurement noise. [CORRECTED ICRF: APOLLO RMS 4.990 cm; Matera RMS 6.523 cm — station/ephemeris controls running]
- [x] Provide reproducible Python analysis scripts and numerical result files. [PORTABLE SLR, APOLLO AND MATERA RAW RESULTS COMMITTED LOCALLY]
- [x] Validate focused native calculations and Python interfaces; preserve numerical accuracy. [SIX NATIVE SUITES: 34 CASES / 11154 ASSERTIONS; 59 PYTHON TESTS; ALL 87 AUTOSUMMARY NAMES PASS]
- [ ] Review for duplicate code, generic invasive changes, redundant helpers and excessive tests. [ONGOING — native Time ephemeris reused; pinned loading reader simplified]
- [x] Push a substantive milestone with analysis script and evidence. [PUSHED — `9ff7bd7f8` on `origin/feature/high-accuracy-laser-ranging`; scripts, raw SLR results and validation report included]

## Current residual checkpoint

- Frozen SLR repeat: identical CSV bytes, 23 points, mean −0.474 mm and RMS 8.843 mm; six fixed exclusions unchanged.
- Precise-Time SLR callbacks: same 23 points, mean −0.215 mm and RMS 8.768 mm. Changing body-table spacing from 60 to 30 s changes computed half-range by at most 9.585 micrometres. Ocean loading omitted in this baseline.
- LLR survey station, cubic Time tables: eleven points, mean +219.550 mm and RMS 223.164 mm.
- Published INPOP17a G+IR APOLLO position/velocity sensitivity: eleven points, mean +79.877 mm and RMS 89.294 mm. This is a station-only comparison with DE440, not a matched fitted solution; no fitting or residual cuts. Matched-kernel control is exact and portable results are committed locally.
- Independent Matera, 2022-08-21: both Apollo 15 points retained, residuals +177.724 and +122.627 mm, mean +150.176 mm and RMS 152.682 mm. Native SLRF2020 history; no applicable official PSD/discontinuity record; no invented NP uncertainty. Ocean loading omitted.
- Historical radial-Up APOLLO loading controls give RMS 235.149 mm (survey) and 100.350 mm (INPOP station). These are preserved unchanged; provider source now supports geodetic Up, independently of residuals.
- The SLR loading run stopped at a real receiver-leg two-cycle: epoch separation 31.752 fs, equivalent to 9.519 micrometres, above the unchanged 7.267-micrometre convergence criterion. An independent 80-digit calculation proves that the old cubic weighted sum introduced artificial 30.518/7.629-micrometre jumps. The same polynomial evaluated from the nearer endpoint correctly rounds both positions. The fix and regression pass the complete existing cubic suite; native N-way tests also pass. Kernel rebuild and affected residual runs remain; no tolerance changes or complete loaded SLR result yet.
- A second week, 7–13 April 2024, is frozen before computing residuals: 134 raw normal points across 17 Graz/Herstmonceux sessions. Same orbit/frame/calibration conventions and predeclared cuts as the first week.
- Precise-Time binding build passes both focused environment tests and all 57 environment autosummary names. Kernel SHA-256: `bd8d893574bd29281188934ff5a62647e4c08494eb4bfb846c17cfc06a062968`.
- Independent lunar degree-2/3 tide calculation agrees with the native model within 0.15 nanometres; complete reflector correction within 0.27 nanometres. No tide implementation error or duplicate constant tide was found.
- Independent LLR clock calculation agrees within 0.094 mm, with geocentric agreement within 0.397 micrometres. The existing 1 ps inverse-conversion tolerance explains nearly all the small discrepancy. No centimetre-scale clock error or duplicate station term was found; converter unchanged.
- Actual-input optical audit reproduces all 48 APOLLO/Matera leg delays to floating-point precision, with correct geodetic geometry and pressure/temperature/humidity units. Matera's documented 2.4 m pressure-sensor height difference could affect range by about 0.8–0.9 mm, but whether the reported pressure is already reduced is unknown; no speculative correction added.
- Current APOLLO CRDv2 release changes the same-epoch ranges by −3.687 to +12.052 mm, shifts one representative epoch by 5 s and adds one normal point. It is archived in quarantine with WGL021 (missing station XYZ for the archive's orbit-fit check). Original frozen observations remain unchanged; per-point differences and source hashes are committed in worker commit `75ca415c4`, pending integration after build provenance is recorded.

## Findings requiring explicit accounting

- Existing CRD structures retain NP times and basic metadata but omit several session correction
  flags. The current range extractor returns one-way c*TOF/2 while native N-way range returns total
  round-trip distance. The new workflow must reconcile them explicitly.
- The full-link correction interface supports target geometry. The implementation corrects each
  leg once using a common weighted response; mean path change and return spread are distinct.
- JPL DE440 reports about 1.3 cm recent LLR RMS. APOLLO reports millimetre measurement accuracy and
  documents ephemeris deficiencies above that scale. Reproducing the former does not establish the latter.
  Sources: https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf and https://arxiv.org/html/2304.11128
- The first SLR case uses the coherent 2024 products; the unresolved 2018 frame-label discrepancy
  is recorded only for the secondary case.
- Existing station spatial relativity does not retime the rotating station to the local terrestrial
  epoch of a barycentric event. This must be accounted for before a millimetre result is claimed.
- JPL's full 2020 LLR reduction model is now cached in `laser-ranging-work/research/JPL-LLR-Model-2020.*`.
  It supplies the lunar spatial scale (LM = 1.4825e-8), surface Love numbers and the reflector
  coordinate convention: published array positions refer to the front face, so internal glass
  path delay is separate from the already-reduced normal-point array profile. Implementation
  adds an explicit optional optical path length to the array model, without guessing a glass law.
  Source: https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/LLR_Model_2020_DR.pdf

## PR 1092 CI (separate from this new branch)

- Cleanup commit `638167cd0` pushed; local working tree clean; both commit hooks passed.
- Build/test run: https://github.com/tudat-team/tudatpy/actions/runs/38797661346
- Last checked: 23:16 UTC. Formatting passed. Linux and macOS compiled successfully and each passed 377/379 native suites. Both fail EarthOrientationCalculator; Linux additionally fails ObservationRoundingErrors, macOS additionally fails ItrsToGcrsRotationModel. Windows remains in progress.
- The first failure compares a TT epoch as though it were UTC when retrieving pole corrections. The matrix and rounding checks are under numerical review; no tolerance changes made.
- Dominic changed polling to once per hour. Next check: 00:17 UTC on 16 October.
- Continue diagnosing and fixing actual CI failures on the rotation/deformation branch.

## Implementation allocation at 21:24 UTC

- desai_impl: `tudatpy-laser-signature`, `work/laser-target-signature`; native scalar CoM/optical-array
  target correction using the existing full-link light-time API. Owns sole compile budget, up to8
  jobs in new `/home/dominic/Tudat/laser-ranging-build` against the clean integration source.
- history_impl: `tudatpy-laser-data`, `work/laser-ranging-data`; CRD session/configuration/calibration
  metadata and ordered precise epochs, plus raw SP3 data input. No compilation allocation.
- rotation_impl: `tudatpy-laser-llr`, `work/laser-ranging-llr`; Mendes-Pavlis optical zenith delay with
  FCULa mapping, and fixed-ephemeris Python LLR analysis. No compilation allocation.
- First SLR case is frozen before residual calculation: LAGEOS1 2024-03-31–2024-04-06, Graz7839 and
  Herstmonceux7840, combined v81 SLRF2020 orbit, matching station/EOP products and DHF2024.02.13 +
  CoM2022.09.15. ASI orbit is a separate comparison. External research notes document conventions.
- First APOLLO case is all11 calibrated observations from2017-03-09 with matched DE440 PA products.
  Precise fitted observatory parameters remain unresolved; initial results must reflect that limitation.
- Optical group-delay support driven by measured meteorology is missing; the existing VMF3o route
  consumes different products. The new optical correction follows the IERS model and will be checked
  against published reference values. No radio-atmosphere substitution.

## Integration and verification at 21:40 UTC

- The three feature implementations and minimal bindings are integrated locally on the laser
  branch through `c71f6f827`: scalar/array target response with explicit internal optical path,
  Mendes-Pavlis/FCULa atmosphere, precise CRD metadata/epochs, raw SP3 input, shared rotation
  time-converter access, and exposure of the existing body-deformation station-motion factory.
- The base dependency build exposed a missing explicit observation-factory include in the existing
  n-way range test. Commit `c71f6f827` adds that include; no numerical behavior or tolerance changed.
- The final native graph builds LaserTargetSignature, LaserRangingData, CrdAndSinexReaders,
  NWayRangeObservationModel and AtmosphereCorrection. One compiler owner, at most eight jobs;
  tests will run serially. Python kernel follows for the actual analysis and exposed interfaces.
- Build directory: `/home/dominic/Tudat/laser-ranging-build`; logs:
  `/home/dominic/Tudat/laser-ranging-work/logs/`. Nothing on the new feature has been pushed yet.
- No residual calculation has completed yet. Source commits and analytical controls are progress,
  not evidence of achieved SLR or LLR accuracy.

## Validation checkpoint at 21:53 UTC

- All five focused native targets compile. The existing range/atmosphere tests required explicit
  factory/collection/shape includes after transitive includes were removed; only includes changed.
- LaserTargetSignature: PASS, four cases / 21 assertions.
- LaserRangingData: FAIL on the archived real SP3 input (incomplete header/samples/EOF diagnostic).
  The data worker owns the parser/input investigation; no numerical tolerance was changed.
- Remaining native suites run serially while the sole compiler worker builds the Python kernel
  with at most eight jobs. Source native content is `d550146ca`, script/docs successor `8b11b118a`.
- Both analysis scripts now live in `analysis/laser_ranging/` on the integration branch.
  Each has a portable provider-URL/hash manifest. Actual residual runs remain pending.
- Static independent-software notes are in
  `/home/dominic/Tudat/laser-ranging-work/research/software-audit/NOTES.md`.
- Public IMLS GOT4.10c harmonics were located for Graz/Herstmonceux; they are a different loading
  realization from ASI's BLQ/HARDISP setup and have not been silently substituted.

## Numerical verification checkpoint at 21:59 UTC

- Native passes retained: LaserTargetSignature 4 cases/21 assertions,
  CrdAndSinexReaders 6/59, NWayRangeObservationModel 5/1019.
- SP3 failure diagnosed: the original combined ILRSA file has all 5040 epochs and
  complete final P/V records, but no literal EOF line. Commit `67c7d49e4` accepts
  physical EOF only after the existing completeness checks; archived bytes and
  numerical assertions are unchanged. Its failed suite will be rebuilt/rerun.
- AtmosphereCorrection passes 14 of 15 cases. Three zenith-delay assertions differ
  from the IERS source header's printed values by up to 3.8 micrometres. Independent
  arithmetic matches our formula; the printed reference implies a seven-metre
  different height. We will execute the unmodified official Fortran routine at its
  stated inputs before accepting any reference correction. The 1 nm tolerance stays.
- The Python kernel is still compiling, with no compiler errors. Afterward, the
  sole compiler owner will execute the small IERS reference and rebuild only the
  changed reader, affected atmosphere test and kernel links.
- A further precision review identified double epoch rounding in direct SPICE/custom
  ephemerides as a possible millimetre BCRS endpoint error. Existing Time-valued
  tabulated ephemerides are being assessed as an analysis-only remedy; no generic
  ephemeris framework or unverified correction has been added.

## Independent reference and runtime preparation at 22:10 UTC

- The initial Python kernel build completed successfully. The only active build is
  the incremental reader/optical-test/kernel graph, owned by desai_impl with at most
  eight jobs. The two failed native suites will run serially after this graph.
- The original, unmodified IERS Fortran at its printed inputs returns total
  1.93522972487606770 m, dry 1.93299597223469832 m, wet
  0.00223375264136944527 m. These agree with the new formula within the unchanged
  1 nm tolerance. The printed header controls are inconsistent; no model formula
  or tolerance was changed to match them.
- Commit `854029719` records those independent controls. The archive now has 830
  unique entries, including all current SLR fixtures plus the licensed original
  Fortran, driver, compiler/command/output and discrepancy note. Every existing
  fixture and new reference was checked byte-for-byte; archive CRC checks pass.
- Both analysis workers have their frozen input selections and scripts ready.
  They will run only after the compiler worker verifies the final kernel import
  path and binary hash. No residual accuracy result is available yet.
- A public GOT4.10c ocean-loading harmonic product is cached for independent
  evaluation. It includes the actual Graz, Herstmonceux and APOLLO sites and
  declares a whole-Earth centre-of-mass origin. Its 17-constituent synthesis differs
  from ASI's BLQ/HARDISP realization; it has not been silently substituted.

## First runtime results and recovery at 22:14 UTC

- All five native suites PASS: 31 cases / 11052 assertions. Only the reader and
  optical suites were rerun after their fixes; earlier unaffected passes retained.
  `logs/native-final-summary.json` identifies all final log/XML artifacts.
- The 59 focused Python tests PASS. Their output is `logs/python-focused.log`.
  Checking every edited autosummary entry against Python remains in progress.
- Verified kernel: `laser-ranging-build/src/tudatpy/kernel.so`, SHA-256
  `c475a3e90a39eed4332442f58c631398cd608f2195d25bc08cd59e3cfd0bd866`,
  built from `854029719`. Runtime provenance is in `logs/final-kernel-provenance.json`.
- All eleven APOLLO points ran without drops or fitted corrections. Initial raw
  one-way mean is +0.219974697 m, RMS 0.223582017 m, scatter 0.0400006395 m.
  This is above the target and uses an approximate station survey with no velocity.
  Original CSV/JSON/log are retained in `results/apollo-20170309-initial/`.
- One computed glass increment differs from its expected size by about 1.095 mm.
  Direct SPICE also loses about 0.6 mm of Earth/Moon motion for two `Time` epochs
  20 ns apart that round to the same double. An independent existing Time-valued
  interpolator recovers that motion within 14/37 micrometres. A controlled LLR
  interpolation experiment will determine the effect on the actual ranges.
- The SLR analysis stopped before any residual: the SINEX reader mistakes an
  indented ` +SNX` INPUT/HISTORY record for a block delimiter. The fix will inspect
  column 1 and extend the existing real-data functionality test; no data stripping
  or numerical workaround is used.
- A server restart removed the previous worker sessions. Commits, logs and results
  survived; no compiler/test process remained. Three replacement Sol HIGH workers
  continue in the same isolated worktrees: `slr_sol` (SINEX and SLR), `llr_sol`
  (LLR precision experiment), `compile_sol` (sole compiler allocation and verification).
