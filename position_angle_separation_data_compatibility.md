# Position-Angle and Separation Observation Data Compatibility

## Intermediate conclusion — 16 September 2026

The literature supports three definite additions to this branch's forward `P/S` model: a selectable celestial reference pole/frame, stellar aberration applied to each sight line, and optical atmospheric refraction applied to each sight line. These must be independently selectable: an archived observation called **apparent** need not still contain refraction. A fourth, conditional extension is a tangent-plane/projected-polar output convention for data that have not been reduced to spherical `P/S`.

This is an implementation conclusion from the papers reviewed below, not an exhaustive census of every NSDB record. I am confident that these cover the principal remaining direction/convention requirements in those papers. I am not claiming that implementing them makes every historical record uniquely interpretable. Missing reduction metadata remain a real limitation. Weighting and estimation strategy are outside this review.

### Compilation and validation status

- The Tudat Python kernel build completed and the new observation model was successfully exercised through `tudatpy.kernel`.
- A focused C++ test now evaluates all 60 `pm0001` observations through `PositionAngleAndSeparationObservationModel` using independent light-time solutions for Pluto and Charon.
- On 15 September 2026, `cmake --build build --target test_observation_models_PositionAngleAndSeparationObservationModel -j6` completed successfully. Running the resulting executable completed all three test cases with no errors.
- The verified Tudat C++ pre-fit RMS is **4.222797 mas in separation** and **4.409502 mas in the transverse position-angle direction**.
- The selectable reference-pole/frame adaptation is now implemented locally. Its focused C++ target compiles and all four tests pass, including the original J2000 Pluto--Charon regression and the new true-of-date Mars-satellite test. The Python kernel also builds and all new enum values and ancillary factories import successfully.
- The stellar-aberration adaptation is implemented locally and the focused real-data test passes. The shared correction is compiled in a `.cpp` file and is used by both the pre-existing pixel observable and the P/S observable. The P/S ancillary settings now select astrometric or aberrated directions independently of the reference pole, with astrometric retained as the default.

## Current validation datasets

### Simple J2000 geometry benchmark, with a photocentre caveat: IMCCE `pm0001`

- Data: <https://nsdb.imcce.fr/obspos/OBS_COLL/P/pm0001.txt>
- Metadata: <https://nsdb.imcce.fr/obspos/OBS_COLL/P/pm0001.html>
- Reference: D. J. Tholen and M. W. Buie, *The Orbit of Charon. I. New Hubble Space Telescope Observations*, Icarus 125, 245–260 (1997).
- Contents: 60 HST Pluto–Charon observations from 1992–1993.
- Coordinates: directly observed separation and position angle.
- Convention: astrometric, J2000, topocentric to HST. Although historical metadata describe a dynamical time convention, the numerical Julian dates in the distributed file match the archived HST exposure midpoints when interpreted as **UTC**, not TT/TDT.
- The data file includes the fitted-orbit residuals published with the observations.

This is a simple match to the current geometry, but calling it a proven **clean** validation set was too strong. The original paper specifies FK5 J2000 mean-equator orientation and centre-of-light image measurements; it discusses albedo-related centre offsets separately. I have not established that the archived values are fully centre-of-body corrected. The photocentre implementation on the other branch cannot be used here, so this limitation must remain explicit. See the original paper, especially pp. 248–249 and its Appendix, pp. 258–259. [Tholen & Buie (1997), author-hosted full paper](https://www.boulder.swri.edu/~buie/biblio/pub024.pdf).

All 60 numerical epochs were matched one-to-one to 60 unique MAST HST exposures from proposals 3848 and 5020. The largest difference between an NSDB epoch and the corresponding exposure midpoint was 5.658 ms. Interpreting the NSDB value as TT instead would shift the 1992 epochs by 58.184 s and is therefore inconsistent with the exposure records.

The compact SHF headers provide geocentric J2000 `POSTNST*` and `VELOCST*` HST states at their scheduled state epoch. Propagating those states linearly to the exposure midpoint and repeating the calculation showed that substituting Earth's centre changes the predicted separation by at most 0.001626 mas and the transverse position-angle component by at most 0.001043 mas. The RMS changes are 0.001016 mas and 0.000425 mas, respectively. Consequently, the C++ regression test deliberately uses Earth's centre and documents a conservative sub-0.002-mas approximation instead of embedding uncertain historical spacecraft states.

The Tudat C++ calculation uses JPL PLU060, separate converged one-way light times, and J2000 astrometric directions. Over all 60 records it gives:

- 4.222797 milliarcseconds RMS in separation.
- 4.409502 milliarcseconds RMS in the transverse position-angle direction, expressed as `rho * wrapped_delta(P)`.
- The residuals from the 1997 fitted orbit stored in the IMCCE file have RMS values of 2.93 and 2.70 milliarcseconds, respectively.
- Omitting light time degrades the separation RMS to about 82 milliarcseconds and severely degrades position angle.

The earlier 4.43/4.49-mas result was produced by incorrectly treating the NSDB Julian dates as TT. It is superseded by the UTC result above. The C++ test performs UTC-to-TDB conversion with SPICE's date-string conversion rather than applying uniform Julian-date arithmetic.

To keep the test fixture small, a 604-KiB Type-9 SPK was sampled from PLU060 every three hours over the validation interval. At every receive epoch and both iterated transmit epochs in this dataset, it agrees with full PLU060 to better than 0.000082 m for Earth, 0.002571 m for Pluto, and 0.002388 m for Charon. Its SHA-256 digest is `846362f034e547fc9a7217929f985d1db045013e15a5bdc8a5e6646473642c16`.

### Validation progress log

- Test data added inside `tests/test_tudat/data.zip`: the original `pm0001.txt`, a provenance README, and `plu060_pm0001_subset.bsp`.
- Focused build command: `cmake --build build --target test_observation_models_PositionAngleAndSeparationObservationModel -j6`.
- Test command: `build/tests/test_observation_models_PositionAngleAndSeparationObservationModel --log_level=test_suite`.
- Result: three test cases passed; the real-data case evaluated 60 records and matched the independently reproduced RMS values to 0.000001 mas.
- Residual definitions: `separation_observed - separation_computed`; wrapped position-angle difference via `atan2(sin(delta), cos(delta))`; transverse residual `separation_observed * wrapped_delta_position_angle`.
- No propagation, estimation, observation partials, relativistic angular deflection, photocentre correction, or catalogue reduction is used.

### Ground-based extension set: IMCCE `mm0012`

- Data: <https://nsdb.imcce.fr/obspos/OBS_COLL/M/mm0012.txt>
- Metadata: <https://nsdb.imcce.fr/obspos/OBS_COLL/M/mm0012.html>
- Reference: D. H. P. Jones, A. T. Sinclair, and I. P. Williams, *Secular acceleration of Phobos confirmed from positions obtained on La Palma*, MNRAS 237, 15P–19P (1989).
- Contents: 166 Phobos/Deimos observations from the Jacobus Kapteyn Telescope on La Palma in 1988.
- Coordinates: separation and position angle, with 0.2-arcsecond quoted noise.
- Convention: astrometric, topocentric, true equator and equinox of date, UTC/UT1.
- Reduction: marked as corrected for phase.

This dataset exercises the new date-dependent position-angle pole. A pure frame rotation does not change separation; position angle requires a date-dependent pole or conversion of the observation into J2000. The implementation uses the observation reception epoch unless an explicit reference epoch is supplied and evaluates the IAU 1976/1980 precession-nutation matrix in TT. Only the pole is required: changing the equinox origin rotates axes about that pole and does not alter local celestial north/east.

The receiver uses the Jacobus Kapteyn Telescope coordinates published by the Isaac Newton Group (28 deg 45 min 40.1 sec N, 17 deg 52 min 41.2 sec W, 2364 m) and Tudat's IAU-2006 GCRS/ITRS rotation with Earth-orientation data. The station position is formed on a WGS84 ellipsoid. At the roughly 0.2-arcsecond precision of these observations, uncertainty in the historical site realization is negligible.

The validation uses a compact 824-KiB Type-9 SPK made from JPL Horizons geometric ICRF vectors over 1988-09-19 through 1988-10-02. Horizons identifies DE441 for Earth and MAR099 for Mars, Phobos, and Deimos. The state vectors are sampled every five minutes and interpolated at degree 15. Against fresh Horizons vectors at 288 independent half-grid epochs per body, the largest position difference is 1.19 m. The kernel SHA-256 is `a504982dfa91452a20e0babdf2f81245f0a8968d9715384a1b7115acc3d85c07`.

All 166 accepted records are evaluated with their documented body ordering, UTC-to-TDB conversion, separate converged light times, and the JKT receiver. Current pre-fit RMS values are:

- Separation: **0.183369789 arcsec**.
- Transverse position angle with a deliberately incorrect J2000 pole: **0.150921158 arcsec**.
- Transverse position angle with the documented IAU-1976/1980 true-of-date pole: **0.152725650 arcsec**.

The true-of-date frame correction itself is about 6 milliarcseconds transverse, far below the individual 0.2-arcsecond uncertainties. It does not reduce this unweighted RMS because the series has a larger positive position-angle bias and its differential-refraction status is unknown. The test therefore compares against independently reproduced numerical results; it does **not** use "RMS becomes smaller" as a validity criterion. An independent Python calculation using ERFA `pnm80`, direct SPICE light-time iteration, and the same JPL vectors reproduced the sign and magnitude (Earth-centre approximation: 0.152827 arcsec true-of-date RMS and a -0.005981-arcsec mean transverse true-of-date-minus-J2000 correction).

The frame is unambiguous in the NSDB metadata, but this is not yet a clean atmospheric-reduction benchmark. The original short paper describes calibration using star trails and a double star, but leaves a fuller reduction analysis to later work. Its text does not establish whether differential refraction was removed. [Jones, Sinclair & Williams (1989), original paper](https://articles.adsabs.harvard.edu/pdf/1989MNRAS.237P..15J).

#### True-of-date implementation log

- Added ancillary selectors for J2000, B1950, IAU-1976 mean-of-date, IAU-1976/1980 true-of-date, IAU-2006 mean-of-date, IAU-2006/2000A true-of-date, and a custom ICRF/J2000 pole.
- Moved non-templated reference-pole evaluation into `src/tudat/astro/observation_models/positionAngleAndSeparationObservationModel.cpp` so future changes do not require recompiling every template consumer.
- Preserved J2000 as the default and separation invariance under every pole choice.
- Added equivalent Python enum and ancillary-settings factories.
- Added `mm0012.txt`, its provenance documentation, and `mar099_de441_mm0012_subset.bsp` to `tests/test_tudat/data.zip`; archive integrity passes `unzip -t`.
- Focused build command: `cmake --build build --target test_observation_models_PositionAngleAndSeparationObservationModel -j6`.
- Focused test command: `build/tests/test_observation_models_PositionAngleAndSeparationObservationModel --log_level=test_suite`.
- Python build command: `cmake --build build --target kernel -j6`; result: `kernel.so` linked successfully. A direct import check created every predefined-frame setting and a custom-pole setting.
- The first build of the new test exposed two `std::make_pair` values inferred with `const char*`; they were changed to explicit `std::string` station names. The next build completed. The first test run rejected an invalid assumption that the correct frame must lower noisy pre-fit RMS; the independent calculation confirmed the implementation and the regression now checks the reproduced physical values instead.

### Apparent-direction extension set: Qiao et al. (1999), Table 4

- Paper and data table: <https://aas.aanda.org/articles/aas/pdf/1999/10/ds7914.pdf>
- Reference: Qiao, Shen, Liu, and Harper, *1994-1996 CCD astrometric observations of Saturn's satellites and comparison with theories*, A&AS 137, 1-5 (1999), DOI 10.1051/aas:1999238.
- Contents: 41 Enceladus/Tethys/Dione/Rhea positions relative to Titan on UTC JD 2450376, observed at Sheshan.
- Convention: apparent topocentric polar coordinates. The paper explicitly says that differential refraction was removed and that stellar aberration and parallax were not removed.
- Receiver: Sheshan at 31 deg 05 min 46.1 sec N, 121 deg 11 min 03.3 sec E, altitude 97 m, on WGS84 with Tudat's IAU-2006 GCRS/ITRS rotation.

This is the unambiguous switch combination needed to validate the next adaptation: true-equator/equinox-of-date north, stellar aberration enabled, no forward optical-refraction correction, and a terrestrial topocentre. Titan is the first/reference line of sight and each listed satellite is the second line of sight, so the model position angle is that of the satellite relative to Titan.

The test uses all 41 published rows and a compact 288-KiB Type-9 SPK generated from five-minute JPL Horizons vectors over 1996-10-18 through 1996-10-21. Horizons identifies DE441 for Earth and SAT441L for Enceladus, Tethys, Dione, Rhea, and Titan. Fresh off-grid Horizons states on the observing night agree at the sub-metre level. The kernel SHA-256 is `bac178e2dc038637e286a4c0f3f548dd34cf8992fee73352c604a24ffbc6c102`.

The resulting Tudat pre-fit values are:

- Astrometric separation RMS: **0.211648598 arcsec**.
- Aberrated/apparent separation RMS: **0.214326363 arcsec**.
- Astrometric transverse position-angle RMS: **0.079196664 arcsec**.
- Aberrated/apparent transverse position-angle RMS: **0.079340229 arcsec**.
- Maximum differential stellar-aberration correction: **0.007612263 arcsec** in separation and **0.000386417 arcsec** transverse in position angle.

The roughly 0.079-arcsecond transverse residual is consistent with the paper's reported approximately 0.08-arcsecond precision. The larger separation RMS contains systematic pre-fit offsets. Aberration is only several milliarcseconds in this close-pair observable and therefore does not have to reduce the unweighted RMS. The regression checks reproduced values and the correction magnitude rather than using an unjustified “RMS improves” criterion.

An independent SpiceyPy/ERFA calculation used direct iteration of each Titan/satellite light time, DE441/SAT441L state vectors, the IAU-1976/1980 `pnm80` pole, and the exact Lorentz aberration formula. With Earth's centre as the receiver it obtained maximum corrections of 0.007498 arcsec in separation and 0.000381 arcsec transverse in position angle, within 0.12 mas of the full Sheshan result. The remaining difference is consistent with the station position and diurnal velocity omitted from that independent calculation.

#### Stellar-aberration implementation log

- Added `PositionAngleDirectionType` with explicit `astrometric` and `aberrated` ancillary choices; existing callers remain astrometric by default.
- Moved the pixel observable's existing aberration formula into the shared non-templated `stellarAberrationCorrection.cpp` implementation and reused it in P/S.
- Applied the correction separately to both independently retarded lines of sight using the common receiver inertial velocity at reception, including annual and diurnal velocity for a ground station.
- Added equivalent Python enum and ancillary-factory arguments.
- Added `qiao1999_table4.txt`, provenance documentation, and `sat441_de441_qiao1999_subset.bsp` to `tests/test_tudat/data.zip`.
- No propagation, estimation, observation partials, relativistic angular deflection, photocentre correction, catalogue reduction, weighting, or atmospheric refraction is used.

### Refraction-on dataset audit: deferred for lack of a clean P/S regression

After completing the aberration milestone, I rechecked the downloaded NSDB P/S metadata specifically for records marked as retaining differential refraction. NSDB distinguishes `Done`, `None`, and `Unknown` for differential-refraction correction in its current schema. The result is restrictive:

- Every P/S collection explicitly labelled “not corrected” or “not excluded” is a visual/micrometer series from **1866–1922**. Examples include [`nm0028`](https://nsdb.imcce.fr/obspos/OBS_COLL/N/nm0028.html) (Triton, Halsted, 1883–1888), [`nm0034`](https://nsdb.imcce.fr/obspos/OBS_COLL/N/nm0034.html) (Triton, Lick, 1898–1902), and the Uranian series [`um0055`](https://nsdb.imcce.fr/obspos/OBS_COLL/U/um0055.html) through `um0066` from Washington/Yerkes. Their metadata provides the site, telescope, apparent topocentric/date-frame convention, and observation time, but not observing wavelength or contemporaneous pressure, temperature, and humidity.
- `mm0001` says “not corrected or no information”, which is explicitly ambiguous and therefore unusable for this validation gate.
- The modern [Shen et al. (2002) Sheshan Uranus work](https://doi.org/10.1051/0004-6361:20020872) is physically informative but not a P/S validation dataset. It specifies a 900 nm filter and applies rigorous Woolard--Clemence refraction corrections (up to about 0.03 arcsec) during reduction, then deliberately publishes **raw pixel coordinates** so readers can redo the calibration. Testing those data would require the camera/plate calibration path that this task excludes.
- Qiao et al. (1999), used for the aberration test above, explicitly removed differential refraction from its published P/S. It is a clean refraction-**off** case, not evidence for a forward refraction implementation.

Consequently, atmospheric refraction is not implemented in this milestone. A standard-atmosphere assumption could produce a plausible correction for the historical data, but it would not be an unambiguous reproduction of the physical conditions under which the observations were made. Adding such a model and fixing its output as a regression value would test our assumption rather than the archived data. This deferral follows the agreed rule that each adaptation must first have a real P/S dataset whose convention and required model are explicitly documented.

A suitable future case needs, at minimum, retained-refraction P/S, a known terrestrial site and time, wavelength or passband, and either measured meteorology or an explicit published standard-atmosphere prescription. If such a reference is supplied, the next implementation should bend both independently retarded apparent directions toward the local zenith before forming P/S; it must remain independent of the aberration and celestial-pole choices.

### Large but heterogeneous set: IMCCE `sm0034`

- Data: <https://nsdb.imcce.fr/obspos/OBS_COLL/S/sm0034.txt>
- Metadata: <https://nsdb.imcce.fr/obspos/OBS_COLL/S/sm0034.html>
- Contents: 17,297 explicitly encoded astrometric position-angle/separation records within a much larger Saturn-satellite catalogue.
- Limitations: the position-angle records are predominantly from 1874–1973, use equator/equinox of date, span many observatories, and include missing or unknown station information.

This is useful as a later robustness set, not as the first implementation validation.

## Remaining model additions

The list below contains only physical or reference-convention capabilities that are still missing from this branch's `P/S` observable. Independent light times, arbitrary receiver link ends, high-accuracy terrestrial stations and Earth rotation, time-scale conversion, subject/reference ordering, units and component ordering, optional components, and residual wrapping are already available or are ingestion bookkeeping rather than new observation models.

Relativistic angular light deflection, phase/photocentre correction, and star-catalogue astrometric reduction are implemented on other branches. They are therefore not additions to make here and cannot be enabled in the present baseline validation test. Their absence must still be respected when selecting data and interpreting residuals.

Line counts below are rough engineering estimates, not measured patches. They cover the named model option and integration/tests, not a general historical re-reduction package or a complete extension of all analytical partials.

- **Differential optical atmospheric refraction**
  - **Confidence:** High about the physics; medium or low for applying it to an individual archived series because wavelength and meteorology are often absent.
  - **Current Tudat support:** Missing for optical angular observations. Existing radio/tropospheric light-time corrections are not the required direction correction.
  - **Estimated size:** 200–400 production lines plus 150–250 binding and test lines.
  - **Physical model:** Apply wavelength- and atmosphere-dependent bending as a function of zenith distance to both sight lines. Absolute refraction largely cancels for a close pair, but its differential angular gradient does not.

### Conditional extension: projected-polar observation geometry

- **When needed:** Only when published `P/S` still represents calibrated focal-plane distances/orientations rather than great-circle separation and celestial position angle. The fact that an observation originated on a CCD is not sufficient evidence that this option is needed.
- **Confidence:** High that the geometries differ; medium that this is necessary for any particular archived polar series until its reduction is checked.
- **Current Tudat support:** The existing `P/S` helper already computes exact spherical geometry. Tudat also has pixel-coordinate projection machinery, but there is no selectable projected-polar representation in this `P/S` model.
- **Estimated size:** Approximately 100–250 production lines plus 100–200 integration/test lines if existing projection utilities are reusable. Instrument-specific distortion/re-reduction is excluded.
- **Geometric model:** Project both directions onto a specified tangent plane, then form their plane-coordinate difference and its angle/length using the documented calibration convention. This needs the optical centre and orientation. It is not an extra spherical-curvature correction to the current exact `P/S` formula. Undocumented plate geometry cannot be reconstructed uniquely.

## Literature: what the forward models actually use

- **Desmars, Vienne & Arlot (2009), COSS08, A&A 493, 1183–1195, Sections 3.3–3.7.** Their Saturn catalogue treats refraction and aberration separately and documents whether each is removed, retained, or only presumed corrected. For unremoved optical refraction they use a Laplace model with atmosphere- and wavelength-dependent coefficients, noting loss of accuracy at large zenith distance. They also warn that some published differential coordinates are actually tangent-plane coordinates and that the optical centre is often unavailable. Their historical default assumptions about missing correction information are not definitive evidence of how an individual series was reduced. **Implication:** independent correction options and explicit unknown metadata are necessary; a generic `apparent` label is insufficient. [Full paper, HAL](https://hal.science/hal-00588233/document). The catalogue schema also records B1950, true-of-date, J2000, and historical year-start mean frames. [CDS COSS08 format](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/A%2BA/493/1183?format=html&tex=true).

- **Yuan et al. (2021), Neptunian-satellite catalogue, A&A 645, A48, Sections 3.3.2, 3.5.1 and 3.6.** They explicitly convert observed position angles between celestial reference systems using the angle between the two local north directions; separation is unchanged under this frame rotation. Their catalogue independently records whether aberration and refraction are included. True-of-date conversion uses precession, nutation and frame bias. Their discussion of FK4 absolute coordinates is more involved than a pole rotation and should not be indiscriminately transferred to every historical position angle. **Implication:** a selectable north/frame has direct precedent in processing archived `P/S`; it is not merely a hypothetical feature. [Full original paper, author-uploaded copy](https://www.researchgate.net/publication/348331406_New_precise_positions_in_2013-2019_and_a_catalog_of_ground-based_astrometric_observations_of_11_Neptunian_satellites_1847-2019_based_on_Gaia-DR2).

- **Qiao, Shen, Liu & Harper (1999), 1994–1996 Saturn CCD observations, A&AS 137, 1–5, Sections 3.3–4.** They form polar coordinates from measured image-coordinate differences using a calibrated scale and orientation. Crucially, their published `P/S` are **apparent topocentric coordinates with differential refraction removed, but stellar aberration retained**. Their computed positions include aberration. The refraction contribution is reported as at most about 20 mas for their observing geometry. **Implication:** there are real polar datasets requiring aberration on and refraction off. Their image-plane definition also motivates checking projected versus spherical polar geometry before assuming an exact match. [Full original paper, author-uploaded copy](https://www.researchgate.net/publication/45678823_1994-1996_CCD_astrometric_observations_of_Saturn%27s_satellites_and_comparison_with_theories).

- **Shen et al. (2002), 1995–1997 Uranian-satellite astrometry, A&A 391, 775–779, Section 3 and appendix.** Their calibration comparison includes aberration and differential refraction in both separation and position angle, with refraction effects reaching approximately 30 mas. They neglect projection nonlinearity after estimating it below approximately 5 mas for their small field. They also publish uncorrected raw pixel positions, which must not be confused with the corrected quantities used during calibration. **Implication:** optical bending affects both polar components, and the published representation—not just the corrections mentioned in a paper—determines the required forward model. [Full original paper, author-uploaded copy](https://www.researchgate.net/publication/252987841_Astrometry_of_five_major_Uranian_satellites_in_1995-1997).

- **Morgado et al. (2016), Galilean mutual approximations, MNRAS 460, 4086–4097, Sections 3–4.** Their detector-distance prediction uses corrected directions and gnomonic projection, including optical refraction and aberration. This is related separation-based astrometry, **not a direct fit of the archived `P/S` observable**, so it supports the direction/projection pipeline but does not justify adding a new timing observable to this task. [Full paper, arXiv](https://arxiv.org/abs/1605.06573).

- **Tholen & Buie (1997), Pluto–Charon HST observations, Appendix.** They predict a retarded relative orbit projected into a north/east sky plane with J2000 orientation. The historical calculation uses a common barycentric light-time epoch and small-angle projection; our independent two-ray calculation uses more exact geometry. These are not grounds for removing Tudat's existing independent light times. More importantly, the image extraction measures centres of light, which is why `pm0001` is not yet certified as a clean centre-of-body test. [Original paper](https://www.boulder.swri.edu/~buie/biblio/pub024.pdf).

### Synthesis for implementation

The recurring practice is to predict the same angular quantity that the observer published, applying only effects still present in those data. There is no universal combination of date-frame, aberration and refraction that represents all natural-satellite `P/S`.

My proposed structure is:

1. Retain the existing two independent retarded sight lines and receiver state.
2. Apply the selected stellar-aberration transform separately to both sight lines.
3. If refraction remains in the observations, bend both directions using the local vertical and atmosphere, with explicit wavelength assumptions.
4. Define north/east in the requested celestial frame and compute spherical `P/S`; alternatively, use explicitly documented projected-polar geometry.

These are direction transformations and an output-frame/geometry choice, not new receiver, station, time-scale or light-time models. Analytical partials must eventually follow the same transformations. Before implementation, test frame invariance of separation, the sign of the position-angle frame change, aberration against an independent standard, and differential refraction in both polar components.

## Recommended implementation order

1. **Completed:** run `pm0001` through Tudat itself as an astrometric J2000 geometry benchmark, retaining the photocentre caveat and quantifying the HST-receiver approximation.
2. **Completed:** add explicit ancillary settings for the celestial reference frame, retaining astrometric J2000 as the default, and validate true-equator/equinox-of-date against all 166 `mm0012` observations from the high-accuracy JKT receiver.
3. **Completed in `66049da80`:** factor stellar aberration into a shared angular-direction correction, add an explicit astrometric/aberrated-direction option, and validate it against all 41 documented apparent topocentric measurements in Qiao et al. (1999), Table 4.
4. **Deferred after dataset audit:** optical refraction remains the lowest priority. The available refraction-on P/S records lack the atmospheric/passband metadata needed for an unambiguous real-data regression; do not implement it until a suitable reference supplies those assumptions.

Integration of the separately implemented angular light-deflection, photocentre, and star-catalogue-reduction capabilities is explicitly outside this branch and baseline test. Assess projected-polar support against a specifically documented dataset before implementing it.

## Important limitation

The aim can be full support for all sufficiently documented observations, but not automatic faithful processing of every record in the archive. Records with unknown stations, ambiguous time scales, arbitrary plate origins, or undocumented apparent-place reductions do not contain enough information for a unique modern observation model.

## Local handoff and continuation

- **Worktree:** `/home/dominic/Tudat/tudat-monorepo/tudatpy`
- **Branch:** `feature/position-angle-and-separation-observation-model-and-partials`
- **Baseline before this validation:** `f76b40736` — `Checkpoint before natural-satellite observation validation`.
- **Earlier work:** `40af0d1c3` consolidated the models; `3b0daf3bb` merged PR #857 residual wrapping into PR #860's branch.
- **Local milestone content:** the selectable-frame and stellar-aberration implementations and bindings, three real-data C++ validation cases, the seven archived files in `tests/test_tudat/data.zip`, and this report. Preserve the unrelated dirty `examples/tudatpy` submodule and other pre-existing untracked files. Nothing has been pushed.
- **Relevant code:** `include/tudat/astro/observation_models/positionAngleAndSeparationObservationModel.h`; `src/tudat/astro/observation_models/positionAngleAndSeparationObservationModel.cpp`; `include/tudat/astro/observation_models/stellarAberrationCorrection.h`; `src/tudat/astro/observation_models/stellarAberrationCorrection.cpp`; `include/tudat/astro/observation_models/pixelCoordinatesObservationModel.h`; and `src/tudatpy/estimation/observations_setup/ancillary_settings/expose_ancillary_settings.cpp`.
- **Build:** No build remains active. The focused P/S target and Python kernel built successfully with `-j6`; all five P/S test cases pass. The existing pixel-coordinate and PSF/pixel stellar-aberration tests also pass after the shared-code refactor. The Python import check verified the astrometric default and both aberrated factory paths. The formatter hook, `git diff --check`, and complete `unzip -t tests/test_tudat/data.zip` check pass. Python environment: `/home/dominic/miniconda3/envs/tudatpy-dev`.
- **Temporary research inputs:** `/tmp/ps-literature/` contains downloaded papers and extracted text; `/tmp/tudatpy-hst-validation/` contains retrieved HST header products; `/tmp/tudatpy-pm0001-subset-20260915-a/` contains the Pluto compact-SPK inputs; `/tmp/tudatpy-saturn-aberration/` contains the Qiao/SAT441L generation and independent-validation scripts. Temporary paths may disappear after reboot; permanent sources and the required regression fixtures are recorded in the repository data archive.
- **Outstanding work:** none for the clean model adaptations that can currently be validated. Optical atmospheric refraction is explicitly deferred until an unambiguous retained-refraction P/S dataset or publication supplies the passband and atmosphere assumptions. No stellar-aberration implementation or validation work remains for this milestone.

The command below resumes this specific saved conversation and supplies this worktree as its working directory. Its syntax was checked using the OpenAI Docs skill against the [official CLI command documentation](https://learn.chatgpt.com/docs/developer-commands?surface=cli) and the installed `codex resume --help`. It does not build, commit or push anything by itself. Run it later, after leaving the current session:

```bash
codex resume -C /home/dominic/Tudat/tudat-monorepo/tudatpy 01a0a586-412e-78c1-91a9-c499b297bc91
```
