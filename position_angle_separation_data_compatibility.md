# Position-Angle and Separation Observation Data Compatibility

## Intermediate conclusion — 15 September 2026

The literature supports three definite additions to this branch's forward `P/S` model: a selectable celestial reference pole/frame, stellar aberration applied to each sight line, and optical atmospheric refraction applied to each sight line. These must be independently selectable: an archived observation called **apparent** need not still contain refraction. A fourth, conditional extension is a tangent-plane/projected-polar output convention for data that have not been reduced to spherical `P/S`.

This is an implementation conclusion from the papers reviewed below, not an exhaustive census of every NSDB record. I am confident that these cover the principal remaining direction/convention requirements in those papers. I am not claiming that implementing them makes every historical record uniquely interpretable. Missing reduction metadata remain a real limitation. Weighting and estimation strategy are outside this review.

### Compilation and validation status

- The Tudat Python kernel build completed and the new observation model was successfully exercised through `tudatpy.kernel`.
- A focused C++ test now evaluates all 60 `pm0001` observations through `PositionAngleAndSeparationObservationModel` using independent light-time solutions for Pluto and Charon.
- On 15 September 2026, `cmake --build build --target test_observation_models_PositionAngleAndSeparationObservationModel -j6` completed successfully. Running the resulting executable completed all three test cases with no errors.
- The verified Tudat C++ pre-fit RMS is **4.222797 mas in separation** and **4.409502 mas in the transverse position-angle direction**.
- This stage validates the existing astrometric J2000 model. It does not yet add the selectable reference-pole/frame adaptation.

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

The frame convention is not directly compatible with the current position-angle implementation because Tudat currently uses a fixed J2000 north pole. A pure frame rotation does not change separation; position angle requires a date-dependent pole or conversion of the observation into J2000. This statement concerns frame compatibility only, not certification of all other corrections.

The original short paper describes calibration using star trails and a double star, but leaves a fuller reduction analysis to later work. Its text alone does not establish all aberration/refraction conventions required for a clean test. [Jones, Sinclair & Williams (1989), original paper](https://articles.adsabs.harvard.edu/pdf/1989MNRAS.237P..15J).

MPC observatory code 950 supplies a usable geocentric site position for La Palma, but a surveyed Jacobus Kapteyn Telescope coordinate and its reference epoch would be preferable for a final high-accuracy result. Tudat's high-accuracy Earth rotation should be used with Earth-orientation data appropriate to 1988.

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

- **Selectable position-angle reference pole and celestial frame**
  - **Confidence:** High for explicit pole rotations and modern mean/true-of-date frames; medium for faithful reproduction of legacy FK4 reductions and undocumented historical formulae.
  - **Current Tudat support:** Partial. Tudat has the underlying precession/nutation and frame-rotation machinery, but the `P/S` observable currently uses a fixed J2000 pole.
  - **Estimated size:** 150–250 production lines plus 150–250 binding, documentation, and test lines.
  - **Physical model:** Define local north/east using the observation's celestial reference frame. Support ICRS, explicitly specified FK5 J2000, B1950 north, mean equator/equinox at a supplied epoch, and true equator/equinox of date, including historical year-start epochs. A user-provided pole/frame function is useful. ICRS and FK5 J2000 should not be silently conflated for precision work. A rigid pole change affects `P`, not `S`. **B1950 north does not automatically imply a full FK4 catalogue transformation:** whether legacy E-terms or other non-rigid conventions are relevant must be established from the reduction. The separate star-catalogue implementation is not to be duplicated here.

- **Stellar aberration / astrometric versus aberrated direction**
  - **Confidence:** High for the required architecture and aberration physics; medium for classifying individual historical datasets whose use of `apparent` is underspecified.
  - **Current Tudat support:** Partial. The retarded astrometric geometry is present. Stellar-aberration code exists for the pixel-coordinate observable but is not a shared correction available to `P/S`.
  - **Estimated size:** 150–300 production lines plus 150–250 binding and test lines if the existing implementation is factored into a shared utility.
  - **Physical model:** Transform each retarded sight line using receiver velocity, including annual and diurnal contributions for a ground receiver. This must act on both directions before calculating `P/S`: differential aberration does not vanish merely because the pair is close. Make it an independent ancillary choice, not a switch that also forces refraction or true-of-date north. The existing pixel implementation uses the inverse of Jacobson's stellar-aberration correction; verify its accuracy against the selected standard before treating it as a general high-accuracy implementation. This is angular aberration, not another light-time solver or angular gravitational-deflection model.

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
2. **Next:** add explicit ancillary settings for direction convention and celestial reference frame, retaining astrometric J2000 as the default.
3. Add true-equator/equinox-of-date support. Use `mm0012` only after resolving its remaining reduction conventions, with a high-accuracy La Palma receiver state.
4. Factor stellar aberration into a shared angular-direction correction and add an explicit apparent-direction option.
5. Add optical refraction only where the observation metadata provide enough information.

Integration of the separately implemented angular light-deflection, photocentre, and star-catalogue-reduction capabilities is explicitly outside this branch and baseline test. Assess projected-polar support against a specifically documented dataset before implementing it.

## Important limitation

The aim can be full support for all sufficiently documented observations, but not automatic faithful processing of every record in the archive. Records with unknown stations, ambiguous time scales, arbitrary plate origins, or undocumented apparent-place reductions do not contain enough information for a unique modern observation model.

## Local handoff and continuation

- **Worktree:** `/home/dominic/Tudat/tudat-monorepo/tudatpy`
- **Branch:** `feature/position-angle-and-separation-observation-model-and-partials`
- **Baseline before this validation:** `f76b40736` — `Checkpoint before natural-satellite observation validation`.
- **Earlier work:** `40af0d1c3` consolidated the models; `3b0daf3bb` merged PR #857 residual wrapping into PR #860's branch.
- **Local changes for this validation:** the focused C++ test, the three archived fixture files in `tests/test_tudat/data.zip`, and this report. Preserve the unrelated dirty `examples/tudatpy` submodule and other pre-existing untracked files. Nothing has been pushed.
- **Relevant code:** `include/tudat/astro/observation_models/positionAngleAndSeparationObservationModel.h`; `include/tudat/simulation/estimation_setup/createObservationModelFactory.h` (`getJ2000NorthPoleDirectionInGlobalFrame`); `include/tudat/astro/observation_models/pixelCoordinatesObservationModel.h` (existing aberration utility).
- **Build:** No build remains active. The focused C++ target built successfully with `-j6`, and its executable passed. Python environment: `/home/dominic/miniconda3/envs/tudatpy-dev`.
- **Temporary research inputs:** `/tmp/ps-literature/` contains downloaded papers and extracted text; `/tmp/tudatpy-hst-validation/` contains retrieved HST header products; `/tmp/tudatpy-pm0001-subset-20260915-a/` contains the generated compact SPK. Temporary paths may disappear after reboot; permanent sources and the required regression fixtures are recorded in the repository data archive.
- **Outstanding work:** implement and validate the selectable position-angle reference pole/celestial frame as the first model adaptation. Stellar aberration follows only after that case is independently validated; optical atmospheric refraction remains lowest priority.

The command below resumes this specific saved conversation and supplies this worktree as its working directory. Its syntax was checked using the OpenAI Docs skill against the [official CLI command documentation](https://learn.chatgpt.com/docs/developer-commands?surface=cli) and the installed `codex resume --help`. It does not build, commit or push anything by itself. Run it later, after leaving the current session:

```bash
codex resume -C /home/dominic/Tudat/tudat-monorepo/tudatpy 01a0a586-412e-78c1-91a9-c499b297bc91
```
