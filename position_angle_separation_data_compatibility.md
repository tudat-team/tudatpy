# Position-angle and separation: information to retain

This note identifies the documentation, validation references, and future work
to retain from the position-angle and angular-separation development work.

## Documentation with literature references

The [P/S model docstring](src/tudatpy/estimation/observable_models_setup/model_settings/expose_model_settings.cpp)
already defines the spherical observables, separate light-time solutions, reference
pole, and astrometric or stellar-aberrated directions. Add a short explanation
of why frame, aberration, and refraction conventions must be considered
separately: an archived observation labelled *apparent* does not uniquely state
which reductions remain. [Desmars et al. (2009)](https://doi.org/10.1051/0004-6361:200810203)
discuss the separate reduction metadata; [Qiao et al. (1999)](https://doi.org/10.1051/aas:1999238)
publish apparent P/S with differential refraction removed but stellar
aberration retained.

Cite [Yuan et al. (2021)](https://doi.org/10.1051/0004-6361/202038776)
for converting position angle between celestial reference systems by changing
local north while leaving angular separation unchanged. State that Tudat's
`aberrated` option applies stellar aberration, not every possible apparent-place
reduction. The current P/S formula gives spherical separation and position
angle; calibrated image-plane polar coordinates may require a different
projection. The B1950 pole choice should not be presented as a complete FK4
catalogue-coordinate transformation; see the [NAIF Frames Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/frames.html).

The explicit API lists in
[`model_settings.rst`](docs/tudatpy/source/estimation/observable_models_setup/model_settings.rst)
and [`ancillary_settings.rst`](docs/tudatpy/source/estimation/observations_setup/ancillary_settings.rst)
still need the P/S functions and ancillary choices, respectively.

## Publication references in validation tests

Keep each publication and the specific observation convention beside the
corresponding decoding or assertion in
[`unitTestPositionAngleAndSeparationObservationModel.cpp`](tests/test_tudat/src/astro/observation_models/unitTestPositionAngleAndSeparationObservationModel.cpp).
The fixture README inside `tests/test_tudat/data.zip` records source files,
compact-kernel provenance, and receiver approximations. Published post-fit
residuals supply broad comparison scales; they are not exact targets for
Tudat's unfit, unfiltered pre-fit residuals.

| Dataset | References and facts to retain |
| --- | --- |
| `pm0001` | [Tholen & Buie (1997)](https://www.boulder.swri.edu/~buie/biblio/pub024.pdf) and the [NSDB record](https://nsdb.imcce.fr/obspos/OBS_COLL/P/pm0001.html): J2000 P/S, the centre-of-light caveat, the UTC interpretation checked against MAST exposure midpoints, and the quantified Earth-centre approximation. |
| `mm0012` | [Jones et al. (1989)](https://doi.org/10.1093/mnras/237.1.15P) and [NSDB](https://nsdb.imcce.fr/obspos/OBS_COLL/M/mm0012.html): true-of-date convention and published accuracy; differential-refraction reduction remains undocumented. |
| `um0027` | Tomita & Soma (1979), [NSDB](https://nsdb.imcce.fr/obspos/OBS_COLL/U/um0027.html), and [Jacobson (2014), Table 3](https://doi.org/10.1088/0004-6256/148/5/76): native B1950 P/S and ET; Jacobson's selected post-fit sample supplies a scale, not the expected pre-fit residuals for all rows. |
| `nm1007` | Veillet (1982), Veillet & Bois (1988), [NSDB](https://nsdb.imcce.fr/obspos/OBS_COLL/N/nm1007.html), and [Yuan et al. (2021)](https://doi.org/10.1051/0004-6361/202038776): X/Y converted to P/S as a secondary B1950 frame check, with unknown atmospheric-reduction status. |
| Qiao Table 4 | [Qiao et al. (1999)](https://doi.org/10.1051/aas:1999238): identify the table, Titan as the reference, retained stellar aberration, removed differential refraction, and the published residual scale. |

The test comments already include many author-year references. Add stable
publication links or DOIs and keep the distinction between published post-fit
results and Tudat pre-fit checks explicit.

## Open points for a future issue

1. **Aberrated P/S estimation partials.** The current
   [`positionAngleAndSeparationPartial.h`](include/tudat/astro/orbit_determination/observation_partials/positionAngleAndSeparationPartial.h)
   forms Jacobians from geometric sight lines without differentiating the
   stellar-aberration transformation. Add the required derivatives and numerical
   partial checks before relying on aberrated P/S for estimation.
2. **RA/Dec conventions.** Extend angular and relative angular position
   observables to selectable date frames and astrometric or aberrated directions.
   A custom RA/Dec frame needs both a pole and an RA origin.
3. **Projected P/S, conditionally.** Support calibrated tangent-plane polar
   measurements when a dataset documents its optical centre, orientation, and
   published coordinate definition. Keep spherical P/S as the default.
4. **Composition with other reductions.** Assess existing light-deflection,
   photocentre/phase, and star-catalogue correction work together with P/S,
   including which effects remain in each dataset.
5. **Atmospheric refraction, conditionally.** This is outside the present PR.
   Revisit it only with retained-refraction P/S, a known site and time, a
   passband, and measured meteorology or a published atmosphere prescription.
6. **Heterogeneous archival records.** Later use collections such as `sm0034`
   to check handling of missing stations, ambiguous timescales, and unknown
   reductions. Such records cannot always be modelled uniquely.

Build history, temporary paths, commit IDs, handoff instructions, and exact
current pre-fit digits do not need to be retained. Validation assertions should
continue to use independently justified ranges.
