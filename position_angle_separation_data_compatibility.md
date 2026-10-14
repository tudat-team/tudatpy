# Position-angle and separation: information to retain

This note identifies the documentation, validation references, and future work
to retain from the position-angle and angular-separation development work.

## Documentation with literature references

Done: the [P/S model docstring](src/tudatpy/estimation/observable_models_setup/model_settings/expose_model_settings.cpp)
now lists the literature references, describes the (classical) stellar-aberration model and its limits, and explains that the
reference pole and direction convention are per-observation ancillary settings rather than model settings. The P/S functions and
ancillary choices are included in `model_settings.rst` and `ancillary_settings.rst`.

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

Done: the test comments give these references (with DOIs or links) and relate each tolerance to the published precision.

## Open points for a future issue

1. **RA/Dec conventions.** Extend angular and relative angular position
   observables to selectable date frames and astrometric or aberrated directions.
   A custom RA/Dec frame needs both a pole and an RA origin.
2. **Projected P/S, conditionally.** Support calibrated tangent-plane polar
   measurements when a dataset documents its optical centre, orientation, and
   published coordinate definition. Keep spherical P/S as the default.
3. **Composition with other reductions.** Assess existing light-deflection,
   photocentre/phase, and star-catalogue correction work together with P/S,
   including which effects remain in each dataset.
4. **Atmospheric refraction, conditionally.** This is outside the present PR.
   Revisit it only with retained-refraction P/S, a known site and time, a
   passband, and measured meteorology or a published atmosphere prescription.
5. **Heterogeneous archival records.** Later use collections such as `sm0034`
   to check handling of missing stations, ambiguous timescales, and unknown
   reductions. Such records cannot always be modelled uniquely.

Build history, temporary paths, commit IDs, handoff instructions, and exact
current pre-fit digits do not need to be retained. Validation assertions should
continue to use independently justified ranges.
