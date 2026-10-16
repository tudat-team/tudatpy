# Fixed-orbit LAGEOS residual reference cases

Research checkpoint: 2026-10-15. No Tudat, Orekit or GROOPS execution, builds,
propagation, estimation or residual-based data selection has been performed.
Downloaded provider bytes are unchanged. `../data/slr/download-manifest.json`
records URL, SHA256, byte count, retrieval time and rights notes for each file.

## Recommended first case

Use **LAGEOS-1, 2024-03-31 through 2024-04-06**, with the final weekly ILRS-A
combined orbit, plus the ASI individual orbit as a separately reported comparison.
Start with **Graz 7839 and Herstmonceux 7840**, selected from documented station
quality and available observations, before looking at any residuals. The seven
daily CRDv2 files contain 19 normal points at Graz and 10 at Herstmonceux. These
are counts before elevation/validity selection; the last Graz session crosses
midnight and the SP3 endpoint must bound usable measurements. A larger fixed
selection can add Yarragadee 7090 (180 normal points in the same daily files).

Public paired data are now local under
`/home/dominic/Tudat/laser-ranging-work/data/slr/2024-04/`:

| Product | Local filename | Provider facts inspected |
|---|---|---|
| Combined postfit orbit | `ilrsa.orb.lageos1.240406.v81.sp3.gz` | `FIT COMB`, `SLR20`, comment explicitly SLRF2020; UTC; 120 s samples; 5040 epochs |
| Single-centre postfit orbit | `asi.orb.lageos1.240406.v80.sp3.gz` | `FIT ASI`; comment explicitly SLRF2020; UTC; 120 s; 5041 epochs |
| Orbit combination report | `ilrsa.orb.lageos1.240406.v81.sum.gz` | AC-versus-combined radial/cross/along differences, editing and fitted frame parameters; these are not observation residuals |
| Actual observations | `lageos1_20240331.np2`, `lageos1_20240401.np2` through `lageos1_20240406.np2` | CRDv2 normal points, system/configuration records, meteorology and calibration information |
| Matching weekly ASI solution | `asi.pos+eop.240406.v80.snx.gz` | Estimated station positions at arc mean epoch, daily EOP, actual analysis description |
| Weekly combined solution | `ilrsa.pos+eop.240406.v80.snx.gz` | Combined station/EOP solution; retain its distinct version from orbit v81 |

Current archive links (the older published directory layout no longer works):

- [Final orbit directory](https://edc.dgfi.tum.de/pub/slr/products/orbits/TRF/lageos1/2024/240406/).
- [Normal-point directory](https://edc.dgfi.tum.de/pub/slr/data/npt_crd_v2/lageos1/2024/).
- [Weekly station/EOP directory](https://edc.dgfi.tum.de/pub/slr/products/pos+eop/weekly/2024/240406/).
- [ILRS orbit product description](https://edc.dgfi.tum.de/en/products/precise-orbits/)
  and orbit dataset DOI **10.5067/SLR/SLR_ILRSORB_001**.

The weekly orbit is fitted to the SLR observations. It is suitable for testing
the observation reduction against a fixed postfit ephemeris, but is not an
independent validation of orbit accuracy using withheld measurements. It is not
a CPF prediction. The SP3 clock entries are missing-value sentinels, not a
satellite clock product needed for this passive two-way measurement.

## Precise conventions and fixed corrections

1. **Frame and scale.** Treat this SP3 as Earth-centred, Earth-fixed SLRF2020,
   with TT-compatible terrestrial coordinates, positions in km and SP3 velocity
   units of dm/s. Read the actual `%c ... UTC` field; do not assume GPS time from
   GNSS SP3 habits. SLRF2020 station coordinates, valid solution intervals,
   velocities, any required PSD, and telescope-reference eccentricities must be
   applied consistently. A generic Earth-fixed label does not supply a TRF
   transformation or a station solution. The local
   `metadata/SLRF2020_POS+VEL_2023.10.02.snx` is pinned, as well as both weekly
   solutions. Distinguish an external station-history experiment from a
   reconstruction using already-estimated weekly station coordinates; publish
   both configurations rather than fitting station offsets here.
2. **Analysis model from the actual ASI SINEX.** Its `FILE/COMMENT` declares
   ITRF2020 a priori positions/velocities, monument coordinates, IERS2010 station
   uplift/pole tide, GOT4.10c ocean loading, no atmospheric-pressure loading,
   Mendes–Pavlis troposphere, 20-degree elevation cutoff, and relativistic light
   time. Positions and daily EOP were estimated by the provider. Root-level
   `README_AC.asi` is obsolete (SLRF2008/GGM02C/GOT4.7), so do not use it in
   preference to the dated solution. Adding all modern loading models is a
   distinct scientific configuration, not exact reproduction of this fit.
3. **Observation epoch and value.** CRD record11 holds the time of flight;
   report `O-C = c*(observed TOF - computed TOF)/2` in one-way metres. Respect
   each epoch-event flag (ground transmit, bounce, receive) and midnight
   rollover. The selected Graz/Herstmonceux records use event2 (transmit).
   H2 time-scale4 is **UTC(GPS)**, not GPS coordinate time; CRD also defines3 as
   UTC(USNO),7 as UTC(BIPM). UTC time tags and the measured interval are different
   quantities: do not subtract large scalar TDB epochs to recover millimetre
   light time. Preserve split/high-precision events through both legs.
4. **Applied flags.** Selected Graz/Herstmonceux H4 records have refraction0,
   CoM0, amplitude0, station-system-delay1, range-type2 and quality0. Thus CoM
   and atmospheric delay still need treatment; calibration delay must not be
   applied twice. Matera2024 in the same files sets amplitude1, illustrating why
   flags cannot be inferred for an entire satellite/day. Optional record12 is
   supplementary information, not proof that its correction is already applied.
   Preserve per-session and per-system identifiers.
5. **Centre of mass / target signature.** The pinned
   `metadata/ILRS_Data_Handling_File_2024.02.13.snx` explicitly requires the
   **2022-09-15** target-signature model. Therefore use pinned
   `metadata/com6.220915.tar.gz`, not today's model or a universal251mm. Its
   LAGEOS-1 entries give Graz532nm **250.4mm**, Herstmonceux532nm **244.5mm**
   for these dates. LAGEOS-2 values differ. The newer `com.231128.tar.gz` is
   retained for provenance comparison, not silently substituted. The table
   `systems.dat` records detector/processing configurations; reported decimal
   precision is not a0.1mm physical accuracy claim. A positive CoM length adds
   to an uncorrected reflector range to refer it to the centre; equivalently
   subtract that length from a centre-of-mass theoretical range. Implement one
   convention and record it once.
6. **Published station corrections.** The same DHF contains mandatory
   wavelength/satellite/time-dependent range, time, meteorological and exclusion
   entries. Its legend says supplied range values are **subtracted from one-way
   observations**. Current open-ended LAGEOS-1 range entries include Graz
   **+3.9mm** and Herstmonceux **-2.2mm**. These are externally published input
   corrections, not a mean fitted to our residuals. Do not mix this DHF with a
   different CoM release, or apply its model bias and an independently fitted
   station bias twice. Handle E/Q/N/exclusion/time records explicitly; report
   unsupported estimation recommendations rather than performing estimation.
7. **Light propagation.** Solve transmit→bounce→receive with moving station and
   satellite, station deformation at the relevant events, optical atmospheric
   delay from recorded wavelength/meteorology on both legs, and Earth Shapiro
   delay. IERS2010§11.2 uses the geocentric formulation for near-Earth ranging.
   If using Tudat's barycentric light-time path instead, station and satellite
   geocentric vectors both require the consistent TT-compatible→TDB-compatible
   transformation of Eq11.19 and consistent measured-interval conversion; a
   raw rotated ECEF orbit plus Earth translation alone is not a complete
   millimetre BCRS model. No propagation force model needs rebuilding when the
   orbit is fixed.

Primary references:

- [CRD v2.01 manual](https://ilrs.gsfc.nasa.gov/docs/2022/crd_v2.01e3.pdf),
  H2, H4, records11/12/20/40; [v1.01 manual](https://ilrs.gsfc.nasa.gov/docs/2009/crd_v1.01.pdf)
  for the2018 secondary case.
- [ILRS CoM guidance](https://ilrs.gsfc.nasa.gov/missions/satellite_missions/current_missions/lag1_com.html)
  explicitly links DHF bias values to the selected CoM model;
  [Rodríguez, Appleby & Otsubo2019](https://doi.org/10.1007/s00190-019-01315-0),
  with [official model archive](https://datos-geodesia.ign.es/SLR/centre_of_mass_models/models/).
- [DHF2024-02-13](https://ilrs.gsfc.nasa.gov/docs/2024/ILRS_Data_Handling_File_2024.02.13.snx),
  opening comments, model/range-bias block and data-handling legend.
- [Station reference-point eccentricities](https://ilrs.gsfc.nasa.gov/network/site_procedures/eccentricity.html).
- IERS2010§11.2 and Eq11.19, cached primary source
  `/home/dominic/Tudat/earth-station-model-work/references/relativistic-stations/iers_chapter11.{pdf,txt}`.
- [Mendes & Pavlis2004](https://doi.org/10.1029/2004GL020308), optical zenith delay,
  together with the IERS optical mapping prescription.

## Honest accuracy targets and independent evidence

The [April2024 ILRS report card](https://ilrs.gsfc.nasa.gov/network/system_performance/global_report_cards/monthly/2024/04/2024_04_Monthly_report_card.html),
Table2, gives average LAGEOS NP RMS from rapid orbit analyses over three months:
Hitotsubashi/JCET report Graz **1.4/1.4mm** and Herstmonceux **0.7/1.3mm**.
Their short-term pass-bias stabilities are **4.9/7.3mm** and **3.8/7.0mm**,
respectively. These are distinct metrics; Table1 single-shot RMS is also distinct.
They support millimetre-level scatter as a scientific aim for good stations,
but do not establish a universal raw-RMS threshold or zero mean for our exact
fixed combined SP3/station pair. The combined SP3 report shows millimetres radial
and centimetres transverse inter-centre differences. A preliminary expectation
of millimetres to a few centimetres for unadjusted residuals is an inference,
not an acceptance tolerance fitted to results.

Report per station/pass and pooled: count, signed mean, uncentred RMS,
mean-centred standard deviation, median, and clearly defined percentile/outlier
counts. Retain the raw residuals and each correction. No pass/time/range bias
fit, demeaning, linear detrending or sigma clipping may be used to claim the
headline fixed-orbit accuracy. An independently prescribed20-degree cutoff is
available from the ASI model description; record all excluded observations.

The strongest exact dated published case found is
[Dequal et al.,2020, Matera100kHz demonstration](https://doi.org/10.1007/s00190-020-01469-2),
§3/Fig7: global fitted LAGEOS orbits over23–30January2018, a LAGEOS-2 pass with
**1.3mm raw one-way residual RMS**, but LAGEOS-1 at centimetres because time/range
biases were not optimized. Fig8's1.2–2.0mm values follow a fitted line and are
**not** a no-tuning benchmark. Signature/slow spin and atmospheric structure
limit precision; the paper's approximately0.1mm internal polynomial statistic
is not fixed-orbit absolute accuracy. The exact author-generated GEODYN orbit
and per-normal-point reference residuals are not publicly supplied by that
paper; the data-availability statement offers additional data on request.

## Secondary2018 archive and explicit limitations

`../data/slr/2018-01/` contains actual25JanuaryCRDv1 LAGEOS-1/2 files, ASI and
combined weekly LAGEOS-1/2 SP3, one NSGF orbit, combination summary and matching
ASI weekly SINEX. Headers give UTC and120s spacing over21–27January. The raw
daily files inspected do not contain the paper's special Matera7941 sessions;
do not call them the exact paper reproduction. In addition, ASI and combined
SP3 comments say SLRF2008/SLR08 while the matching dated ASI SINEX saysITRF2014
and official v70 product documentation saysSLRF2014. This is likely stale header
text, but not silently resolved here. Keep this case secondary until the exact
orbit-frame realization and historical CoM/DHF versions are established.

The pre-existing local ILR file
`/home/dominic/Software/ILR/trunk/ILR/ILR/InputOutput/UnitTests/lageos.sp3`
is a DEOS fitted2002GPS-time/ITR08 product; its original publication provenance
and matching observations are unestablished. It is useful parser history, not
the new scientific truth dataset.

## Source-reading reuse notes; no external software run

GROOPS source cache revision `7e1bd6dc8bf2fb1c9180f2f9366efa7e7b4477ee`:

- `source/slr/slrObservation.cpp:20` iterates both legs and station motion;
  `:185` reduces to half the path sum, `:191` applies Earth Shapiro.
- `source/programs/conversion/slr/crd2NormalPoints.cpp:146` converts event flags
  to transmit time, UTC to its internal GPS time and TOF to half-range. The
  converter only reads a subset of H4 and is not a substitute for auditing
  correction-applied flags in the selected products.
- `source/programs/conversion/slr/slrComModel2RangeBiasStationSatellite.cpp:95`
  uses published station/date/wavelength CoM tables. Reuse the convention,
  not a fitted-bias workflow.

Orekit has typed CRD headers/supplements and a dedicated optical
`MendesPavlisModel`; useful source/API inspection references are
[CRD.RangeSupplement](https://www.orekit.org/site-orekit-12.0.2/apidocs/org/orekit/files/ilrs/CRD.RangeSupplement.html)
and [MendesPavlisModel](https://www.orekit.org/site-orekit-12.2/apidocs/org/orekit/models/earth/troposphere/MendesPavlisModel.html).
Neither codebase was run. Official formats/papers determine the convention.

Current Tudat `src/tudat/io/readCrdFile.cpp:123–168` retains station ID, dates
and wavelength but discards H2 system occupancy/time reference and H4
correction-applied/quality flags. This is a concrete input gap for a safe
millimetre residual script; reading the raw records in the analysis script or
a narrow reader extension must precede assuming all corrections are absent.

## Rights and remaining decisions

All selected data are public official scientific downloads; manifests retain
original URLs/hashes. No explicit permissive per-file redistribution licence
was found in SP3, CRD, SINEX or CoM-table headers. Keep originals external for
now, cite ILRS, orbit DOI and CoM authors, and resolve redistribution terms
before vendoring them. Public access is not itself a licence claim.

Before the implementation claims reproduction, settle: chosen station solution
(SLRF2020 trajectory versus provider weekly estimates), exact EOP/release and
ocean-loading model, propagation-frame/time formulation, and which fixed DHF
entries are supported. The combined and ASI orbit variants must be reported
separately. An exact author residual-vector reference is not yet available;
published report-card/paper RMS values are context, not fabricated frozen
references. The first deliverable can be a transparent fixed-orbit residual
analysis with a correction ledger, exposing remaining model differences.

## Pending pre-run precision refinement: absolute BCRS ephemeris time

Static review on 2026-10-15, before the first observation reduction, found a
larger timing limit than the script's relative-orbit-only callback bound.
`SpiceEphemeris` uses the base extended-time fallback, and
`createEphemeris.h` constructs custom callbacks as `CustomEphemeris<double,double>`.
Both therefore round arbitrary `Time` queries to double TDB. At April 2024
(~7.65e8 seconds since J2000), a double spacing is 2^-23 seconds, so a half-ULP is
59.6046 ns. At 30 km/s the absolute Earth translation error can reach 1.788 mm
at one queried event. Common translation may cancel partly between endpoints,
but a submillimetre end-to-end claim cannot assume this cancellation. This is
an a priori interface limit, not a diagnosis obtained by fitting residuals.

A minimal analysis-only remedy can reuse existing ephemeris machinery:

1. Configure Earth and Moon with `ephemeris.interpolated_spice`, supplying
   exactly representable integer TDB start/end/step values and the same SSB/J2000
   frame. `create_body_ephemeris` and `create_system_of_bodies` use
   `<double,Time>`; `createStateInterpolatorFromSpice` samples SPICE only at the
   chosen integer knots and constructs a `Time`-keyed interpolator. The override
   `TabulatedCartesianEphemeris<double,Time>::getCartesianStateFromExtendedTime`
   interpolates the original precise query directly. No extended-propagation
   build or new time wrapper is required.
2. Sample the already specified SP3-to-BCRS transformation at integer TDB knots
   surrounding the selected observation arcs, before residual calculation.
   Feed that state map to `ephemeris.tabulated` and
   `environment_setup.create_body_ephemeris`; the ordinary map's exactly
   representable double integer keys are copied to `Time` in the native factory.
   Runtime range iterations then bypass the custom callback's narrowing.
   A predeclared 10-second spacing and 60-second padding around each observed
   session are a plausible initial choice, with an independent 5-second
   comparison/withheld-integer-knot check. Keep each query inside the padded
   interval, including the range solver's light-time extent; do not interpolate
   a scientific query across the gaps between observation arcs. This is a small
   analysis data preparation step, not an ephemeris framework.
3. Check tabulation independently of observed residuals: compare intermediate
   integer TDB epochs directly to the original source transformation, and compare
   5/10-second tables at the same precise off-knot query epochs. Retain the
   independent SP3 source-node holdouts and 8/10-point source interpolation checks.
   The secondary table must not be chosen to minimise the measurement residuals.

Important existing API limits: **do not use `tabulated_from_existing`** for this
precision route. Its current factory invokes `getTabulatedEphemeris` without
explicit template arguments, which defaults to double time. `ephemeris.tabulated`
hardcodes six-point Lagrange interpolation but preserves native `Time` queries;
`interpolated_spice` accepts an explicit existing interpolation setting. The
body-level `state_in_base_frame_from_ephemeris(Time)` binding preserves `Time`.
The Python `TabulatedCartesianEphemeris.interpolator` property is registered only
for the double-time subclass and is not needed for this route.

A further generic Lagrange limitation must remain explicit. It recognises a node
within `10*epsilon(ScalarType)*abs(epoch)`. For a `Time` interpolator the scalar is
long double: at this epoch the threshold is about 0.83 ns on x86 extended long
double (~25 micrometres at 30 km/s), but about 1.7 microseconds on platforms where
long double equals double (~5 cm) extremely near an absolute-time node. Earth/Moon
`interpolated_spice` can select the existing cubic spline interpolator to avoid
that node check. The fixed-SP3 `tabulated` factory does not expose its interpolator
choice. No generic interpolation fix is part of this analysis proposal; do not
claim unconditional cross-platform nanosecond response. The initial run remains
useful as a separately labelled double-callback baseline.

Double BCRS state storage also has a remaining position ULP of roughly 30
micrometres at an AU, and the Python rotation matrix API is double-time (Earth
rotation times the half-ULP is of order tens of micrometres for these radii).
Those limits are smaller than the original absolute-time translation error but
still belong in the precision ledger. No source changes or kernel execution
were performed for this review.

Further static check of the existing spline installation route: Python registers
`environment.TabulatedEphemeris` only for native `<double,double>`, with a
`.interpolator` property (there is no `reset_interpolator` Python method).
Consequently the already available
`create_one_dimensional_vector_interpolator_time_object` can create a Time-valued
cubic/Hermite spline, but cannot currently install it into the `<double,Time>`
ephemerides returned by the configured factories. An optional, small later binding
correction could register that existing class/property against configured
`STATE_SCALAR_TYPE`/`TIME_TYPE`, without adding a native evaluator or API. It must
wait until after the initial run and explicit coordination with the build owner.
A cubic spline is simpler than Hermite here: all six state components can be
interpolated without inventing missing acceleration records for Hermite's
six-component derivative input.

Table spacing must be set from an interpolation budget, not assumed harmless.
A simple LAGEOS circular-frequency estimate, radius ~1.2e7 m and angular frequency
~4.65e-4 rad/s, gives `r*omega^4*h^4/384` about 1.2 mm for a 30-second cubic
midpoint, about 0.07 mm at 15 seconds, and about 0.014 mm at 10 seconds (the
usual rigorous cubic error constant may be larger). These are order-of-magnitude
predictions, not measured errors or adopted tolerances. A 10/5-second source-state
holdout/convergence comparison near observed arcs is therefore preferable to
claiming that a full-week 30-second cubic table already supplies submillimetre
accuracy. Earth/Moon barycentric tables have much lower orbital frequencies;
choose and validate their step independently. No residuals were used in these
spacing proposals.
