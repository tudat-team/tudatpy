# Fixed-orbit LAGEOS-1 observation reduction

`fixed_orbit_slr.py` defaults to the frozen week 2024-03-31 through 2024-04-06.
The separately frozen consecutive week 2024-04-07 through 2024-04-13 is selected
with `--week-ending 2024-04-13` and its own data root.
It uses actual ILRS CRDv2 normal points and a final fitted SP3 orbit. It does not
propagate, estimate, fit an offset, remove a trend, clip residuals, or use a CPF
prediction as truth. Orbit and station configurations are separate experiments;
do not select one because it has a smaller residual mean.

The default combined-orbit/SLRF2020 configuration has been executed; its measured
baseline and raw results are recorded below. Run with a kernel containing the
SP3/CRD readers, optical correction, body-deformation station motion factory and
selected Earth-orientation converter property. No mean or offset is fitted.

## Inputs and provenance

The pinned research cache is `laser-ranging-work/data/slr`, next to the research
notes `research/SLR-REFERENCE-CASES.md`. The script takes its location explicitly:

```
python analysis/laser_ranging/fixed_orbit_slr.py \
  --data-root /path/to/laser-ranging-work/data/slr \
  --output /path/to/new-result-directory
```

The output directory must not exist; each configuration keeps its own artifacts.
The script reads, but does not download or modify, provider data. The adjacent
`dataset.json` records official URLs, hashes, retrieval dates and rights notes
for both frozen weeks and the optional loading product. It checks pinned
input hashes against that manifest, and records actual input hashes
including the selected packaged C04 EOP file. Provider gzip files are decompressed
unchanged into the output directory. Public availability is not an invented
redistribution licence: retain each provider's citation and rights notes in the
manifest and the test archive's `slr/README.md`. The manifest also pins the
optional HARPOS product described below; its bytes are kept outside the repository.

Required products are:

- `2024-04/ilrsa.orb.lageos1.240406.v81.sp3.gz` (default combined), or
  `asi.orb.lageos1.240406.v80.sp3.gz` with `--orbit asi`;
- all seven `2024-04/lageos1_YYYYMMDD.np2` daily observation files;
- `metadata/SLRF2020_POS+VEL_2023.10.02.snx` (default station solution), or
  `2024-04/asi.pos+eop.240406.v80.snx.gz` with `--station-solution asi-weekly`;
- `2024-04/ilrsa.pos+eop.240406.v80.snx.gz` for the selected zero XYZ telescope
  eccentricities, required to cover the entire week;
- `metadata/com6.220915.tar.gz` and
  `metadata/ILRS_Data_Handling_File_2024.02.13.snx`.

For the second week, the exact seven daily filenames are derived from the selected
week ending; the orbit is `ilrsa.orb.lageos1.240413.v81.sp3.gz` and eccentricities
come from `ilrsa.pos+eop.240413.v80.snx.gz`. Its combined postfit UTC/SLRF2020
orbit has 5040 position/velocity samples at 120-second cadence, through April 13
23:58. The same dated station, CoM and DHF releases are reused. This week supplies
the combined-orbit/SLRF2020 configuration; the ASI alternatives remain confined
to the original week. Data live in a separate root to retain the original frozen
seven-file selection. The [preflight record](results/frozen_2024-04-07_13/preflight.json)
pins the selection before residual evaluation: 58 raw Graz and 76 raw Herstmonceux
normal points across 17 selected sessions. The same 20-degree, weather-coverage
and SP3-interior criteria apply; no residual-dependent selection is introduced.

The station selection is fixed in advance: Graz 7839 / 11001S002 and
Herstmonceux 7840 / 13212S001 (19 and 10 normal points before physical selection).
Selected station solutions must be unique and have the correct DOMES identity.
The external SLRF2020 experiment uses its supplied velocities, with a Julian-year
conversion performed by the native reader. The ASI weekly experiment uses the
provider's already estimated weekly positions; missing velocities remain absent
and are explicitly reported as a position snapshot. Neither is an independent
withheld-data orbit validation. No extra PSD or fitted seasonal station terms
are supplied for these two selected monuments in this analysis.

## Coordinate and time conventions

The provider explicitly declares UTC, SLR20 and SLRF2020. Raw SP3 UTC epochs are
converted to TT using the converter owned by the selected GCRS-to-ITRS rotation.
Lagrange interpolation uses small TT differences from the first knot, SI position
and velocity, and no extrapolation. Missing/predicted orbit samples, maneuvers,
changing cadence and gaps are rejected. A common interior margin is used for all
8- and 10-point comparisons. A two-second interior guard covers each range
iteration and the velocity difference stencil. This deliberately excludes usable
source endpoints rather than extrapolating a high-order polynomial there.

The observation model operates in BCRS with J2000 axes and an SSB origin. The SP3
callback evaluates

`TT_site = TT_geocentre(TDB) - v_E dot r_GCRS / c^2`,

iterates the position-dependent offset, rotates at that retimed TT instant, and
applies the spatial terms of IERS 2010 Eq. 11.19 before adding Earth's BCRS
position. The solar potential and TT-to-TDB spatial scale match the existing
native relativistic station model. The latter omits simultaneity, so a small
analysis-local custom station displacement supplies `-delta_t * v_rot`, expressed
in terrestrial axes. It uses the nominal rigid station; retiming of tectonic and
tidal displacements, acceleration terms and higher relativistic orders are not
included. This is a documented approximation, not a relabelling of TT as TDB.
The native custom station factory accepts a position correction and pads its
velocity with zeros; that velocity is not used by the range observable.

CRD H2 code 4 denotes UTC(GPS), not GPS coordinate time; codes 3/7 are also UTC
realisations. Extracted native `Time` tags retain the original event and duplicate
normal points. Only the selected transmit-epoch event 2 is supported. UTC tags
are converted to TDB with the nominal station position for the local clock term.
Native `two_way_range` uses `time_scale_for_observable=tai_scale`, transferring
the full two-leg interval to station TAI once. Do not add a second time-scale
range bias. The native observable's default converter uses the same default C04
profile; orbit and external tag conversions use the actual rotation-owned object.

The script selects the packaged `eopc04_20_IAU2000.62-now.txt`, IAU 2006 rotation,
and the standard INPOP19a SPICE ephemeris with its constants, PCK and leap-second
kernel (all actual files are hashed). These are recorded analysis choices, not the
provider's estimated daily EOP. SP3 and stations both start in the same terrestrial
frame, making their common rotation less sensitive to this difference than an
inertial-orbit/terrestrial-station comparison. The script does not claim exact
reproduction of the provider's orbit determination.

`--earth-orientation baseline` retains the conventional IERS 2010 ocean EOP
terms and linear daily interpolation. The explicit alternative
`--earth-orientation desai-sibois-four-point` replaces those ocean terms with the
159-term [Desai and Sibois (2016)](https://doi.org/10.1002/2016JB013125)
series and uses four-point Lagrange interpolation of the same daily C04 data in
TAI. Conventional libration terms remain included. The model evaluates the CIO,
short-period EOP and TDB–TT directly, without additional interpolation tables.
The manifest pins the installed coefficient tables and the unchanged C04 bytes;
the output records the selected profile and actual hashes. This option changes
Earth orientation only; it does not select extra station deformations.

The rotation-owned converter supplies external time tags, the SP3 callback and
native station deformations. The native observable's TAI interval correction
still uses its default IERS 2010 / linear C04 converter. That distinction is
recorded explicitly for the alternative profile. These settings are available in
the current kernel; no high-fidelity residual result is claimed here yet.

The current default uses existing native cubic-spline `Time` ephemerides for
Earth, Sun and Moon on integer 30-second SPICE knots, padded one day beyond the
week. The satellite uses `environment.CustomEphemerisTimeObject`: its callback
keeps TDB and converted TT as `Time`, converting only local SP3 time differences
to double. `--body-ephemeris direct --satellite-epoch double` selects the original
frozen baseline route. Use `--satellite-epoch double` alone for a body-only change
and `--body-step 60` for a table-spacing convergence check. Routes are recorded
explicitly; no route is selected by its residual mean.

## Physical corrections and selection

- Native Sun/Earth relativistic light time and Mendes-Pavlis zenith delays with
  FCULa mapping are iterated on both legs. Failure to converge raises an error.
  The native fractional light-time tolerance is `1e-12` (about ten micrometres per
  LAGEOS leg); its existing Python argument is named `absolute_tolerance`.
- CRD pressure is converted from hPa to Pa, temperature is kelvin, and relative
  humidity is converted from percent to a fraction. Native meteorological
  conversion supplies vapour pressure. CRD record 20 origin 0 denotes measured
  values; origin 1 denotes station-interpolated values at the printed epoch,
  as defined in [CRD version 2.01, section 3.4](https://ilrs.gsfc.nasa.gov/docs/2022/crd_v2.01e3.pdf).
  Both supplied origins are retained in the session summary. Weather is linearly
  interpolated in UTC only within supplied sample coverage. Unsorted weather records are assigned to the
  day nearest the midpoint of the short session. The C0 system's actual wavelength
  is used (Graz 532 nm; Herstmonceux 532.080 nm). Native optical CO2 is 375 ppm.
- H4 flags must say refraction, CoM and amplitude corrections are unapplied,
  station calibration applied, and passive spacecraft delay unapplied. Calibration
  is retained as metadata and never reapplied. Unsupported flags fail explicitly.
- The pinned CoM table supplies Graz +250.4 mm and Herstmonceux +244.5 mm. The
  table's own Fortran uses inclusive integer-day bounds and a 50 nm wavelength
  match; this analysis uses those rules and rejects multiple active corrections.
- The dated DHF supplies fixed one-way biases Graz +3.9 mm and Herstmonceux
  -2.2 mm, subtracted from the observations. Any other active DHF instruction
  (including quarantine, exclusion, humidity, pressure or timing corrections)
  stops this bounded analysis rather than being silently ignored.
- Native IERS2010 solid and pole tides act on conventional tide-free station
  coordinates. The current native shared mean-pole convention is retained.
  No atmospheric-pressure loading is applied, matching the dated ASI description.
- **Exact GOT4.10c ocean loading is absent by default.** The summary labels this
  omission. `--got410c-blq /path/to/product.blq` accepts a user-declared compatible
  product, records its hash, and uses native HARDISP synthesis. Its station labels
  must match the selected DOMES identifiers and its reference origin must already
  match the station frame; no label/frame/origin transformation is inferred.
  Do not silently substitute another ocean model or call the default run a
  complete provider-model reconstruction.
- `--got410c-harpos /path/to/toc_got410c_harmod.hps` selects a separate, pinned
  GOT4.10c realization from the International Mass Loading Service. It is mutually
  exclusive with BLQ loading and leaves the default unchanged. Exact case-sensitive
  provider sites are `Graz` and `Herstmon`; the nearby `GRAZ` and `HERS` monuments
  are not substituted. The summary and manifest retain product hash, provider
  XYZ, separation from the selected station and the 3 km applicability radius.
  The helper uses 17 supplied TT harmonic terms, geodetic Up/East/North,
  whole-Earth centre-of-mass origin and STW105/REAR response. Up follows the product
  narrative and provider's consumer code; the basis section below documents the
  conflicting formal HARPOS wording. Ocean-pole and equilibrium
  terms are absent. This is an alternate realization, not exact ASI 11-constituent
  BLQ/HARDISP reproduction. One native custom station-motion callback supplies
  the displacement; the terminal relativistic map then sees the complete loaded
  station position. The existing nominal-only simultaneity approximation remains.
- The provider's 20-degree elevation cutoff is applied at both ground events,
  using the observed half-flight time to locate the passive bounce for selection.
  No cut is made on a computed residual. Weather and orbit coverage exclusions
  are listed separately. The final Graz session crosses the weekly orbit end.

## Outputs and precision assessment

`residuals.csv` retains every accepted point separately, exact split UTC fields,
source/session/configuration identifiers, raw TOF, CoM and DHF contributions,
source table lines, full computed two-way range, NP scatter, elevation and
8-versus-10-point terrestrial-position differences. The sign and units are

`O-C [one-way metres] = c*observed_TOF/2 + CoM - DHF_bias - computed_two_way/2`.

`summary.json` records input hashes, station coordinates/velocities, H4 flags,
calibration/supplement counts, every exclusion, and signed mean, uncentred RMS,
centred standard deviation and median per session, station and total. None of
these statistics is fitted back into the model.

Source-node holdouts omit one central sample in a local window and compare the
interpolated position to the withheld provider sample at orders 8 and 10. These
and order convergence quantify interpolation independently of SLR residuals;
the SP3 position print resolution is 1 mm. The default custom callback retains
the native observable's `Time` epoch. The explicit legacy double route reports
its largest callback half-ULP and a position bound from relative orbital/rotational
speed; zero callback truncation in the Time route is not an overall error bound.
The rotation interface remains double, and direct SPICE also rounds epochs;
none of these diagnostics establishes a submillimetre end-to-end floor.
The 0.5 s velocity stencil is only for a consistent custom state; range uses its
position. Absolute BCRS position subtraction can add tens of micrometres of
floating-point rounding.

Official April 2024 report cards give station NP precision near 0.7–1.4 mm and
pass-bias stability several mm. These describe a complete analysis, not a raw
zero-mean acceptance threshold for this SP3 reduction with omitted ocean loading.
The combined product's orbit-comparison RMS is not a measured-range RMS. Report
all results and correction omissions before interpreting any difference.

## Measured frozen baseline

The 2026-10-16 run against the unchanged combined orbit and SLRF2020 station
release accepted 23 normal points: 13 Graz and 10 Herstmonceux. Signed one-way
O−C mean was −0.474 mm, uncentred RMS 8.843 mm and centred standard deviation
8.830 mm. Graz RMS was 11.526 mm; Herstmonceux RMS was 2.675 mm. No fitted bias,
demeaning, residual clipping or parameter estimation was applied. Three further
points fell below the declared 20-degree cutoff and three lacked the common SP3
interpolation interior. Exact GOT4.10c ocean loading was omitted.

The [raw residuals](results/slrf2020_combined_2024-03-31_2024-04-06/residuals.csv)
and [summary](results/slrf2020_combined_2024-03-31_2024-04-06/summary.json) retain
every accepted point, exclusions, input/software hashes and correction metadata.
Only machine-specific paths were replaced for portability. Independent 8/10-point
SP3 node holdouts had position RMS 1.412/1.481 mm and maxima 5.818/5.986 mm.
The custom callback's 59.605 ns half-ULP gives a 0.415 mm relative-position bound;
this is not yet an end-to-end range precision measurement. Body SPICE evaluation
and the satellite custom callback retain their original double epoch boundaries
in this baseline.

An isolated follow-up replaced only the Earth/Sun/Moon direct-SPICE ephemerides
with existing native cubic-spline `Time` tables on integer 60/30-second knots.
The satellite custom callback, physical corrections, 23 accepted points and six
exclusions stayed identical. Maximum absolute computed half-range changes were
0.773/0.780 mm relative to the baseline; 60-to-30-second convergence was 0.0117 mm
maximum and 0.00461 mm RMS. These are route differences, not fitted residual
improvements. Equal-double table/direct position controls differed by at most
0.063 mm for Earth/Moon and 0.0005 mm for the Sun. The remaining satellite and
rotation double boundaries prevent an end-to-end submillimetre claim.

The [precision summary](results/body_time_precision/summary.json) and
[per-point half-range changes](results/body_time_precision/half_range_changes.csv)
retain the controls and hashes. The original experiment wrapper is retained in
source commit `91ee4afbd`; route selection now belongs to the main script.
For a controlled comparison on a current kernel, first run that script into
separate result directories with the desired explicit route options, then use:

```
python analysis/laser_ranging/body_time_precision.py \
  --baseline /path/to/current-kernel-legacy-result \
  --candidate /path/to/current-kernel-time-result \
  --output /path/to/new-comparison-directory
```

The comparison requires identical input, script and kernel hashes, physical
configuration, accepted point keys and exclusions. It reports computed half-range
changes without fitting or comparing residual improvement.

The subsequent precise-callback experiment used one fresh kernel and script for
all four explicit routes. The legacy route reproduced every original frozen CSV
field exactly. Body-only and satellite-callback changes reached 0.780 and
1.118 mm in computed half-range respectively; their combined effect reached
1.647 mm. With the Time callback, cubic 60-to-30-second convergence was 0.00959 mm
maximum and 0.00344 mm RMS. All routes retained the same 23 points and six
exclusions. The precise 30-second no-loading result had one-way O−C mean
−0.215 mm and RMS 8.768 mm; no mean or offset was fitted.

The [precise raw CSV](results/precise_time_slrf2020_combined_2024-03-31_2024-04-06/residuals.csv),
[run summary](results/precise_time_slrf2020_combined_2024-03-31_2024-04-06/summary.json),
[route changes](results/precise_time_slrf2020_combined_2024-03-31_2024-04-06/route_changes.csv)
and [precision controls](results/precise_time_slrf2020_combined_2024-03-31_2024-04-06/precision_summary.json)
retain the numerical evidence and hashes. Callback epoch truncation is zero in
the Time route. Rotation still rounds to double: its 59.605 ns half-ULP gives a
0.0537 mm satellite-position bound using the maximum source radius and Earth
rotation rate, not a measured full-range error. Absolute BCRS state arithmetic
also retains tens of micrometres of rounding. The measured table convergence,
1 mm SP3 print resolution and approximately 1.48 mm source-node holdout RMS do
not establish an end-to-end precision floor below these remaining limitations.

The first separate HARPOS experiment stopped at Herstmonceux session 29, point 0.
With the original 50-iteration, `1e-12` fractional criterion, receiver epochs
entered an exact two-cycle 31.752 fs apart: 9.519 µm of light distance, above the
7.267 µm criterion for that leg. The Earth BCRS coordinates alternated by
30.518/7.629/0 µm while the bounce epoch stayed fixed and loading evaluations
were identical. No tolerance, physical cut or acceptance rule was changed.
The same-script no-loading control reproduced every precise-baseline CSV field.
The [portable failure record](results/got410c_harpos_radial_failure/failure.json)
preserves exact epoch fractions, states, selected point and input/software hashes.
This failed run used the original radial-Up helper; the subsequent provider
clarification to geodetic Up is a separate correction. No complete ocean-loading
residual result is claimed from this attempt.

The [same-polynomial proof](results/got410c_harpos_radial_failure/cubic_same_polynomial_proof.json)
reconstructs the original cubic output exactly and independently solves the
natural spline at 80-digit precision. Both failed receiver epochs round to the
same reference position; their actual separation is below 1 nm in each component
and does not cross a double rounding midpoint. Evaluating the same polynomial
from the nearer endpoint removes the artificial large-offset weighted-sum jitter
while retaining both exact endpoints. A five-node analytic-motion regression in
the existing cubic suite exercises this numerical issue without changing any
light-time tolerance or physical model.

## Stable-cubic results for both frozen weeks

The corrected cubic evaluation let the original failing radial-basis point and
all declared full-week runs converge with the original 50-iteration, `1e-12`
criterion. First-week accepted keys and all six exclusions stayed unchanged;
the historical and date-parameterized no-loading runs matched every CSV field
on the same kernel. The generic arithmetic change shifted no-loading computed
half-ranges by at most 0.00959 mm relative to the preceding kernel.

| Week ending | Accepted / raw | No-loading mean / RMS (mm) | Geodetic HARPOS mean / RMS (mm) |
| --- | ---: | ---: | ---: |
| 2024-04-06 | 23 / 29 | −0.214 / 8.769 | −0.057 / 7.221 |
| 2024-04-13 | 131 / 134 | +5.258 / 18.752 | +7.659 / 17.435 |

All statistics are signed one-way O−C with uncentred RMS. The second week retains
55 Graz and 76 Herstmonceux points; three Graz points fall below the unchanged
20-degree cutoff in both loading variants. The historical radial first-week
control had mean −0.052 mm and RMS 7.222 mm. Provider source evidence determined
the geodetic basis. No fitting, demeaning, clipping or residual-based model choice
was applied; no loading remains the default.

The [compact combined summary](results/stable_cubic_two_weeks/summary.json) links
five portable raw CSVs, computed half-range changes, actual input/software hashes
and the preserved failure/proof. HARPOS changes computed half-range by up to
6.737 mm in week one and 14.007 mm in week two. This alternate CM/STW105
realization still differs from the provider's exact BLQ/HARDISP reduction.

## GOT4.10c loading basis

`ocean_loading.py` evaluates the pinned 17-constituent GOT4.10c HARPOS product,
including its supplied nodal sidebands, in TT seconds since J2000. Its callback
uses the selected Earth-rotation converter at the geocentre and returns terrestrial
XYZ metres. Exact case-sensitive sites are `Graz`, `Herstmon` and `Apollo`;
`GRAZ` and `HERS` are different nearby monuments. The nominal station must lie
within the provider's 3000 m applicability radius.

Up follows the product narrative and NASA's actual HARPOS consumer: it is normal
to the reference ellipsoid. Tudat's existing Cartesian-to-geodetic conversion
uses the consumer's `REF_ELL` constants, `a = 6378137 m`, `f = 1/298.257`, at
the provider XYZ. This is not the exact WGS84 flattening. The pinned consumer
`MAKE_UEN_TO_XYZ` disables its radial branch and constructs Up/East/North from
geodetic latitude; `VTD_SET_HARPOS` applies that matrix to the amplitudes.
The [basis provenance](ocean_loading_basis.json) links the primary source files
at the audited NASA revision.

The formal HARPOS specification instead describes radial Up. The exact 2017
generator binary/source has not been recovered, so the current consumer supports
this choice without proving the historical generator's basis. This change follows
source evidence, without selecting a basis by residual size. Existing radial
result files remain unchanged as historical controls; their metadata identifies
the former radial basis. They must not be relabelled as geodetic results.

Only the 17 supplied gravity-tide terms are included; ocean-pole and equilibrium
terms are absent. Do not duplicate these gravity tides with a BLQ realization.
The STW105/REAR whole-Earth-centre-of-mass realization differs from ASI's
11-constituent BLQ plus HARDISP and JPL's FES2004 loading, so it is not an exact
replay of either provider's fitted reduction. The frozen baseline above still
omits ocean loading.
