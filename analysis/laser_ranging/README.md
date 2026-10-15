# Fixed-orbit LAGEOS-1 observation reduction

`fixed_orbit_slr.py` is a bounded analysis of 2024-03-31 through 2024-04-06.
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
The script reads, but does not download or modify, provider data. The adjacent `dataset.json` records official URLs, hashes, retrieval dates and
rights notes for all fourteen selected provider products. It checks pinned
input hashes against that manifest, and records actual input hashes
including the selected packaged C04 EOP file. Provider gzip files are decompressed
unchanged into the output directory. Public availability is not an invented
redistribution licence: retain each provider's citation and rights notes in the
manifest and the test archive's `slr/README.md`.

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
the SP3 position print resolution is 1 mm. The custom ephemeris API accepts double
TDB, while native observable epochs use `Time`. The largest callback half-ULP
and an upper position bound from relative orbital/rotational speed are reported.
The full native SPICE and rotation paths also have double interfaces; none of
these diagnostics claims a submillimetre end-to-end floor without runtime checks.
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
