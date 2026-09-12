# Eros catalogue-debiasing sign experiment

These fixtures test the optical reader's catalogue-correction sign against real
MPC observations of (433) Eros and a frozen JPL Horizons reference orbit. They
contain **111 observations in two intervals**. No orbit is fitted and no equations
of motion are integrated. Every selected row is retained, including outliers.

| Case | UTC interval | Station | Catalogue code | Measurements |
| --- | --- | --- | --- | ---: |
| `eros_711_1996` | 1996-04-23 – 1996-06-25 | 711, McDonald Observatory, Fort Davis | i (GSC 1.1) | 20 |
| `eros_t05_2016` | 2016-06-13 – 2016-12-22 | T05, ATLAS-HKO, Haleakala | L (2MASS) | 91 |

Selection keeps all CCD observations with a supported catalogue code for each
station/calendar year. The intervals were selected after screening for correction
size and astrometric scatter. This is a targeted regression, not a blind assessment
of the debiasing scheme's performance. A correction derived from average catalogue
errors need not improve every observation or series. Horizons may have assimilated
these same MPC observations; its orbit is an external reference, not an independent
holdout solution.

The first run through the PR kernel also exposed an ancillary-conversion defect:
metadata-only input produced an empty, non-null simulation-settings object, which
the angular model rejects. The accompanying factory fix checks the converted
settings rather than the input metadata maps. Both the small optical regression
and this complete residual calculation cover that path.

## Reproduce

Build the `kernel` target from the PR branch and import that build together with
the branch's Python sources. Verify `tudatpy.kernel.__file__` before running:

```bash
cmake --build build --target kernel --parallel 4
python tests/test_tudatpy/validation/eros_catalogue_debiasing.py --output /tmp/eros-check
python -m pytest -q tests/test_tudatpy/test_data_eros_catalogue_debiasing.py
```

The validation script writes residual CSVs, two-panel PDF plots, and a JSON summary
with RMS values, model sensitivities, and SHA-256 hashes of the kernel, EOP data,
and bias map. Tests and residual calculation make no network requests. They need
the usual Tudat SPICE/Earth-orientation resources and
`~/.tudat/resource/star_catalog_biases/debias_2018/bias.dat`. The Python environment
also needs NumPy, pandas, Astropy, astroquery, astropy-healpix, and (for plots) matplotlib.

To explicitly obtain fresh fixtures in a separate directory:

```bash
python tests/test_tudatpy/validation/download_eros_catalogue_fixtures.py --output /tmp/eros-fixtures
```

That script queries MPC and Horizons and independently samples the JPL bias map.
New service responses may change as observations or the Horizons solution change;
inspect differences before replacing the frozen fixtures. Per-case `metadata.json`
records acquisition time, station coordinates, selection, and conventions.

## Observation model and conventions

- Eros, Earth, Sun, and Moon translations all come from Horizons geometric vectors
  relative to SSB, in the ICRF/J2000 equatorial frame. The NPZ files store TDB seconds
  since J2000 and SI states, hourly with at least two days of boundary padding.
  Tudat uses its tabulated ephemeris interpolation. The script compares hourly and
  two-hour grids to measure interpolation sensitivity.
- Angular observations use retarded Eros positions at Earth-station reception
  epochs, with converged light time and first-order solar Shapiro delay. The CSVs
  also contain Horizons quantity 1, astrometric ICRF RA/DEC, queried directly at each
  MPC station and UTC epoch with extra precision. These are compared with Tudat's
  predictions. Apparent of-date coordinates, stellar aberration, and atmospheric
  refraction are inappropriate for this catalogue-frame comparison.
- Earth uses WGS84 geometry and IAU 2006 GCRS-to-ITRS rotation, including Tudat's
  default IERS C04 Earth-orientation corrections (UT1, polar motion, celestial-pole
  offsets and subdaily rotation terms). The available daily series covers both
  selected intervals. UTC-to-TDB conversion is performed by the tracking-data
  collection factory with the station position.
- Station ITRS positions are calculated directly from MPC longitude and geocentric
  parallax constants, with equatorial radius 6378137 m. IERS 2010 solid Earth tides
  and pole tides are applied through station deformation. **Ocean loading tides
  are omitted.** Subdaily ocean terms in Earth rotation are a separate model.
  MPC coordinates do not provide a surveyed velocity history; the script therefore
  measures the angular sensitivity to a 100 m Cartesian offset instead of assuming
  exact historical station positions.
- Residuals are observed minus computed, in arcseconds. RA differences are wrapped
  at the meridian and multiplied by cos(observed DEC). Combined RMS is
  `sqrt(mean((delta_RA*cos(DEC))**2 + delta_DEC**2))`, not the per-component RMS.
  Raw, PR-corrected, and reversed-sign cases use exactly the same predictions.

## Results with the rebuilt PR kernel

| Case | Raw RMS | PR subtraction RMS | Reversed sign RMS | Correction/raw RMS |
| --- | ---: | ---: | ---: | ---: |
| 711/1996 | 0.725733 | 0.414097 | 1.277745 | 84.4% |
| T05/2016 | 0.219203 | 0.139551 | 0.320136 | 51.9% |

RMS values are arcseconds. Maximum Tudat/Horizons observer-angle disagreement was
0.0000141 arcsec. Changing the ephemeris grid to two hours changed predictions by
at most 0.00000112 arcsec; omitting solid/pole tides or solar Shapiro delay changed
them by less than 0.00000034 arcsec. A 100 m station offset changed predictions by
at most 0.000103 arcsec. These sensitivity checks are much smaller than the
catalogue corrections; they are not a claim of absolute station-coordinate accuracy.

The tests require agreement with Horizons within 1 milliarcsecond, a correction
exceeding 40% of raw RMS, at least 20% improvement with the reader's correction,
and at least 30% worsening with its reversed sign. They also check map values,
row retention, and identical predictions between the three correction cases.

The recorded run used EOP SHA-256
`cc838b27572e89ff6d36ee7cbeb6a8f104f1f8e59ddf769e8c93e6515276cd32`
and bias-map SHA-256
`ef6a50830bb83d8b1161e7acb391eed988ea9708e9fbcccc0c834da6109476e3`.

## Older observations and alternative intervals

Screening used the 18,506 available MPC Eros observations retrieved on 2026-09-12.
The table below reports direct Horizons observer-table screening, before the
separate Tudat model comparison. Units are arcseconds; no residual-based row
rejection was applied within these candidate series.

| Station/year | N | Catalogue | Raw RMS | Subtract bias | Add bias | Bias RMS |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| 532/1898 | 84 | p | 3.513 | 3.634 | 3.410 | 0.275 |
| 000/1900 | 139 | m | 0.393 | 0.486 | 0.366 | 0.176 |
| 608/1995 | 11 | c | 0.884 | 0.935 | 0.967 | 0.351 |
| **711/1996** | **20** | **i** | **0.726** | **0.414** | **1.278** | **0.613** |
| A26/2003 | 29 | c | 0.724 | 0.596 | 1.041 | 0.442 |
| C95/2012 | 158 | u | 0.241 | 0.233 | 0.271 | 0.077 |
| G45/2016 | 84 | R | 0.444 | 0.389 | 0.507 | 0.085 |
| H21/2016 | 119 | u | 0.204 | 0.193 | 0.328 | 0.175 |
| **T05/2016** | **91** | **L** | **0.219** | **0.140** | **0.320** | **0.114** |

The 1996 interval has the largest correction relative to raw residual RMS among
these useful candidates (84%). The 2016 T05 interval offers lower residual scatter
and a 52% correction fraction. A26/2003 and H21/2016 are possible additional Eros
cases, but give weaker improvement. The 1898/1900 catalogue-coded observations
explicitly do **not** favour subtraction. They also predate the available daily
1962-onward Earth-orientation series and require long proper-motion extrapolations,
so they are documented as counterexamples rather than asserted as high-accuracy
sign regressions.

## Sources

- [MPC observations via astroquery](https://astroquery.readthedocs.io/en/latest/mpc/mpc.html)
  and [MPC observatory codes](https://docs.minorplanetcenter.net/mpc-ops-docs/apis/obscodes/).
- [JPL Horizons manual](https://ssd.jpl.nasa.gov/horizons/manual.html), vector tables
  and observer quantity 1. Horizons quantity 36 formal uncertainties are retained
  in the CSVs; they describe the orbit prediction, not MPC measurement noise.
- [JPL debiasing archive](https://ssd.jpl.nasa.gov/ftp/ssd/debias/debias_2018.tgz),
  `bias.dat` and its README: position/proper-motion biases relative to Gaia DR2,
  NSIDE 64 RING. The README prescribes subtracting the epoch-adjusted biases.
  The CSV bias columns are independent direct map samples, including RA*cos(DEC).
- Tudat examples `estimation/improved_estimation_with_mpc.py` and
  `estimation/eros_space_astrometry_residuals.py` informed the setup. The latter's
  Horizons/residual approach is used without its optional propagation machinery.
