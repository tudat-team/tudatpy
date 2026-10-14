# Independent short-period EOP benchmarks

`expected.csv` was produced by executing unmodified, checksum-pinned IERS
reference software. No Tudat output is an expected value. The 35 UTC MJD epochs
span 1962–2025, published validation epochs, negative J2000 times, midnight and
fractional days, and samples over an 18.6-year nodal cycle. Normal C++ tests are
offline; the generated header avoids runtime file access for expected vectors.

## Sources, versions and licensing

`sources.json` pins exact URLs and SHA-256 values for every compiled reference
source, the relevant chapter PDF, and the ocean-table routine's explanatory
readme. Sources were retrieved on 14 October 2026. Downloaded sources and binaries
stay in an external local cache; they are not redistributed in this repository.

- IERS Chapter 5 software and correction history:
  https://iers-conventions.obspm.fr/chapter5.php
- `PMSDNUT2.F`: 13 October 2011 revision, announced 14 October 2011. Its corrected
  published check at MJD 54335 is `[24.83144238273364834,-14.09240692041837661]`
  microarcseconds. `IBAND=1` includes only ten quasi-diurnal terms; fifteen long
  periods and secular drift are deliberately excluded.
- `UTLIBR.F`: 23 June 2010 coefficient revision, eleven semidiurnal terms. Outputs
  are incremental UT1 and excess LOD, in microseconds. The header calls LOD units
  “microseconds per day”; this is a day-length increment, converted to seconds,
  not radians/second or a UT1 clock drift rate.
- `FUNDARG.F`: IERS lunisolar fundamental arguments based on Simon et al. (1994).
- IERS Chapter 8 `ORTHO_EOP.F`/`CNMTX.F`:
  https://iers-conventions.obspm.fr/chapter8.php
  `ORTHO_EOP` includes its 19 March 2010 Fortran compatibility revision; `CNMTX`
  includes its 17 March 2010 table-layout revision. They return oceanic
  `[xp,yp,UT1]` in microarcseconds/microseconds using an independent ortho-weight
  formulation.
- Matching 71-wave ocean-table implementation, Observatoire de Paris:
  https://hpiers.obspm.fr/eoppc/models/interp.f
  The reference calls only `PMUT1_OCEANS`, not its interpolation or older PM_GRAVI
  routines. Its header credits Ch. Bizouard (2002), McCarthy and D. Gambis (1997);
  the file identifies a September 2007 correction and was listed by the source
  server as updated 14 June 2022. This supplied file has no embedded license
  notice; it is executed locally and neither modified nor redistributed.
- `DAT.F` and `CAL2JD.F`: unmodified SOFA dependencies distributed with IERS Chapter
  7, including the 2017 leap second. They establish TT epochs independently of
  Tudat. Their post-version warning is accepted for later fixtures; the supplied
  leap-second table gives TAI−UTC=37 seconds for the 2017–2025 fixture dates.

The results acknowledge software supplied by the IERS Conventions Center,
A. Brzezinski's libration models, R. Eanes's ocean model, and the Paris Observatory
matching-table implementation. `LICENSE_IERS.txt` reproduces the IERS notice.
The reference driver, generator and fixtures are Tudat materials, neither
supplied nor endorsed by the Center. No original reference routine is modified.
The separately named driver combines calls and reproduces their documented GMST
polynomial solely to export independent arguments. Reference source notices
remain intact in the external cache; SOFA dependencies retain their own notices.
The small UT1/LOD coefficient extracts describe their IERS origin in the generator.

## Reproduction

Use Python 3 and GNU Fortran, keeping all sources/build products outside Git:

```sh
python3 generate.py --reference-dir /absolute/local/reference-cache --download
python3 generate.py --reference-dir /absolute/local/reference-cache --check
```

Use `--compiler /absolute/path/to/gfortran` if required. The fixtures were generated
with GNU Fortran 15.2.0 (conda-forge build 15.2.0-19), flags
`-O0 -ffixed-line-length-none`, one compiler process. Checksums are verified before
compilation. Published PMSDNUT2, UTLIBR and ORTHO checks must reproduce before any
artifact is written. `--check` verifies the decimal CSV, C++ table and extracted
UT1/LOD coefficients byte for byte. Compiler-dependent last-digit differences
should be assessed numerically against those independently published checks.

## Explicit time, phase and coefficient conventions

The original programs take a single MJD; they use that same epoch for the old
GMST polynomial and lunisolar centuries. For these benchmarks that input is
explicitly UTC MJD, with UT1=UTC as a controlled approximation. TT epochs are
computed independently using the supplied SOFA DAT routine and TT=TAI+32.184 s.
The production default uses TT lunisolar arguments and UTC approximately UT1 for
GMST06, with its existing TDB approximately TT API convention. Tests distinguish
strict coefficient/phase checks using exported original arguments from full
production-default checks with these physical time scales. They do not claim the
single-MJD reference programs independently represent an exact two-scale model.

Arguments are `[GMST+pi,l,l-prime,F,D,Omega]`, in radians, and each harmonic is
`sine_amplitude*sin(argument) + cosine_amplitude*cos(argument)`. Polar motion
amplitude tables use microarcseconds; UT1 tables use microseconds. No additional
Doodson phase or sign changes are applied. The matching ocean routine's original
amplitudes are REAL*4; its solar-anomaly cubic term differs in sign from FUNDARG.
Those tiny reference details are retained, not silently patched.

The strict libration checks use 1e-5 microarcseconds for PM, 1 ps for UT1, and
10 ps for original UTLIBR LOD. Matching ocean-table checks use 0.001
microarcseconds and 0.1 ns UT1 to cover its original REAL*4 and phase-polynomial
choices. Independent SOFA/Fortran diagnostics found maximum production-default
versus single-MJD differences of 0.0067 microarcseconds for libration PM,
0.98 ns libration UT1 after the explicit coefficient delta below, 0.148
microarcseconds ocean PM and 20.5 ns ocean UT1. Default checks use 0.01/0.2
microarcseconds and 2/30 ns, with sums checked against independent summed oracles.
The alternative ortho-weight formulation differs by up to 1.64 microarcseconds
and 107 ns over this epoch set; its separate 2-microarcsecond/120-ns checks
supplement the tighter matching-table oracle.

### Published Table 5.1b versus UTLIBR source inconsistency

The accepted Chapter 5 Table 5.1b (PDF page 11, zero-based page 10) gives the T2
wave `[2,0,-1,-2,2,-2]` UT1 sine coefficient as **0.05 microseconds**.
The downloaded unmodified UTLIBR code uses **0.04**. Tudat's accepted table remains
unchanged. The C++ test explicitly requires

`Tudat_UT1 - original_UTLIBR = 0.01 microseconds * sin(T2_phase)`

using the independently generated original arguments and a strict residual
bound. This prevents a permissive tolerance from hiding a coefficient mismatch.
The PDF also gives N2 LOD cosine −4.2 microseconds where the reference code uses
−4.1. The LOD test is labelled as the original UTLIBR convention and uses its
original coefficient extract; there is no claim that a production LOD correction
API implements the table's alternate coefficient.

## Composition and cutoff regressions

The default correction models include only ten quasi-diurnal PM libration and
71 diurnal/semidiurnal ocean terms. Tests reject long-period/secular terms by
checking their independent multiplier scope, then add a controlled daily EOP
baseline exactly once in the real PolarMotionCalculator path. Zonal long-period
UT1 corrections are not included in these short-period calculators.

A synthetic observed CPO containing an FCN-like 430-day signal is applied through
the production PrecessionNutationCalculator and compared to direct SOFA Xys06a
plus the observed offsets. This checks that no extra FCN signal is added when
CPO already contains it. It does not claim to validate a separate empirical FCN
model, which this path does not call.

Positive-amplitude cutoff tests exercise scalar UT1 and two-component PM, known
retained waves, all-filtered outputs, and empty matrices. They protect the
production fix to a reversed amplitude/multiplier argument order that previously
caused filtering to use integer multipliers and invalid Eigen block dimensions.
