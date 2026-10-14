# Independent four-point EOP interpolation reference

`expected.csv` and `expected.h` are outputs from the unmodified Paris Observatory
`LAGINT` routine in `interp.f`, with unmodified SOFA `IAU_DAT` and `IAU_CAL2JD`
distributed by the IERS Conventions Center. No Tudat result supplies an expected
value. The original reference sources remain in an external cache and are neither
modified nor redistributed here.

`sources.json` pins the URLs and SHA-256 checksums retrieved on 14 October 2026:

- https://hpiers.obspm.fr/eoppc/models/interp.f
- https://hpiers.obspm.fr/eoppc/models/interp.readme
- https://iers-conventions.obspm.fr/content/chapter7/software/dehanttideinel/DAT.F
- https://iers-conventions.obspm.fr/content/chapter7/software/dehanttideinel/CAL2JD.F

The Paris readme identifies the revised IERS Gazette No. 13, 29 January 1997,
by Dennis McCarthy, Daniel Gambis and Ch. Bizouard. The routine credits Bizouard's
November 2002 implementation and September 2007 correction. The downloaded
Paris files have no embedded license notice; they are executed locally without
redistribution. The SOFA sources retain their own license notices in the external
cache. The original driver, generator and small controlled input/output fixtures
are Tudat materials and are not supplied or endorsed by those institutions.

The three controlled eight-sample datasets cover ordinary daily epochs in 2020,
the December 2016 inserted second, and the 1962 UTC drift era. The generator
specifies the raw UTC MJD inputs, polar motion and CPO in arcseconds, DUT1 in
seconds, and LOD in seconds. It calls original SOFA DAT for the sample and query
epochs, forms continuous relative TAI seconds, subtracts DAT from DUT1, and passes
those arrays to original LAGINT. Thus the oracle independently validates the
accepted four-point algorithm on the model's continuous time axis. It does not
call `INTERP`, whose combined correction and old scalar-UTC handling are outside
this interpolation benchmark.

The reference uses the four points surrounding the bracketing interval, with
the window clamped to the first/last four samples near coverage boundaries.
Affine changes of its independent coordinate preserve the Lagrange polynomial;
the driver uses relative TAI seconds to retain leap and drift timing precision.
The original routine leaves its window index undefined when the query is exactly
the final sample. That branch is never executed for a reference output. C++ tests
check the final sample directly against the raw fixture input and label this
separate provenance. No silent patch is applied to the original reference.

Reproduce with Python 3 and GNU Fortran, using an external cache:

```sh
python3 generate.py --reference-dir /absolute/reference-cache --download
python3 generate.py --reference-dir /absolute/reference-cache --check
```

Use `--compiler /absolute/path/to/gfortran` if needed. Generation used GNU
Fortran 15.2.0 (conda-forge 15.2.0-19), `-O0 -ffixed-line-length-none`, and one
compiler process. All source checksums are verified before execution. Normal
C++ tests use the small committed header offline. Independent tests of a known
cubic polynomial validate analytic slopes separately from this value oracle.
Another test deliberately has increasing raw UT1 samples but a negative interior
cubic derivative, verifying that the model rejects an invalid inverse map.
