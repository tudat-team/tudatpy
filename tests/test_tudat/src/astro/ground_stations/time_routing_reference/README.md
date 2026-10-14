These clock inputs are generated independently of Tudat by the original C driver
`sofa_clock_reference.c`, linked directly with unmodified SOFA release 2023-10-11.
It scans calendar 2020 at 6-hour intervals and refines each seasonal TT/TDB offset
extremum by ternary subdivision. The two-part TT dates and offset literals are
used in `unitTestEarthDeformationTimeRouting.cpp`; those dates are reference
inputs, rather than an assertion that the extrema are unique across all years.
Leap inputs independently exercise SOFA UTC -> TAI -> TT.

Reproduce with one compiler process, directing the binary outside the checkout:

```
cc -I "$SOFA_PREFIX/include" sofa_clock_reference.c "$SOFA_PREFIX/lib/libsofa.a" -lm -o /tmp/sofa_clock_reference
/tmp/sofa_clock_reference > /tmp/sofa_clock_reference.txt
diff -u expected.txt /tmp/sofa_clock_reference.txt
```

SOFA source is not vendored or modified. The reference library was installed as
`sofa-cmake-1.0.1-heb58013_0` (IAU release is given by its public header):
https://conda.anaconda.org/tudat-team/linux-64/sofa-cmake-1.0.1-heb58013_0.conda
Package SHA256: cc97c9a22dd6d931d90525a6e1696fe0f90dc24b3b88a89ff8a2260ffd715322
`lib/libsofa.a` SHA256: ca7408bcb8aabc6e5ca21e17b4a51ada37fea66c24066b3103b9a5954a8772f3
`include/sofa/sofa.h` SHA256: 78765fb63762ca64dbeaf277202d19216d57a406267ac1f64769d11ea42369d2

The library carries the IAU SOFA Software License, copyright 2023 IAU SOFA Board;
its actual license is in the installed header, regardless of the package's
coarser BSD metadata. Only our driver and numerical outputs are retained here.
See https://www.iausofa.org/2023-10-11c and the SOFA time-scale cookbook
https://www.iausofa.org/s/sofa_ts_c.pdf section 3.6.5 for the geocentric TDB/TT
model and approximately 1.7 ms annual offset. The driver uses zero observer
coordinates: UT1 does not enter its geocentric result.

The routed displacement tests reuse the original IERS DEHANTTIDEINEL and
HARDISP fixtures with their existing pinned provenance. No expected displacement
is generated from a Tudat production model. Mixed-state tests separately check
convergence of complete routed displacement derivatives; they are consistency
checks rather than an external scientific displacement oracle.
