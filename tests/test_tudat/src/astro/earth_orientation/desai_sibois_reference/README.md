# Desai–Sibois ocean-tide EOP reference

The production coefficient table embeds the exact NASA IVS ASCII mirror of the
Desai and Sibois (2016) supplement, DOI `10.1002/2016JB013125`. The pinned table
has 159 lines: 86 diurnal and 73 semidiurnal. `sources.json` identifies the
original table, paper and NASA working-group memo, including byte counts and
SHA-256 hashes. The 2017 memo's model-format table identifies gamma as GMST+pi,
the next five arguments as the Delaunay arguments, and xp/yp and UT1/LOD units
as microarcseconds and microseconds. Its model-composition paragraph explicitly
retains IERS libration alongside this ocean model.

`generate.py` reads the original pinned table, checks its hash, and performs an
independent scalar harmonic sum. It imports no Tudat modules and uses no Tudat
reader, summation engine or rate routine. Fifty-digit Decimal arithmetic
evaluates the IERS 2003 Delaunay polynomials and the SOFA ERA00/GMST06 equations;
only the reduced phases are converted to double for Python's sine/cosine.
Analytic polynomial derivatives supply the independent harmonic-rate oracle.
The 14 epochs span 1962–2025, including four fractions of a day near J2000 and
several phases across a lunar nodal cycle. TT and UT1 are explicitly separate,
with UT1=TT-69.584 seconds and dUT1/dTT=0.99999998, a controlled clock fixture
that does not claim historical observed UT1. Values are in µas, µas, µs, µs;
rates are in those units per SI second.

Regenerate `expected.h` from a local copy of the original supplement:

```sh
python generate.py /path/to/desai_model_jgrb51665-sup-0002-ds01.txt
```

Native tests subtract the separately retained conventional libration component
from the composite calculators and check ocean values and rates against these
fixtures. Value tolerances are 2e-10 µas and 2e-11 µs. The rate tolerances are
2e-10 µas/s and 2e-11 µs/s, allowing double phase-reduction noise and the
production five-point harmonic stencil. A controlled daily-EOP factory test
also checks independent SOFA CIO matrix assembly, native UT1 rate, full rotation
derivative closure, unchanged default selection, and ambiguous custom-file or
unit combinations. Published LOD amplitudes are preserved as raw metadata;
the clock derivative comes from the selected UT1 function.
