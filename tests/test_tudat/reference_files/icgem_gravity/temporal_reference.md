Independent temporal ICGEM references
====================================

The standalone temporal-model tests use the existing verbatim
`eigen6s4_v2_c20_1950_2006.gfc` extract. Its original source, full-file SHA256,
line ranges, fixture SHA256 and CC-BY-4.0 attribution are pinned in `sources.json`.
Cite Foerste et al. (2016), EIGEN-6S4 Version 2, GFZ Data Services,
https://doi.org/10.5880/icgem.2016.008.

No source clock or year definition is inferred. These tests explicitly declare
EIGEN printed dates as **TAI**, and years as elapsed TAI SI seconds / 31557600.
This is a controlled consuming convention, not a claim about the provider's
otherwise unstated clock. The invalid literal-HHMM token `20041226.0060` receives
an exact-token caller override to **2004-12-26 01:00:00 TAI**; that value is a test
choice, not an inferred rollover or fractional-day interpretation. Omitting the
override must reject. Comparisons at 2003-01-01 and the override boundary check
the exact original-double TDB image and its immediately preceding image.
The boundary belongs to the right interval, with no epsilon snapping. A query
before the first baseline or at the final half-open end rejects.

The 2000 and 2001 C20 baselines each retain the annual/semiannual terms' own
1950 reference origin and validity through 2003. Their phase is not restarted
at each annual baseline. At the 2003 shared boundary both the selected baseline
and periodic interval change. The fixed consuming nominal C20 is `-4.84165e-4`;
it is subtracted once from the complete published temporal coefficient.

`controlled_temporal_c21.gfc` is an explicitly synthetic analytic input, not a
provider product. It has an unnormalized C21/S21 baseline, independent trend
origin at `20000101.1234`, and annual/semiannual cosine and sine terms with origin
`20000101.0000`. Its source dates are explicitly **TT**. The source/output GM
ratio is 4/5, radius ratio 8/7, and conversion of unnormalized degree-two/order-one
coefficients contributes `sqrt(3/5)`, so the complete conversion factor is
`sqrt(3/5)*(4/5)*(8/7)^2`. The consuming nominal C21/S21 is `1e-6/-2e-6`.
Static-only selected coefficients produce zero anomaly. Both year conventions
are evaluated in leap year 2000 and the following year: elapsed TAI Julian years,
and fractional Gregorian source-calendar years. A separate compact format-1.0
case has date-only midnight origins inherited by its undated trend/periodic
records; its quarter-Julian-year value follows directly from sin(pi/2)=1 and
cos(pi/2)=0.

Reproduce the scalar references offline with Python's standard library:

```
python tests/test_tudat/reference_files/icgem_gravity/freeze_temporal_reference.py
```

The script checks both input hashes, selects the cited EIGEN source rows
explicitly, and evaluates scalar equations with 80-digit Decimal Taylor
trigonometry. It imports no Tudat module or numerical Python package. It does
not implement the runtime model's interval search or clock conversion. Frozen
results include computational budgets, absolute coefficient rates and, separately,
half-final-printed-unit bounds for the provider amplitudes. The latter are not
used as test tolerances: tests verify the formula applied to the supplied printed
numbers more tightly than the physical source precision.

Computational budgets use conservative binary64 epsilon `2^-52` and
`gamma(k)=k*epsilon/(1-k*epsilon)`. They include `gamma(32)` times the absolute
coefficient/nominal magnitudes for source parsing, normalization, GM/radius
scaling, additions and baseline subtraction. Periodic terms add their absolute
amplitude times `[gamma(12)*abs(phase)+4*epsilon]` for phase/trigonometric rounding.
The controlled calendar case additionally bounds subtraction of calendar-year
coordinates near year 2002 by `epsilon*2002*conversion_factor*1e-6`.
At runtime, each comparison adds the recorded absolute derivative times four
input-TDB ULPs to cover scalar conversion rounding and the requested neighboring
image. This changes a coefficient comparison budget only; it never changes
interval classification or the query epoch. No tolerance was fitted by running
the C++ implementation.
