These are numeric provider records, not external model source code.

`manifest.json` records the complete source SHA256, exact one-based inclusive
source ranges, extract SHA256 and edits (none). The coefficient ranges refer to
the decompressed IERS source; its source SHA256 refers to the pinned gzip bytes.
The positive table retains all degree/order 1..180 rows. Degree-one coefficients
remain raw published data; the body-centered model selects degrees >=2.

The IERS reference is the separate **366-row gravity reference for 2004**.
Its header fixes the physical constants, gamma2 and k'_2. Tests use its supplied
m1/m2 columns (arcseconds converted to radians), not a modern mean-pole model.
The coefficient comparison budget is half the last printed digit of each
seven-significant-digit wobble/output propagated with the published linear
response; 1e-20 covers the finer printed constants. No fitted tolerances.

Four COST-G force rows span the day 2008-07-03. Their exact ITRF positions,
ICRF-to-ITRF rotation matrices, corrected xp/yp and GPS MJD epochs are retained.
The benchmark secular pole is xp_bar=0.055+0.001677*t,
yp_bar=0.3205+0.003460*t arcseconds, t=(MJD-51544.5)/365.25.
The exact supplied Gegout97 load Love numbers are indexed by degree, starting
with k'_0 on line 1. IERS forcing constants are used; the baseline EIGEN-6C4
header gives output GM=3.986004415e14 m^3/s^2, R=6378136.460 m.
Force comparisons use the existing Tudat normalized harmonic force utility and
rotate its output with the transpose of the supplied matrix. The fixed budget
is 1e-11 m/s^2. This is a matched historical-forcing check; production time/pole
routing is a separate factory responsibility.

Legacy audit: the local old ILR oceanPoleTideGravityFieldVariations implements
the same Eq. 6.23 A/B linear combination and precomputed degree response. That
algorithm structure is retained. Its uninitialized-member base dimensions,
max-index/size confusion and shifted k'_n denominator/application indices are
corrected; the modern backend also requires selected coverage, omits degree one,
and applies explicit source/output GM/R conversion. No external GPL program
implementation is used.

Primary equations: IERS Conventions (2010), updated Chapter 6, section 6.5,
https://iers-conventions.obspm.fr/content/chapter6/icc6.pdf .
COST-G conventions: https://ftp.tugraz.at/pub/ITSG/COST-G/softwareComparison/00README_simulation.txt .

Offline extraction from the pinned cache, without downloads:

```
python regenerate_provider_extracts.py /path/to/references/gravity
```

The script verifies complete source hashes before extracting, then verifies the
extract hashes before writing. It never changes provider records.
