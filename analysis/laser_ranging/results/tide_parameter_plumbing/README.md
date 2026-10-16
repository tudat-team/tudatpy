The bounded audit finds no lunar radius, Love-number order, degree normalization
or duplicate-deformation error in the current instantiated LLR environment.
It evaluates native station positions at three saved March 9 reflection epochs,
without solving observations, fitting parameters or averaging an arc.

| Check | Maximum vector difference |
| --- | ---: |
| Native actual tide versus independent JPL equation and original constants | 0.145 nm |
| Native degree 2 plus degree 3 versus actual combined model | 0.261 nm |
| New probe station versus existing intrinsic reflector | exactly zero |
| Final reflector versus one tide plus one spatial mapping | 0.214 nm |
| Doubled reference radius, degree 2: expected factor 16 | 1.993 nm |
| Doubled reference radius, degree 3: expected factor 32 | 4.093 nm |
| Doubled shape radius with unchanged explicit deformation radius | exactly zero |
| Doubled shape radius with NaN deformation radius: expected factor 16 | 1.993 nm |

The larger errors in the scaled checks amplify the rounding incurred when
subtracting metre-scale motions from million-metre Cartesian positions. No
observation or solver tolerance is adjusted. Four reflectors already present in
the production environment are checked directly; a fifth native probe station
uses the same body deformation model at the published L2 Cartesian coordinates.

The original DE440 ASCII header independently gives AM = 1738 km and
RE = 6378.1366 km. GMS, GMB, EMRAT and AU yield the Earth, Moon and Sun GMs
used by the instantiated bodies, within floating representation roundoff.
The lunar displacement factory receives Earth/Sun, radius 1738000 m and
(h2,l2) = (0.0423,0.0107), (h3,l3) = (0.0234,0.0030). These are displacement
Love/Shida numbers multiplying unnormalized Legendre polynomials, without a
spherical-harmonic gravity coefficient normalization or potential Love number.
The independent calculation uses direct geocentric SPICE vectors, the original
header constants, and polynomial gradients; it does not read the native model's
cached tidal components.

At these epochs the complete degree-2 vector norms are about 228–509 mm and
degree-3 norms 0.47–1.07 mm across the five reflectors. This is consistent with
the physical scale in [JPL2020 section 3.12 and Appendix A2](https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/LLR_Model_2020_DR.pdf).
Appendix A2 tabulates analytic constant components, not instantaneous tides.
The full equilibrium tide, including its constant component, is retained here;
no finite-arc mean is used as an exact A2 reference or subtracted. The native
once-only reconstruction also excludes an extra permanent-tide subtraction.

The argument path is Python `basic_solid_body_tidal` →
`BasicSolidBodyDeformationSettings::bodyReferenceRadius` → the native factory →
`BasicTidalBodyDeformation`. An explicit finite radius is passed unchanged.
The actual NaN fallback uses the body's shape-model average radius. Python's
current docstring instead claims a gravity-field radius fallback; this is a
documentation mismatch, bypassed by the explicit radius in the LLR script.
Doubling the shape radius empirically confirms both paths. Moon reflector
coordinates are Cartesian and are not projected onto the chosen shape sphere.

Earth's IERS factory fixes its tide radius at 6378136.6 m, with principal
(h2,l2) = (0.6078,0.0847) and (h3,l3) = (0.292,0.015), plus the existing
IERS latitude/frequency terms. The analysis Earth shape is the distinct WGS84
ellipsoid, and the selected station coordinates are Cartesian. Factory GMs
come from the instantiated Earth/Moon/Sun gravity fields. Static assembly
contains one IERS model and one pole-tide model; each station has one
body-deformation motion entry. Intrinsic and final stations are separate
evaluations, not successive additions to an already deformed nominal position.
The native deformation-model list is not exposed to Python; the empirical
lunar scaling and position checks provide the runtime evidence instead of
claiming direct object-member inspection.

`native_displacement_controls.csv` preserves all 105 vector evaluations.
Its independent columns use the original nominal radius; the deliberately
scaled native variants are interpreted using the factors in the table above.
`summary.json` pins settings, body GMs, source/header/kernel hashes and errors.
The portable probe is `../../audit_llr_tide_plumbing.py`; for example::

    python analysis/laser_ranging/audit_llr_tide_plumbing.py \
      --data-dir "$LLR_DATA" --kernel-dir "$LLR_KERNELS" \
      --event-diagnostics "$FROZEN_EVENT_DIAGNOSTICS" --output "$NEW_OUTPUT"

The data directory also needs the unchanged [original DE440 header](https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/header.440),
SHA256 `0f4636b663e6f00af9efa2b605ef23635489f56b23dd1a5e3821d3acc6ac3359`.
The frozen event JSON is identified by its hash in the summary.

The epoch and apparatus audit supports using each published tag unchanged.
[APOLLO's format](https://tmurphy.physics.ucsd.edu/apollo/norm_pts.html) labels it
UTC launch time. [Its reduction paper, section 2](https://doi.org/10.1088/1538-3873/ACF787)
constructs a representative time, fits the residual distribution with the
reflector response, then restores prediction and trend at that time. Its
fiducial path is referenced to the telescope axis intersection. A new average
of raw RTTs over the NP duration, a half-duration tag shift, or another
fiducial-to-axis correction would not reproduce that definition.

[Grasse's instrumentation paper, section 2](https://doi.org/10.1029/2019EA000785)
likewise describes a representative launch time and RTT after consolidating
individual shots. [CRDv1 section 3.2](https://ilrs.gsfc.nasa.gov/docs/2009/crd_v1.01.pdf)
defines event 2 at the station system reference point and expects NP RTT to
include station system delay. The frozen records have that correction flag
set. Bin RMS and peak-minus-mean describe the return distribution; they do not
authorize another centroid correction without the source reduction procedure
and raw shots. Apparatus zero points and return-profile weighting remain
metrology limitations, not independently determined model corrections here.

The isolated +5.43 m NP has one additional archive cross-check:
the EDC daily CRDv1 product repeats exactly the frozen monthly NP record
(daily line 15 versus monthly line 171), including epoch and RTT. This rules
out a monthly-aggregation transcription difference. The V2 URL returned no
CRD product, so no corrected alternative is claimed. `daily_source_check.json`
records the URL, daily hash and exact record. The discrepancy remains included
and unexplained; no new evidence supports rejecting or changing it.
