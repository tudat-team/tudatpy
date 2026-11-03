# Pluto--Charon position-angle/separation validation data

`pm0001.txt` is the 60-observation Pluto--Charon HST series published by
Tholen and Buie (1997) and distributed by the IMCCE Natural Satellites Data
Base:

- https://nsdb.imcce.fr/obspos/OBS_COLL/P/pm0001.txt
- https://nsdb.imcce.fr/obspos/OBS_COLL/P/pm0001.html

The columns are UTC Julian date at exposure midpoint, separation (arcsec),
position angle (degrees, J2000 north through east), and the two residuals
published with the original fitted orbit.  The UTC interpretation is confirmed
by matching all 60 epochs to their MAST HST exposure midpoints to within 6 ms.

`plu060_pm0001_subset.bsp` is an offline test subset derived from JPL's PLU060
satellite ephemeris:

- https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/plu060.bsp

It contains only NAIF bodies 399 (Earth), 999 (Pluto), and 901 (Charon), each
relative to the solar-system barycentre in J2000.  Type-9 segments were sampled
from PLU060 every three hours from 1992-05-20 through 1993-08-20 TDB, using
degree-15 interpolation.  At all 60 reception epochs and both independently
iterated transmission epochs, its maximum position differences from PLU060
were 0.000082 m for Earth, 0.002571 m for Pluto, and 0.002388 m for Charon.

The unit test deliberately uses Earth's centre as the receiver.  Archived HST
J2000 state vectors from the corresponding 60 MAST standard-header packets
show that this approximation changes the predicted separation by at most
0.0017 mas and the transverse position-angle component by at most 0.0011 mas,
well below both the data precision and the residual threshold.  This avoids
embedding a second, independently interpolated spacecraft ephemeris in the
test while keeping the approximation explicit.

## Mars-satellite true-of-date validation

`mm0012.txt` contains 166 Mars-satellite CCD observations made with the
Jacobus Kapteyn Telescope on La Palma between 1988-09-21 and 1988-09-30.  Jones,
Sinclair, and Williams (1989) published the series, which is distributed by
the IMCCE Natural Satellites Data Base:

- https://nsdb.imcce.fr/obspos/OBS_COLL/M/mm0012.txt
- https://nsdb.imcce.fr/obspos/OBS_COLL/M/mm0012.html

NSDB identifies the values as topocentric, astrometric separation and position
angle referred to the true equator and equinox of date.  The observations are
corrected for satellite phase; their differential-refraction status is not
documented.  The unit test therefore uses the series specifically to validate
the selectable true-of-date position-angle pole, not atmospheric refraction.

`mar099_de441_mm0012_subset.bsp` is an offline Type-9 SPK constructed from JPL
Horizons geometric ICRF state vectors for Earth (399), Mars (499), Phobos
(401), and Deimos (402), relative to the solar-system barycentre.  It covers
1988-09-19 through 1988-10-02 TDB, sampled every five minutes and interpolated
with degree-15 polynomials.  Horizons reports DE441 as the source for Earth
and MAR099 as the source for Mars, Phobos, and Deimos.  Independent Horizons
states at 288 half-grid epochs per body differ from the subset by at most
1.19 m in position.  The subset SHA-256 is
`a504982dfa91452a20e0babdf2f81245f0a8968d9715384a1b7115acc3d85c07`.

The JKT geodetic coordinates used by the test are 28 deg 45 min 40.1 sec N,
17 deg 52 min 41.2 sec W, and 2364 m, from the Isaac Newton Group telescope
page: https://www.ing.iac.es/astronomy/telescopes/jkt/

## Saturn-satellite apparent-direction validation

`qiao1999_table4.txt` contains the 41 apparent topocentric measurements in
Table 4 of Qiao et al. (1999), *1994-1996 CCD astrometric observations of
Saturn's satellites and comparison with theories*, Astronomy and Astrophysics
Supplement Series 137, 1-5, DOI 10.1051/aas:1999238:

- https://aas.aanda.org/articles/aas/pdf/1999/10/ds7914.pdf

Columns are UTC Julian date at exposure midpoint, NAIF ID of the measured
satellite, position angle (degrees), and separation (arcseconds). Titan (NAIF
606) is the reference object in every row. The paper identifies these as
apparent topocentric observations, states that differential refraction was
removed, and states that stellar aberration and parallax were not removed. The
observing site is Sheshan (E 121 deg 11 min 03.3 sec, N 31 deg 05 min 46.1 sec,
altitude 97 m).

`sat441_de441_qiao1999_subset.bsp` is a compact test-only Type-9 SPK generated
from JPL Horizons J2000 state vectors relative to the solar-system barycentre
from 1996-10-18 through 1996-10-21 at five-minute intervals. Earth uses DE441;
Enceladus, Tethys, Dione, Rhea, and Titan use SAT441L. It contains NAIF bodies
399, 602, 603, 604, 605, and 606 with degree-15 interpolation. Fresh Horizons
queries at off-grid epochs on the observation night agree at the sub-metre
level in position. The subset is for this regression test only, not as a
general-purpose ephemeris.
Its SHA-256 digest is
`bac178e2dc038637e286a4c0f3f548dd34cf8992fee73352c604a24ffbc6c102`.

## Nereid B1950-reference-pole validation

`nm1007.txt` contains 17 photographic positions of Nereid relative to Neptune
from Mauna Kea, McDonald Observatory, and La Silla between 1977 and 1984:

- https://nsdb.imcce.fr/obspos/OBS_COLL/N/nm1007.txt
- https://nsdb.imcce.fr/obspos/OBS_COLL/N/nm1007.html

NSDB unambiguously identifies these observations as topocentric, astrometric,
B1950, and UTC.  The atmospheric-reduction status is not documented.  No
atmospheric-refraction model is used: this case tests the much larger B1950
reference-pole signal in position angle, while separation remains independent
of the pole.  It is retained as a secondary converted-coordinate case because
a screen of the 140 NSDB collections whose metadata mention separation,
position angle, or distance found no native B1950 P/S series in which the
frame signal is similarly large relative to the measurement scatter.  The
native `um0027` P/S case below remains the primary B1950 residual validation.
Published modern re-reductions of the two Nereid observing series give
coordinate RMS values of approximately 0.17--0.22 arcsec.

`nep105_nm1007_subset.bsp` is sampled from JPL/NAIF's NEP105 kernel and
contains Earth, Neptune, and Nereid relative to the solar-system barycentre
from 1977 through 1984.  It is a degree-15 Type-9 SPK sampled every six hours.
Half-grid checks against the source kernel have maximum position errors of
0.003 m for Earth, 11 m for Neptune, and 0.043 m for Nereid.  Its SHA-256 is
`771bd4d037958e97266c525c8d065cbf8d780042e7afd638a7b986b6ae6c419a`;
the detailed check is in `nep105_nm1007_subset_provenance.txt`.  No spacecraft
tracking or Voyager observation data are included or used.

## Native Uranian-satellite B1950 P/S validation

`um0027.txt` contains 156 photographic observations of Ariel, Umbriel,
Titania, and Oberon relative to Uranus, made at Tokyo-Okayama from 1964 to
1977 by Tomita and Soma (1979):

- https://nsdb.imcce.fr/obspos/OBS_COLL/U/um0027.txt
- https://nsdb.imcce.fr/obspos/OBS_COLL/U/um0027.html

These are native separation (arcseconds) and position-angle (degrees)
observations; no differential-coordinate conversion is involved.  NSDB
explicitly identifies them as topocentric, astrometric, B1950, and ET.  The
test therefore treats the recorded calendar fields as the uniform ephemeris
time scale (TDB is the implemented modern realization; TT/TDB differences are
negligible at this precision).  A separate screening calculation also tried
UTC and shifts of plus or minus half the listed exposure time.  None changes
the approximately 0.8-arcsec residual scale enough to overturn the explicit
ET metadata.  Atmospheric refraction is not modeled.

`ura184_de441_um0027_subset.bsp` is a compact degree-15 Type-9 SPK generated
from five-minute Horizons geometric J2000 vectors on the 11 observing runs.
Horizons identifies Earth as DE441 and Uranus, Ariel, Umbriel, Titania, and
Oberon as `ura184_merged`.  Its SHA-256 is
`11e66b36b31273995300756d4ed7256eb2285f509286c96d4cc05f9971a1e75f`;
all source queries are in `ura184_de441_um0027_subset_provenance.txt`.  Fresh
hourly Horizons queries offset by two minutes from the interpolation nodes
agree within 0.113 km in position for all six bodies.
