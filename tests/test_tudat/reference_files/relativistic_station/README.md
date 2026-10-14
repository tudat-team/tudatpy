# IERS spatial station equation fixtures

Primary source: IERS Conventions (2010), Chapter 11, §11.2, equations 11.18 and 11.19:
https://iers-conventions.obspm.fr/content/chapter11/tn36_c11.pdf
PDF SHA256: `ca147e67fd88e924b5603cea217541a53612493fa12f65bd485409ba2879ea2a`.
The PDF is cached outside the source tree at
`earth-station-model-work/references/relativistic-stations/iers_chapter11.pdf`.

Equation 11.18 maps the unscaled geocentric spatial offset to the barycentric
spatial offset, omitting terms not present in this conventional approximation:

```
r_B = (1 - U/c²) r - v (v·r)/(2c²)
```

Equation 11.19 maps a TT-compatible offset to a TDB-compatible offset:

```
r_B = (1 - U/c² - L_C) r - v (v·r)/(2c²)
```

The correction fixtures in `unitTestRelativisticStationPosition.cpp` were
computed independently with Python's Decimal arithmetic, precision 70 decimal
digits, using exact decimal input strings and scalar component operations:

```
r = (4075578.385, 931852.890, 4801570.154) m
v = (27000, -14000, 6200) m/s
x_host = (0,0,0) m at the fixture epoch
x_sun = (149597870700,0,0) m; GM_sun = 1.32712440018e20 m³/s²
x_moon = (0,384400000,0) m; GM_moon = 4.9028e12 m³/s²
v_sun = v_moon = (0,0,0) m/s
c = 299792458 m/s
L_C = 0 or 1.48082686741e-8
U = GM_sun / 149597870700 + GM_moon / 384400000
Udot = GM_sun * 27000 / 149597870700² - GM_moon * 14000 / 384400000²
delta_i = -(U/c² + L_C)*r_i - v_i*sum_j(v_j*r_j)/(2c²)
delta_dot_i = -Udot*r_i/c²
```

The host has constant velocity, the station offset is stationary in GCRS, and
orientation is the identity, so no acceleration or station-motion terms enter
these fixture derivatives. All GM values are constant.

The separate exact Lorentz fixture uses a stationary comoving point and constant
host velocity. At equal barycentric time the component along the velocity is
contracted by `1/gamma = sqrt(1-v²/c²)`. Its exact correction is evaluated with
`-beta²/(1+sqrt(1-beta²))` to avoid subtraction cancellation. The tolerance budgets
the deliberately omitted fourth-order term. A moving GCRS station would also
require its position-dependent event time; this test makes no such comparison.

Finite differences of moving/displaced stations validate the derivative of the
declared geocentre-epoch approximation. They do not establish accuracy of a
complete IAU event transformation or barycentric-time station worldline. These
fixtures are original equation evaluations, not copied implementation code or
cross-software station trajectories.

Reproduce and check the frozen literals without Tudat imports:

```sh
python3 tests/test_tudat/reference_files/relativistic_station/generate_fixtures.py --check
```
