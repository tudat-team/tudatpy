The exact finite-endpoint, constant-velocity first-post-Minkowskian monopole
calculation changes the selected LLR leg lengths by about 1.7–2.3 micrometres.
It does not explain the centimetre-scale temporal signature or the isolated
5.43 m discrepancy. No observation is solved again and no correction or default
is changed.

| Frozen observation | Uplink change (µm) | Downlink change (µm) | Half-sum at fixed endpoints (nm) |
| --- | ---: | ---: | ---: |
| APOLLO A15, 2017-03-09 01:19:15 UTC | +2.280324 | −2.280473 | −0.0745 |
| APOLLO A14, 2017-03-09 02:10:30 UTC | +2.281810 | −2.281955 | −0.0722 |
| Grasse A15, 2017-03-09 18:34:15.297622200200 UTC | +1.709192 | −1.709468 | −0.1379 |

Changes are sums of the Sun, Earth and Moon monopoles. The half-sums are
algebraic combinations at fixed endpoints, not recomputed observables. Their
sub-nanometre size is below the 3.23 nm reconstruction comparison with saved
native Shapiro totals; it is not an observable accuracy claim. The individual
leg effects are also below the existing micrometre-scale numerical limits of
the complete range model. A full residual rerun is not warranted by this check.

The primary expression is [Hees, Bertone and Le Poncin-Lafitte 2014](https://arxiv.org/pdf/1406.6600),
equations 26, 30, 34, 44–45. It uses a Lorentz factor γL, distinct from the PPN
parameter γ = 1 used here. Let β = v/c, γL = (1−β²)^−1/2,
R = |xB−xA| and N = (xB−xA)/R. For a body tangent worldline with reference
position x0 at t0, define

    UX = (xX−x0) + γL²/(1+γL) β[β·(xX−x0)] − γL v(tB−t0)
    A′ = UA + γL β R,  B′ = UB,  L′ = γL(1−β·N)R
    Δmoving = 2GM/c² γL(1−β·N) log[(|A′|+|B′|+L′)/(|A′|+|B′|−L′)]

The body's velocity is constant in this diagnostic. The reference epoch is the
native frozen-body epoch: Earth at its surface endpoint, Moon at reflection,
Sun at the leg midpoint. Δ is a length; delay in seconds is Δ/c. The native
comparison freezes the same body's position at that epoch and uses the static
logarithm, without a velocity prefactor or transformed endpoints.

The finite-endpoint change is not the heuristic |v|/c times the static delay.
For example, the solar effect reaches 0.622 µm per leg whereas that scale is
0.315 µm. The Earth reaches 2.858 µm and Moon 0.047 µm. The CSV separately
records the signed prefactor-only contribution, the transformed-endpoint
remainder and the unsigned heuristic scale. A velocity prefactor alone would
omit part of the finite-path result.

The first-PM calculation uses the straight null reference segment from equation
22: spatial points xB−λ(xB−xA), coordinate time tB−λR/c. The actual stored
event interval, which includes other delays and double epoch rounding, is not
inserted as a new null constraint. Its difference from R/c is recorded
(12–97 ns). This is the GR monopole approximation, exact in constant velocity;
acceleration, spin, multipoles and higher orders in G are outside this control.
Atmospheric and glass delays remain part of the frozen solved events, not new
terms in the vacuum null-segment formula.

The endpoint reconstruction uses the existing native body state, rotation and
full deformed station position at each saved event double. It is the same route
used by the earlier fixed-event Shapiro audit. For the four APOLLO legs, the
static calculation agrees with saved native correction components within
3.23 nm. Grasse has saved solved epochs but no separately saved static component;
its reconstruction uses the original full 14-NP/14-weather session environment.
It is not described as a direct extraction of native internal endpoint states.

All input binary doubles are retained through round-trip JSON serialization
and then evaluated with 80 decimal digits. Equation 23's independent potential
integral agrees with equation 45 within 6.2×10^−78 m; 50 versus 80 digit results
agree within 6.2×10^−48 m. Setting velocity to zero reproduces the static formula
exactly at this precision. These arithmetic checks do not improve the native
endpoint reconstruction or establish an error bound for omitted physics.

`endpoint_inputs.json` preserves all six legs, event epochs, BCRS endpoints,
body positions/velocities, evaluation epochs, c and GMs. The 18 separate body
contributions and precision controls are in `body_delay_comparison.csv`.
`summary.json` pins source/kernel/event/paper hashes, reconstruction limitations
and fixed-endpoint sums. `reproduce.py` repeats the saved-input arithmetic with
mpmath, without Tudat or new observation solves:

    python analysis/laser_ranging/results/moving_monopole_control/reproduce.py
