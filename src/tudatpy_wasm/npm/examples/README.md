# Tudat WASM Examples

This directory contains JavaScript examples that demonstrate how to use the Tudat WASM library. These examples cover all major functionality available in the WASM module.

## Node.js Examples

Run these examples with Node.js after installing the package:

```bash
npm install @tudat/tudatpy-wasm
node examples/01_keplerian_satellite_orbit.js
```

### Available Examples

#### Propagation Examples

| Example | Description | Concepts |
|---------|-------------|----------|
| `01_keplerian_satellite_orbit.js` | Basic two-body orbit propagation | Kepler problem, orbital period, vis-viva |
| `02_perturbed_satellite_orbit.js` | Orbit with J2 perturbations | Secular drift, nodal regression, apsidal precession |
| `04_thrust_satellite.js` | Continuous thrust with mass propagation | Ion propulsion, Tsiolkovsky equation |
| `05_solar_system_propagation.js` | Multi-body planetary propagation | N-body problem, planetary motion |

#### Element Conversion and Coordinates

| Example | Description | Concepts |
|---------|-------------|----------|
| `06_unit_conversions.js` | Physical constants and units | SI units, astronomical constants |
| `07_orbital_elements.js` | Orbital element conversions | Keplerian, Cartesian, anomalies, MEE |
| `18_coordinate_systems.js` | Coordinate transformations | Spherical, geodetic, RSW, ENU, PQW |

#### Gravity and Perturbations

| Example | Description | Concepts |
|---------|-------------|----------|
| `08_gravity_models.js` | Spherical harmonics gravity | J2-J6 effects, Legendre polynomials, geoid |
| `09_aerodynamics.js` | Atmospheric models and drag | Exponential, NRLMSISE-00, drag coefficient |
| `10_radiation_pressure.js` | Solar radiation pressure | Cannonball model, solar sails, SRP effects |

#### Mission Design

| Example | Description | Concepts |
|---------|-------------|----------|
| `03_lambert_targeting.js` | Lambert problem solver | Interplanetary transfers, TOF, porkchop plots |
| `11_gravity_assist.js` | Flyby maneuvers | Turn angle, B-plane, V∞ leveraging, Oberth effect |
| `17_escape_capture.js` | Escape and capture trajectories | SOI, hyperbolic excess, C3, capture burns |
| `19_interplanetary_transfer.js` | Interplanetary mission design | Hohmann transfers, synodic periods, phase angles |

#### Three-Body Dynamics

| Example | Description | Concepts |
|---------|-------------|----------|
| `13_cr3bp.js` | Circular restricted 3-body problem | Lagrange points, Jacobi constant, halo orbits |

#### Numerical Methods

| Example | Description | Concepts |
|---------|-------------|----------|
| `12_numerical_integrators.js` | Integration algorithms | RK4, RKF7(8), DOPRI8(7), adaptive step size |

#### Time and Ephemerides

| Example | Description | Concepts |
|---------|-------------|----------|
| `14_ephemerides.js` | Keplerian and tabulated ephemeris | Kepler propagation, interpolation, planetary positions |
| `15_time_conversions.js` | Time scales and Earth orientation | UTC, TAI, TT, Julian date, sidereal time, EOP |
| `16_tle_propagation.js` | Two-Line Element propagation | TLE format, SGP4, satellite tracking |

#### Orbit Determination

| Example | Description | Concepts |
|---------|-------------|----------|
| `20_orbit_determination.js` | State estimation fundamentals | Least squares, Kalman filter, covariance, observability |

#### Advanced Propagation

| Example | Description | Concepts |
|---------|-------------|----------|
| `21_reentry_trajectory.js` | Atmospheric reentry simulation | Drag, heating, dynamic pressure, entry corridor |
| `24_rotational_dynamics.js` | Satellite attitude dynamics | Euler equations, quaternions, torque-free motion |
| `25_multi_satellite.js` | Multi-satellite propagation | Relative motion, CW equations, differential drag |
| `27_earth_moon_transfer.js` | Earth-Moon transfer with thrust | Low-thrust, mass depletion, multi-body |
| `28_two_stage_rocket.js` | Two-stage rocket ascent on Mars | Staging, dual-thrust, atmospheric drag |
| `29_juice_flybys.js` | JUICE mission flybys | Jovian system, moon flybys, B-plane |
| `31_coupled_dynamics.js` | Coupled translational-rotational dynamics | Gravity gradient torque, attitude coupling |
| `32_impact_manifolds.js` | Libration point orbit manifolds | CR3BP, stable/unstable manifolds, halo orbits |
| `33_separation_diff_drag.js` | Satellite separation via differential drag | Formation flying, drag modulation |

#### Estimation and Analysis

| Example | Description | Concepts |
|---------|-------------|----------|
| `22_covariance_propagation.js` | Uncertainty propagation | State covariance, state transition matrix, error growth |
| `23_linear_sensitivity.js` | Sensitivity analysis | Parameter derivatives, J2 drift sensitivity, ground track |

#### Advanced Mission Design

| Example | Description | Concepts |
|---------|-------------|----------|
| `26_mga_optimization.js` | Multiple gravity assist optimization | Cassini-like trajectory, EVVEJSA, global optimization |
| `30_hodographic_shaping.js` | Low-thrust trajectory shaping | Hodographic method, basis functions, electric propulsion |
| `34_cassini_mga.js` | Cassini-1 MGA trajectory design | VVEJSA sequence, Lambert targeting, flyby altitude |
| `35_earth_mars_window.js` | Earth-Mars transfer window analysis | Porkchop plots, C3, arrival V-infinity |
| `36_low_thrust_window.js` | Low-thrust Earth-Mars transfers | Edelbaum approximation, electric propulsion |
| `41_mga_trajectories.js` | MGA trajectory types | Unpowered legs, DSMs, low-thrust overview |

#### Estimation Examples

| Example | Description | Concepts |
|---------|-------------|----------|
| `37_full_estimation.js` | Full state estimation | Batch least squares, range/range-rate observations |
| `38_tle_estimation.js` | TLE-based state estimation | TLE parsing, mean elements, SGP4 propagation |
| `39_galilean_moons.js` | Galilean moons state estimation | Multi-body estimation, Laplace resonance |
| `42_covariance_analysis.js` | Covariance estimation analysis | Correlation coefficients, uncertainty ellipsoids |
| `43_estimation_dynamics.js` | Different dynamical models | Truth vs estimation model, model mismatch |
| `45_mpc_observations.js` | MPC observation handling | 80-column format, observatory codes, astrometry |
| `46_tle_ephemeris.js` | TLE-based ephemeris | TLE format, SGP4 propagation, mean elements |

#### Optimization Examples

| Example | Description | Concepts |
|---------|-------------|----------|
| `40_asteroid_optimization.js` | Asteroid orbit optimization | Evolutionary algorithms, multi-objective optimization |
| `44_himmelblau_optimization.js` | Himmelblau function optimization | GA, DE, multimodal optimization |

## Browser Demo

Open `browser_demo.html` in a web browser to interact with Tudat WASM through a graphical interface. Features include:

- **Orbital Element Converter**: Convert between Keplerian and Cartesian elements
- **Orbit Visualization**: See orbits rendered on a canvas
- **Lambert Solver**: Compute interplanetary transfer trajectories

### Running the Browser Demo

1. Build the WASM module:
   ```bash
   cmake --build build-wasm --target tudatpy_wasm
   ```

2. Serve the files with a local web server:
   ```bash
   cd src/tudatpy_wasm/npm/examples
   python -m http.server 8080
   ```

3. Open http://localhost:8080/browser_demo.html in your browser

## Example Structure

Each example follows this pattern:

```javascript
const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    // Initialize the WASM module
    const tudat = await createTudatModule();

    // Create vectors and matrices
    const state = new tudat.Vector6d();
    state.set(0, 7000000);  // x position [m]
    // ... set other elements

    // Call Tudat functions
    const result = tudat.astro.element_conversion.keplerian_to_cartesian(
        state,
        3.986004418e14  // Earth GM
    );

    // Use results
    console.log('Position:', result.get(0), result.get(1), result.get(2));

    // IMPORTANT: Clean up WASM objects
    state.delete();
    result.delete();
}

main();
```

## Memory Management

WebAssembly objects must be explicitly deleted when no longer needed:

```javascript
const vector = new tudat.Vector6d();
try {
    // Use the vector
    const result = someFunction(vector);
    // ... process result
    result.delete();  // Clean up result
} finally {
    vector.delete();  // Always clean up input
}
```

## API Reference

The JavaScript API mirrors the Python tudatpy structure:

### Constants
```javascript
tudat.constants.GRAVITATIONAL_CONSTANT
tudat.constants.SPEED_OF_LIGHT
tudat.constants.ASTRONOMICAL_UNIT
```

### Astrodynamics
```javascript
// Element conversion
tudat.astro.element_conversion.keplerian_to_cartesian(kepler, GM)
tudat.astro.element_conversion.cartesian_to_keplerian(cartesian, GM)

// Two-body dynamics
tudat.astro.two_body_dynamics.compute_kepler_orbit_period(a, GM)
tudat.astro.two_body_dynamics.propagate_kepler_orbit(state, dt, GM)

// Time representation
const dt = new tudat.astro.time_representation.DateTime(2024, 1, 15, 12, 0, 0);
const epoch = dt.epoch();
```

### Trajectory Design
```javascript
tudat.trajectory_design.transfer_trajectory.solve_lambert_problem(
    pos1, pos2, tof, GM, isRetrograde
)
```

## TypeScript Support

TypeScript definitions are included. Use with full type safety:

```typescript
import createTudatModule, { TudatModule, Vector6d } from '@tudat/tudatpy-wasm';

async function main(): Promise<void> {
    const tudat: TudatModule = await createTudatModule();
    const state: Vector6d = new tudat.Vector6d();
    // ... fully typed API
}
```

## Topics Covered

These examples comprehensively cover:

- **Orbital Mechanics**: Two-body problem, Keplerian motion, orbital elements
- **Perturbations**: J2, atmospheric drag, solar radiation pressure, third-body
- **Coordinate Systems**: Cartesian, spherical, geodetic, orbital frames
- **Time Systems**: UTC, TAI, TT, Julian dates, sidereal time
- **Numerical Methods**: Runge-Kutta, adaptive integration, state propagation
- **Mission Design**: Hohmann transfers, Lambert problem, gravity assists
- **Three-Body Dynamics**: CR3BP, Lagrange points, libration point orbits
- **Orbit Determination**: State estimation, Kalman filtering, covariance
- **Ephemerides**: Keplerian propagation, tabulated data, TLE/SGP4

## More Information

- [Tudat Documentation](https://docs.tudat.space/)
- [Python tudatpy Examples](https://github.com/tudat-team/tudatpy/tree/master/examples)
- [API Documentation](https://tudat-team.github.io/tudatpy-wasm/)
