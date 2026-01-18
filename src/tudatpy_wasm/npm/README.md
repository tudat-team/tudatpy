# @tudat/tudatpy-wasm

Tudat astrodynamics library compiled to WebAssembly for browser and Node.js environments.

This package provides the same API as the Python [tudatpy](https://docs.tudat.space/) library, enabling orbital mechanics and astrodynamics simulations directly in JavaScript/TypeScript applications.

## Installation

```bash
npm install @tudat/tudatpy-wasm
```

## Quick Start

### Basic Usage (Node.js)

```javascript
const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    const tudat = await createTudatModule();

    // Convert Keplerian elements to Cartesian state
    const kepler = new tudat.Vector6d();
    kepler.set(0, 7000000);    // semi-major axis [m]
    kepler.set(1, 0.01);       // eccentricity
    kepler.set(2, 0.5);        // inclination [rad]
    kepler.set(3, 0);          // argument of periapsis [rad]
    kepler.set(4, 0);          // RAAN [rad]
    kepler.set(5, 0);          // true anomaly [rad]

    const GM = 3.986004418e14;  // Earth gravitational parameter
    const cartesian = tudat.astro.element_conversion.keplerian_to_cartesian(kepler, GM);

    console.log('Position [m]:', cartesian.get(0), cartesian.get(1), cartesian.get(2));
    console.log('Velocity [m/s]:', cartesian.get(3), cartesian.get(4), cartesian.get(5));

    // Clean up (important for memory management)
    kepler.delete();
    cartesian.delete();
}

main();
```

### ES Module Usage

```javascript
import createTudatModule from '@tudat/tudatpy-wasm';

const tudat = await createTudatModule();
// ... use tudat module
```

### TypeScript Usage

```typescript
import createTudatModule, { TudatModule, Vector6d } from '@tudat/tudatpy-wasm';

async function computeOrbit(): Promise<void> {
    const tudat: TudatModule = await createTudatModule();

    const state: Vector6d = new tudat.Vector6d();
    // ... type-safe usage
}
```

### Browser Usage

```html
<script type="module">
import createTudatModule from './node_modules/@tudat/tudatpy-wasm/dist/tudatpy_wasm.mjs';

async function init() {
    const tudat = await createTudatModule({
        // Optional: customize WASM file location
        locateFile: (path) => `/wasm/${path}`
    });

    // Use tudat module
    const period = tudat.astro.two_body_dynamics.compute_kepler_orbit_period(
        7000000,  // semi-major axis [m]
        3.986004418e14  // Earth GM [m³/s²]
    );

    console.log('Orbital period:', period, 'seconds');
}

init();
</script>
```

## API Overview

The API mirrors the Python tudatpy structure:

### Constants

```javascript
tudat.constants.GRAVITATIONAL_CONSTANT
tudat.constants.SPEED_OF_LIGHT
tudat.constants.ASTRONOMICAL_UNIT
tudat.constants.EARTH_GRAVITATIONAL_PARAMETER
// ... more physical constants
```

### Astrodynamics (`tudat.astro`)

```javascript
// Element conversion
tudat.astro.element_conversion.keplerian_to_cartesian(kepler, GM)
tudat.astro.element_conversion.cartesian_to_keplerian(cartesian, GM)
tudat.astro.element_conversion.mean_to_eccentric_anomaly(e, M)

// Two-body dynamics
tudat.astro.two_body_dynamics.compute_kepler_orbit_period(a, GM)
tudat.astro.two_body_dynamics.propagate_kepler_orbit(state, dt, GM)

// Time representation
const dt = new tudat.astro.time_representation.DateTime(2024, 1, 15, 12, 0, 0);
const epoch = dt.epoch();  // seconds since J2000
```

### Dynamics (`tudat.dynamics`)

```javascript
// Environment setup
const bodySettings = tudat.dynamics.environment_setup.get_default_body_settings(
    ['Earth', 'Moon', 'Sun'],
    'Earth',
    'J2000'
);
const bodies = tudat.dynamics.environment_setup.create_system_of_bodies(bodySettings);

// Propagation setup
const integratorSettings = tudat.dynamics.propagation_setup.integrator.runge_kutta_fixed_step_size(60.0);
const terminationSettings = tudat.dynamics.propagation_setup.propagator.time_termination(86400.0);

// Run simulation
const simulator = tudat.dynamics.simulator.create_dynamics_simulator(
    bodies, integratorSettings, propagatorSettings
);
```

### SPICE Interface (`tudat.interface`)

```javascript
tudat.interface.spice.load_kernel('/path/to/kernel.bsp');
const state = tudat.interface.spice.get_body_cartesian_state_at_epoch(
    'Earth', 'Sun', 'ECLIPJ2000', 'NONE', epoch
);
```

## Memory Management

WebAssembly objects must be explicitly deleted when no longer needed:

```javascript
const vector = new tudat.Vector6d();
// ... use vector
vector.delete();  // Free memory
```

For complex simulations, consider using try/finally:

```javascript
const state = new tudat.Vector6d();
try {
    // ... simulation code
} finally {
    state.delete();
}
```

## Vector and Matrix Types

| Type | Description | Size |
|------|-------------|------|
| `Vector3d` | 3D position/velocity vector | Fixed 3 |
| `Vector6d` | State vector (position + velocity) | Fixed 6 |
| `VectorXd` | Dynamic-size vector | Variable |
| `Matrix3d` | 3×3 rotation matrix | Fixed 3×3 |

```javascript
// Creating vectors
const v3 = new tudat.Vector3d();
const v6 = new tudat.Vector6d();
const vx = new tudat.VectorXd();
vx.resize(10);

// Accessing elements
v6.set(0, 7000000);
const x = v6.get(0);

// Converting to JavaScript arrays
const arr = v6.toArray();  // [x, y, z, vx, vy, vz]
```

## Browser Considerations

### WASM File Location

By default, the module looks for `tudatpy_wasm.wasm` in the same directory as the JavaScript file. Customize this with `locateFile`:

```javascript
const tudat = await createTudatModule({
    locateFile: (path, prefix) => {
        if (path.endsWith('.wasm')) {
            return '/assets/wasm/' + path;
        }
        return prefix + path;
    }
});
```

### Memory

The module starts with 256MB of memory and can grow as needed. For large simulations, consider monitoring memory usage:

```javascript
console.log('Memory pages:', tudat.HEAP8.length / (64 * 1024));
```

## Documentation

Full API documentation is available at:
- [TypeScript API Docs](https://tudat-team.github.io/tudatpy-wasm/)
- [Python tudatpy Docs](https://docs.tudat.space/) (API is equivalent)

## Requirements

- **Node.js**: 16.0.0 or later
- **Browser**: Modern browser with WebAssembly support (Chrome 57+, Firefox 52+, Safari 11+, Edge 16+)

## License

BSD-3-Clause

## Links

- [Tudat Documentation](https://docs.tudat.space/)
- [GitHub Repository](https://github.com/tudat-team/tudatpy)
- [Issue Tracker](https://github.com/tudat-team/tudatpy/issues)
