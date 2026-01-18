# Tudat WASM Test Suite

This directory contains the WebAssembly test suite for Tudat, including a browser-based visualization demo with 3D orbital mechanics.

## Screenshot

The web UI features a modern space-themed design with real-time 3D orbital visualization and interactive charts:

![Tudat WASM Test Runner](docs/screenshot.png)

## Quick Start

### One-Command Build and Run

```bash
# From the tudat root directory
cmake -P wasm.cmake
```

This will:

1. Download and install Emscripten SDK (if needed)
2. Configure the WASM build
3. Build the visualization demo
4. Deploy WASM files to test directory
5. Start a local HTTP server on port 8832
6. Open your browser to `http://localhost:8832`

### Full Build (with API tests)

```bash
cmake -DFULL_BUILD=1 -P wasm.cmake
```

This includes:

- Full `tudatpy_wasm` Embind API module
- Node.js API and integration tests

### Manual Build Steps

If you prefer manual control:

```bash
# From the tudat root directory
cd build-wasm

# Configure with Emscripten (if not already done)
emcmake cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release

# Build the visualization demo
ninja tudat_wasm_web

# Or build the full Embind API
ninja tudatpy_wasm

# Start the server manually
cd ../tests/wasm/web
python3 start_server.py
```

## Build Targets

| Target              | Description                                    |
| ------------------- | ---------------------------------------------- |
| `tudat_wasm_web`    | Visualization demo for browser                 |
| `tudatpy_wasm`      | Full Embind API (mirrors Python bindings)      |
| `tudat_wasm_test`   | Node.js version for CLI testing                |
| `tudat_wasm_serve`  | Build web version and start HTTP server        |

## Web UI Features

The browser-based demo includes:

- **3D Orbital Visualization** using CesiumJS
  - Real-time ISS orbit from TLE data
  - Earth with accurate lighting and atmosphere
  - Orbit trail visualization
- **Force Model Comparison Chart**
  - J2-only vs Full Force Model divergence
  - Shows effect of J3/J4, Sun/Moon, drag, and SRP perturbations
  - Interactive time cursor synced with 3D view
- **Simulation Controls**
  - Play/pause and time warp controls
  - Configurable propagation duration
- **Console Output** with real-time logging
- **Dark space theme** with modern UI

## Directory Structure

```text
tests/wasm/
├── src/
│   └── wasmBindings.cpp  # Emscripten bindings for visualization
├── CMakeLists.txt        # Build configuration
├── README.md             # This file
├── test_wasm_api.js      # Node.js API tests
├── test_wasm_integration.js  # Integration tests
└── web/                  # Browser UI files
    ├── index.html        # Main HTML page
    ├── app.js            # JavaScript application
    └── start_server.py   # Development server
```

## Visualization Modules

The WASM bindings expose:

### Orbit Propagation

- `propagateKeplerOrbit()` - Analytical Keplerian propagation
- `propagateRK4Orbit()` - Fixed-step RK4 numerical integration
- `propagateRK78Orbit()` - Variable-step RK78 integration
- `propagateJ2vsFullForce()` - Compare J2-only vs full force model

### Astrodynamics

- `computeLibrationPoints()` - Earth-Moon L1-L5 points
- `propagateCR3BP()` - Circular Restricted 3-Body Problem
- `solveLambertProblem()` - Lambert targeting (Izzo algorithm)
- `computeAtmosphereDensity()` - Exponential atmosphere model

## Force Models

The full force model comparison includes:

| Perturbation | Description |
| ------------ | ----------- |
| J2 | Earth oblateness (dominant) |
| J3 | North-South asymmetry |
| J4 | Higher-order oblateness |
| Sun | Third-body gravitational |
| Moon | Third-body gravitational |
| Drag | Exponential atmosphere model |
| SRP | Solar radiation pressure |

## Requirements

- CMake 3.20+
- Git (for Emscripten SDK download)
- Python 3 (for web server)
- Modern web browser with WebAssembly support

The build script automatically downloads and installs Emscripten SDK 3.1.51.

## Port Configuration

The development server uses port **8832** (TUD on phone keypad + 2).

If the port is in use, either:

- Kill the existing server: `lsof -ti :8832 | xargs kill`
- Or modify `start_server.py` to use a different port

## Troubleshooting

### Build fails with "command not found: cmake"

Ensure CMake is in your PATH, or use the full path:

```bash
/opt/homebrew/bin/cmake -P wasm.cmake
```

### Server won't start (address in use)

```bash
lsof -ti :8832 | xargs kill
cmake -P wasm.cmake
```

### WASM module not loading in browser

Check browser console for errors. Ensure:

- The server is running with correct CORS headers
- SharedArrayBuffer is enabled (requires COOP/COEP headers)
- The .wasm file was deployed correctly

## Test Categories

The full API test suite covers:

1. **Basic Astrodynamics** - Unit conversions, physical constants
2. **Orbital Mechanics** - Keplerian/Cartesian conversions (NASA ODTBX benchmarks)
3. **Anomaly Conversions** - True/eccentric/mean anomaly
4. **Numerical Methods** - RK4/RK78 integration, interpolation
5. **Propagation** - CR3BP, two-body, mass propagation
6. **SPICE Interface** - Time conversions, frame rotations
7. **TLE/SGP4** - Two-line element propagation (Vallado benchmark)

## Accuracy Requirements

All tests meet or exceed the original Tudat native test requirements:

| Test              | Requirement    | WASM Result       |
| ----------------- | -------------- | ----------------- |
| TLE/SGP4 Position | < 50 m         | 15.6 m            |
| TLE/SGP4 Velocity | < 0.05 m/s     | 0.021 m/s         |
| Orbital Elements  | 1e-14 relative | Meets requirement |
