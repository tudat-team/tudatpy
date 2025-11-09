# Solar Longitude Debug Investigation

## Purpose
This debug script investigates why the Stokes coefficient method produces different densities than the polynomial method at non-zero solar longitudes.

## What It Does
- Tests at solar longitude 90° with a single query point (radius=8km, body_lon=0°, body_lat=0°)
- Creates both polynomial and Stokes ComaModels with identical rotation functions
- Prints detailed debug output showing:
  - Solar longitude calculation from rotation matrix
  - Stokes dataset longitude grid
  - SH coefficients from polynomial evaluation
  - SH coefficients from Stokes interpolation
  - Final density values from both methods

## How to Build and Run

```bash
cd /Users/markusreichel/PhD/tudatpy/cmake-build-debug
make ComaAnalysis_debug_solar_longitude
./bin/ComaAnalysis_debug_solar_longitude > ../examples/tudat/tudat/applications/coma_analysis/debug_solar_longitude/debug_output.txt 2>&1
```

## Debug Output Added to Code

Debug output has been added to:
1. **comaModel.cpp::calculateSolarLongitude()** - Shows rotation matrix, sun direction, and solar longitude computation
2. **comaModel.cpp::computeNumberDensityFromPolyCoefficients()** - Shows polynomial evaluation and resulting SH coefficients
3. **comaModel.cpp::computeNumberDensityFromStokesCoefficients()** - Shows interpolation point, Stokes grid, and interpolated coefficients

## What to Look For

Compare the debug output sections:
1. Check if `calculateSolarLongitude()` returns the same value for both methods
2. Compare the SH coefficients (C(0,0), C(1,0), C(1,1), S(1,1)) between polynomial and Stokes
3. If coefficients differ, this is where the bug is
4. Check if the solar longitude for interpolation matches the Stokes dataset grid

## Expected Behavior
If working correctly:
- Both methods should calculate the same solar longitude (90° = 1.5708 rad)
- The Stokes interpolation should use 1.5708 rad as the interpolation coordinate
- Both methods should produce identical SH coefficients
- Both methods should produce identical final densities

## Cleanup
To remove debug output after investigation, search for "DEBUG" in comaModel.cpp and remove those std::cout statements.
