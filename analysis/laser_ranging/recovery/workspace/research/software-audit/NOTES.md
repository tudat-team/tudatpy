# Static inspection of independent laser-ranging implementations

Inspected 15 October 2026. No external solver was run. Exact inspected files,
retrieval times and SHA-256 hashes are in `SOURCES.json` next to this note.
These observations inform the design; they do not validate Tudat numerically.

## Orekit optical atmosphere

[MendesPavlisModel.java](https://raw.githubusercontent.com/CS-SI/Orekit/develop/src/main/java/org/orekit/models/earth/troposphere/MendesPavlisModel.java)
uses the wavelength-dependent hydrostatic and water-vapour terms and applies
one FCUL mapping factor to their sum. Pressure is supplied in Pa. Its mapping
coefficients use station latitude, altitude and measured temperature. This agrees
with the structure of the new Tudat implementation. Tudat reuses its existing
meteorological object and mapping correction rather than adopting Orekit's class
hierarchy. The native numerical controls come from the official IERS routines,
not from Orekit output. Water-vapour pressure conversion must be recorded when
relative humidity is supplied; different conversions are a possible independent
source of small differences.

## GROOPS station and range organization

[slrStation.cpp](https://raw.githubusercontent.com/groops-devs/groops/main/source/slr/slrStation.cpp)
keeps station displacement and epoch corrections separate from satellite geometry.
This supports using Tudat's existing station-motion machinery and a local target
correction instead of adding another general station class.

[slrObservation.cpp](https://raw.githubusercontent.com/groops-devs/groops/main/source/slr/slrObservation.cpp)
iterates distinct transmission, reflection and reception positions and applies
the reflector offset during the reflection-position calculation. It compares
one-way observed range with half the two geometric legs, then subtracts the
Earth gravitational delay. Its shown geometric iterations terminate at a
0.1 mm position change. The implementation is geocentric; it is not a direct
reference for Tudat's BCRS clock and coordinate transformations. Tudat instead
uses its existing coupled light-time solution, with explicit atomic interval
conversion and one shared target response for both retarded viewing directions.
The inspected GROOPS path also permits residual-dependent accuracy and later
normalization. Our fixed-orbit diagnostics retain raw residuals and supplied
uncertainties; no such reweighting is used.

## Consequences for this branch

The comparison supports the current local additions: raw CRD/SP3 data access,
one optical atmosphere correction, and one target-response correction. It gives
no reason to add another orbit integrator, estimator, generic model hierarchy or
external-solver dependency. Independent published formulas, pinned real inputs
and raw residual evidence remain the required checks.
