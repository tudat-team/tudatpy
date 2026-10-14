# Raw EOP reader fixtures and standalone verification

These small, unchanged extracts contain real IERS C04 and USNO finals2000A
records. `sources.json` records exact URLs, retrieval dates, complete downloaded
file SHA-256 hashes, selected rows, and extract hashes. Original header and
fixed-width padding is intentionally retained, with a local whitespace attribute. Complete downloads remain
outside the repository under the local reference workspace. The C04 14 archive
is the legacy 16-column format; the generic `eopc04.1962-now` URL now serves
modern 21-column C04 and should not be treated as a stable legacy-format URL.
The 12h C04 metadata identifies ITRF 2020; no more specific realization is inferred
from the download's filename.

The tests check published columns and independently specified SI conversions.
Normal tests are offline. Five original data rows from finals2000A include three
observed rows, one prediction with missing LOD and CPO, and a future date-only
placeholder. Only Bulletin A columns are used. Bulletin B values are never used
as replacements. Completely blank data portions of full-width future records are
skipped; they do not extend raw coverage. Partial mandatory values, truncated
fixed-width records, unpaired optional coordinates, non-finite numbers, invalid
calendar dates, inconsistent MJD/hour, and duplicate/unsorted epochs raise errors
with source and line number. Optional values stay missing rather than becoming
zero. C04 records with dPsi/dEps headers are rejected rather than silently
reinterpreted as dX/dY.

The raw reader target uses only standard C++17. It has no Tudat dynamics,
estimation, interpolation, SOFA, Eigen, Boost, or serialization library dependency.
Boost headers are needed by the test executable only. Build its isolated test
from the repository root with one thread:

```sh
cmake -S tests/test_tudat/src/io/earth_orientation_fixtures/standalone -B /absolute/local/eop-build
cmake --build /absolute/local/eop-build --parallel 1
ctest --test-dir /absolute/local/eop-build --output-on-failure
```

If needed, pass `-DBOOST_TEST_INCLUDE_DIR=/path/to/boost/include` during configure.
The same C++ suite is registered in the normal Tudat test tree.

The compatibility EOPReader delegates parsing to this layer but temporarily keeps
its historical one-hour buffer. The raw records retain exact coverage and no
artificial epochs. Python exposes these records only through
`data_input.environment_data.earth_orientation`; it constructs no physical model
or interpolator.
