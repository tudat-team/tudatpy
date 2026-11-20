# Installed C++ gravity consumer

This is a separate CMake project. It must use an actual installed Tudat package,
not build-tree exports or source include paths. Configure/build/install Tudat
with the normal dependency/compiler options and a dedicated prefix. For a
production-only installation, configure with `TUDAT_BUILD_TESTS=OFF`, then build
`kernel` before installing. This avoids requiring unbuilt test-support archives
from an installation configured to include tests.

For each backend setting, use separate build and install directories:

```bash
cmake -S . -B build-consumer-on \
  -DTUDAT_BUILD_WITH_OPENMP=ON -DTUDAT_BUILD_TESTS=OFF \
  -DCMAKE_INSTALL_PREFIX=/tmp/tudat-consumer-install-on
cmake --build build-consumer-on --target kernel --parallel 4
cmake --install build-consumer-on

cmake -S tests/installed_gravity_consumer -B /tmp/tudat-consumer-build-on \
  -DTudat_DIR=/tmp/tudat-consumer-install-on/lib/cmake \
  -DCMAKE_FIND_USE_PACKAGE_REGISTRY=OFF
cmake --build /tmp/tudat-consumer-build-on --parallel 2
ctest --test-dir /tmp/tudat-consumer-build-on --output-on-failure
```

Supply your dependency prefixes and compiler to both projects as usual. Repeat
with `TUDAT_BUILD_WITH_OPENMP=OFF` and distinct `off` paths. The consumer does not
add source/build include directories, OpenMP flags, or manual Eigen/Cereal
workarounds. It links installed targets. It checks agreement between the
installed config and target feature definitions, preservation of compiler
OpenMP support, finite acceleration, exact parity over 20 changing positions,
and rejection of multi-thread requests when the backend is unavailable. The
printed observed team size belongs to the internal loop helper used by the
installed headers, not instrumentation of the model's own force regions.

On 15 November 2026 both installed configurations configure, compile, link and
pass CTest with CMake 4.2 / GCC 14.4. The ON consumer observes four helper threads;
the OFF consumer stays sequential and rejects a request for two threads. Build
and link command inspection confirms that headers/libraries come from the two
install prefixes rather than the Tudat checkout. The installed prefixes were
not moved; this does not test package relocation.

This check exposed missing installed dependency discovery and propagation:
`FindCSpice.cmake` must ship with the config; Eigen and NRLMSISE-00 must be found;
Eigen must propagate through public targets; gravity headers must expose their
Cereal dependency. These are installation/export fixes, not numerical changes.
Local logs are under `local-profiles/20261115/installed-consumer/`.
