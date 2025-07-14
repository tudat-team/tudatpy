mkdir build
if errorlevel 1 exit 1
cd build
if errorlevel 1 exit 1
cmake ^
    -G "MinGW Makefiles" ^
    -DCMAKE_CXX_STANDARD=14 ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DCMAKE_PREFIX_PATH=%LIBRARY_PREFIX% ^
    -DCMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
    -DPREFIX=%LIBRARY_PREFIX% ^
    -DTUDAT_BUILD_STATIC_LIBRARY=on ^
    -DTUDAT_BUILD_TUDAT_TUTORIALS=off ^
    -DTUDAT_BUILD_TESTS=off ^
    -DTUDAT_BUILD_WITH_SOFA_INTERFACE=on ^
    -DTUDAT_BUILD_WITH_SPICE_INTERFACE=on ^
    -DTUDAT_BUILD_WITH_JSON_INTERFACE=off ^
    ..
if errorlevel 1 exit 1
cmake --build . --config RelWithDebInfo --target install -- -j2
if errorlevel 1 exit 1
:: ctest --verbose
if errorlevel 1 exit 1
