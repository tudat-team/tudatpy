mkdir build
if errorlevel 1 exit 1
cd build
if errorlevel 1 exit 1
cmake ^
    -G "Visual Studio 15 2017 Win64" ^
    -DCMAKE_CXX_STANDARD=17 ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DCMAKE_PREFIX_PATH=%LIBRARY_PREFIX% ^
    -DCMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
    -DPREFIX=%LIBRARY_PREFIX% ^
    -DTUDAT_CONDA_BUILD=on ^
    -DTUDAT_BUILD_STATIC_LIBRARY=off ^
    -DTUDAT_BUILD_TUDAT_TUTORIALS=off ^
    -DTUDAT_BUILD_WITH_SOFA_INTERFACE=on ^
    -DTUDAT_BUILD_WITH_SPICE_INTERFACE=on ^
    -DTUDAT_INSTALL=on ^
    -DTUDAT_TEST_INSTALL=off ^
    ..
if errorlevel 1 exit 1
cmake --build . --config RelWithDebInfo --target install
if errorlevel 1 exit 1
ctest
if errorlevel 1 exit 1
