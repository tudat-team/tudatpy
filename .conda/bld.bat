mkdir build
if errorlevel 1 exit 1
cd build
if errorlevel 1 exit 1
cmake ^
    -G "Visual Studio 15 2017 Win64" ^
    -DCMAKE_CXX_STANDARD=17 ^
    -DBoost_NO_BOOST_CMAKE=ON ^
    -DCMAKE_PREFIX_PATH=%LIBRARY_PREFIX% ^
    -DCMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
    ..
if errorlevel 1 exit 1
cmake --build . --config RelWithDebInfo --target install
if errorlevel 1 exit 1
ctest
if errorlevel 1 exit 1
