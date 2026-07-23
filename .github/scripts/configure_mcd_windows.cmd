@echo off
setlocal

for /f "usebackq delims=" %%i in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSINSTALL=%%i"
if "%VSINSTALL%"=="" exit /b 1

set "VCVARS64=%VSINSTALL%\VC\Auxiliary\Build\vcvars64.bat"
if not exist "%VCVARS64%" exit /b 1
echo VCVARS64=%VCVARS64%>>"%GITHUB_ENV%"

call "%VCVARS64%" >nul 2>nul
set "PATH=%PATH%;%CONDA_ENV%\Library\bin;%CONDA_ENV%\Scripts;%CONDA_ENV%"
set "CMAKE_PREFIX_PATH=%CONDA_ENV%\Library"

"%CONDA_ENV%\Library\bin\cmake.exe" -B "%GITHUB_WORKSPACE%\build" -S "%GITHUB_WORKSPACE%" -G Ninja ^
  -DCMAKE_BUILD_TYPE=%BUILD_TYPE% ^
  -DCMAKE_C_COMPILER=cl ^
  -DCMAKE_CXX_COMPILER=cl ^
  -DCMAKE_Fortran_COMPILER="%CONDA_ENV%\Library\bin\flang-new.exe" ^
  -DTUDAT_BUILD_GITHUB_ACTIONS=%1 ^
  -DTUDAT_FORCE_DYNAMIC_RUNTIME=ON ^
  -DTUDAT_BUILD_WITH_TEST_PCH=ON ^
  -DTUDAT_BUILD_WITH_PCH=ON ^
  -DTUDAT_BUILD_WITH_MCD_INTERFACE=ON ^
  -DTUDAT_BUILD_FOR_REDUCED_COMPILE_TIME=ON ^
  -DCMAKE_UNITY_BUILD=ON
