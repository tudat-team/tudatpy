$buildScript = Join-Path $env:RUNNER_TEMP "build_tudat.cmd"
$buildCommands = @"
@echo off
call "%VCVARS64%" >nul 2>nul
set "CONDA_ENV=$env:CONDA_ENV"
set "PATH=%PATH%;%CONDA_ENV%\Library\bin;%CONDA_ENV%\Scripts;%CONDA_ENV%"
"%CONDA_ENV%\Library\bin\cmake.exe" --build "%GITHUB_WORKSPACE%\build" --config "$env:BUILD_TYPE" -j3
"@
Set-Content -Path $buildScript -Encoding ASCII -Value $buildCommands

# Conda's compiler activation hook prints the full Visual Studio discovery script
# for many Ninja rules. Hide that block while preserving all compiler output.
$filteringActivation = $false
& cmd.exe /D /E:ON /V:OFF /S /C "`"$buildScript`"" 2>&1 | ForEach-Object {
    $line = $_.ToString()
    if ($line -match '^[A-Z]:\\.*>SET DISTUTILS_USE_SDK=1\s*$') {
        $filteringActivation = $true
        return
    }
    if ($filteringActivation) {
        if ($line -match '^\[\d+/\d+\] ') {
            $filteringActivation = $false
            Write-Output $line
            return
        }
        if ($line -match '^\[vcvarsall\.bat\] Environment initialized for:') {
            $filteringActivation = $false
        }
        return
    }
    Write-Output $line
}
exit $LASTEXITCODE
