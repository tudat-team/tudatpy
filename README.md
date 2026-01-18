# Tudatpy

The **TU Delft Astrodynamics Toolbox in Python**, or **Tudatpy**, is a library that primarily exposes a powerful set of [C++  
libraries](https://tudat.tudelft.nl/) aiming at accelerating the implementation of simulations, real-data processing and analysis, and quality education in the field of Astrodynamics.
See the [documentation](https://tudat-space.readthedocs.io) for more.

For nominal usage, the use of our distributed **conda package** is recommended. For more details on the project, please refer to the [project website](https://docs.tudat.space/en/latest/) and the [project's Github page](https://github.com/tudat-team).

## Structure of the `Tudatpy` Repository

The `Tudatpy` repository contains both the source code and the binding code, together with the respective documentation and examples folders.
The next steps outline how to get to a working version of Tudatpy. First we list some prerequisites, and then we show how to set it up.

## Prerequisites

- [**Windows Users**] Windows Subsystem for Linux ([WSL](https://docs.microsoft.com/en-us/windows/wsl/install))
  - All procedures, including the following prerequisite, assume the use of WSL. Power users who wish to do otherwise,
    must do so at their own risk, with reduced support from the team.
  - Note that WSL is a, partially separated, Ubuntu terminal environment for Windows. Anaconda/Miniconda, Python and any other dependencies you require while **executing code** from the `Tudatpy` repository, must be installed in its Linux version via the Ubuntu terminal. This does not apply to PyCharm/CLion however, which can be configured to compile and/or run Python code through the WSL.
  - Note that, to access files and folders of WSL directly in Windows explorer, one can type `\\wsl$` or `Linux` in the Windows explorer access bar, then press enter.
  - At the opposite, please follow [this guide](https://docs.microsoft.com/en-us/windows/wsl/wsl2-mount-disk) to access Windows file trough WSL.
  - [This guide from Microsoft](https://docs.microsoft.com/en-us/windows/wsl/setup/environment) contains more information on the possibilities given trough WSL.
  - In the Ubuntu terminal environment under WSL, run the command `sudo apt-get install build-essential` to install the necessary compilation tools
- Anaconda/Miniconda installation ([Installing Anaconda](https://tudat-space.readthedocs.io/en/latest/_src_first_steps/tudat_py.html#installing-anaconda))
- CMake installation
  - Inside the Ubuntu terminal, install CMake by calling `sudo apt install cmake`.

## Setup

1. Clone the repository and enter directory

````
git clone https://github.com/tudat-team/tudatpy
cd tudatpy
````

2. Clone the `examples/tudatpy` submodule:

````
git submodule update --init --recursive
````

> **Note** \
> Submodules "allow you to keep a Git repository as a subdirectory of
> another Git repository" (from [the Git guide](https://git-scm.com/book/en/v2/Git-Tools-Submodules)). In particular,
> This "sub-repository" has its own branches and functions separately from `Tudatpy`. This is why the previous step is needed.

3. Switch `Tudatpy` to a new or an already existing branch using:

````
git checkout develop
````

> **Note**\
> Although you could virtually choose any branch, we recommend working with the `develop` branch, as it receives frequent updates and are the ones used to build the Conda packages.

4. Install the contained `environment.yaml` file to satisfy dependencies, then activate it:

````
conda env create -f environment.yaml
conda activate tudatpy-dev
````

> **Note**\
It is possible that the creation of the environment will 'time out'. A likely reason for this is that the packages required cannot be found by the current channel, `conda-forge`. It is then advisable to add the channel `anaconda` to ensure a proper creation of the environment.
>

5. Build TudatPy

```
python build.py -h                    # Show help and available flags
python build.py -j <number-of-cores>  # Compile Tudatpy
python build.py --tests               # [optional] To verify with ctest (see below)
```
This script compiles Tudatpy. It will take some time to execute, but you can speed up the process by increasing the number of cores used with the `-j` flag. If you wish to verify your installation with `ctest` (see below), add the `--tests` flag.
Once the project is built, all the build output is dumped by default in a directory called `build`, which is not tracked by Git.

6. Install

```
python install.py -h                 # Show help and available flags
python install.py -e                 # Install in "editable mode"
```

> **Note**\
> This script installs Tudatpy in your active conda environment. If you install with the `-e` flag, you will not have to re-install every time you update the source code of the library.
> And that's it! The next step shows you what to do if you want to uninstall the libraries.

7. Uninstall

```
python uninstall.py -h                # Show help and available flags
python uninstall.py                   # Uninstall Tudatpy
```
> **Note**\
> This script will remove Tudatpy from your Conda environment, but it will not delete the build directory.
>
>
## Verify your build

### Running `tudatpy` tests

1. Within the `tudatpy` directory, run `pytest`  (packaged with CMake)

````
pytest
````
Desired result:
````
=========================================== 6 passed in 1.78s ============================================
````
### Running `tudat` tests

2. Enter the `tudatpy/build` directory and run the tests using `ctest`
````
cd build
ctest
````

Desired result:
````
..
100% tests passed, 0 tests failed out of 224
Total Test time (real) = 490.77 sec
````
> **Note**\
> To speed up the tests, you can optionally use multiple cores as follows:
> `ctest -j <number_of_cores>`

## Building for WebAssembly

Tudatpy can be compiled to WebAssembly using Emscripten, enabling orbital mechanics simulations to run directly in web browsers.

### Prerequisites

- CMake 3.20+
- Git
- Ninja build system
- Node.js (for running tests)
- Python 3 (for web server)

### Quick Start

The simplest way to build for WASM is using the provided script:

```bash
cmake -P wasm.cmake
```

This single command will:
1. Download and install the Emscripten SDK (if not present)
2. Configure the WASM build
3. Build the WASM test suite
4. Run tests via Node.js
5. Start a local web server at http://localhost:8080

### Full Build Instructions

For development or more control over the build process, use the manual approach:

#### 1. Configure the Build

```bash
# Configure using the toolchain file (auto-downloads Emscripten SDK if needed)
cmake -B build-wasm -G Ninja \
  -DCMAKE_TOOLCHAIN_FILE=cmake_modules/toolchain-emscripten.cmake \
  -DCMAKE_BUILD_TYPE=Release
```

#### 2. Build the Project

```bash
# Build all WASM targets
cmake --build build-wasm

# Or use ninja directly for more control
cd build-wasm
ninja tudat_wasm_test
```

#### 3. Run the Test Suite

```bash
# Run the comprehensive WASM test suite (468 tests)
node build-wasm/tests/wasm/tudat_wasm_test.js
```

Expected output:
```
=== Tudat WASM Test Suite ===
[PASS] Unit conversion: degrees to radians
[PASS] Unit conversion: radians to degrees
...
=== Test Results ===
[INFO] Tests run:    468
[INFO] Tests passed: 468
[INFO] Tests failed: 0
[PASS] *** ALL TESTS PASSED ***
```

### Test Suite Coverage

The WASM test suite (`tests/wasm/`) validates core Tudat functionality across multiple domains:

| Category | Tests | Description |
|----------|-------|-------------|
| Basic Astrodynamics | 50+ | Unit conversions, constants, Eigen operations |
| Orbital Elements | 30+ | Keplerian, Cartesian, spherical conversions |
| Propagation | 20+ | CR3BP, two-body, mass, custom state propagation |
| SPICE Interface | 15+ | Time conversions, frame rotations, TLE propagation |
| Gravitation | 25+ | Spherical harmonics, libration points, Jacobi energy |
| Aerodynamics | 15+ | Exponential and NRLMSISE-00 atmosphere models |
| Mission Segments | 20+ | Lambert targeting, gravity assists, escape/capture |
| Electromagnetism | 15+ | Radiation pressure, luminosity models |
| Integrators | 15+ | RK4, RK78, RKF45, RKDP87, Adams-Bashforth-Moulton |
| Ephemerides | 15+ | Kepler, tabulated, constant, rotational ephemerides |
| Earth Orientation | 30+ | Time scales, EOP, polar motion, leap seconds |

### Rebuilding After Changes

For incremental rebuilds during development:

```bash
cd build-wasm
ninja tudat_wasm_test && node tests/wasm/tudat_wasm_test.js
```

For a clean rebuild (required after toolchain or CMake changes):

```bash
rm -rf build-wasm
cmake -B build-wasm -G Ninja \
  -DCMAKE_TOOLCHAIN_FILE=cmake_modules/toolchain-emscripten.cmake \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-wasm
```

### Configuration Options

Specify a different Emscripten version:
```bash
cmake -B build-wasm -DCMAKE_TOOLCHAIN_FILE=cmake_modules/toolchain-emscripten.cmake -DEMSDK_VERSION=3.1.52
```

### Managing Emscripten SDK

Update the Emscripten SDK:
```bash
cmake --build build-wasm --target update-emscripten
```

List available Emscripten versions:
```bash
cmake --build build-wasm --target list-emscripten-versions
```

### External Dependencies

All dependencies are automatically handled for WASM builds:

| Dependency | Notes |
|------------|-------|
| **Boost** | Headers automatically provided via Emscripten port |
| **Eigen3** | Automatically downloaded if not found |
| **CSpice** | Automatically downloaded and built with Emscripten |
| **nrlmsise00** | Automatically downloaded and built with Emscripten |
| **SOFA** | Required, automatically fetched and built |

Dependencies are downloaded to the `_deps/` folder inside the build directory.

### Data Files

WASM builds embed essential data files directly into the binary. These are mounted at `/tudat_data/` in Emscripten's virtual filesystem and include:

- Earth orientation parameters (EOP)
- Leap second tables
- NRLMSISE-00 atmosphere tables
- Gravity field models
- SPICE kernels (LSK, PCK)

For custom data files in your JavaScript code:

```javascript
// Access embedded data
const eopPath = '/tudat_data/earth_orientation/';

// Or mount additional data
Module.FS.mkdir('/custom_data');
Module.FS.writeFile('/custom_data/my_file.txt', myData);
```

### Web Server Testing

To test WASM builds in a browser environment:

#### Using the Quick Start Script

```bash
cmake -P wasm.cmake
# Opens browser at http://localhost:8080
```

#### Manual Web Server

```bash
# Start a Python web server in the build directory
cd build-wasm
python3 -m http.server 8080

# Or with specific CORS headers for SharedArrayBuffer
python3 -c "
from http.server import HTTPServer, SimpleHTTPRequestHandler
class Handler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Cross-Origin-Opener-Policy', 'same-origin')
        self.send_header('Cross-Origin-Embedder-Policy', 'require-corp')
        super().end_headers()
HTTPServer(('', 8080), Handler).serve_forever()
"
```

Then open `http://localhost:8080/tests/wasm/` in your browser.

### Build Targets

| Target | Description |
|--------|-------------|
| `tudat_wasm_test` | Main test executable (468 tests) |
| `tudat` | Core Tudat library compiled for WASM |
| `cspice` | CSPICE library compiled for WASM |
| `nrlmsise00` | NRLMSISE-00 atmosphere library |
| `sofa` | SOFA astronomy library |

### Troubleshooting

**Build fails with missing Emscripten:**
The toolchain file automatically downloads and configures the Emscripten SDK. Just ensure you're using the toolchain:
```bash
cmake -B build-wasm -G Ninja \
  -DCMAKE_TOOLCHAIN_FILE=cmake_modules/toolchain-emscripten.cmake
```

**Tests crash with memory errors:**
```bash
# Increase memory limit (default is sufficient for most tests)
cmake -B build-wasm -G Ninja \
  -DCMAKE_TOOLCHAIN_FILE=cmake_modules/toolchain-emscripten.cmake \
  -DCMAKE_EXE_LINKER_FLAGS="-s INITIAL_MEMORY=256MB -s ALLOW_MEMORY_GROWTH=1"
```

**Browser shows CORS errors:**
Use the CORS-enabled Python server shown above, or configure your web server to send appropriate headers.

### Notes

- The Emscripten SDK is installed to `.emsdk/` in the project root (gitignored)
- Standard Boost.Test-based tests and tutorials are disabled for WASM builds
- A lightweight WASM-specific test suite is built instead using `wasmTestFramework.h`
- Data files are embedded at compile time and available at `/tudat_data/`
