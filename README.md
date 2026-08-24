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
- Anaconda/Miniconda installation ([Installing Anaconda](https://docs.tudat.space/en/latest/getting-started/use-of-tools/conda.html))
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

5. Install `pre-commit` hooks

This repository uses [pre-commit hooks](https://pre-commit.com) to automatically apply consistent formatting to all C++ and Python files.

Run 

```
pre-commit install
```

to install the pre-commit hooks.
After this, anything you commit will be automatically formatted using `clang-format` and `black`, without requiring your attention.

6. Build TudatPy

```
python build.py -h                   # Show help and available flags
python build.py -j <number-of-cores>  # Compile Tudatpy
```
This script compiles Tudatpy. It will take some time to execute, but you can speed up the process by increasing the number of cores used with the `-j` flag.
Once the project is built, all the build output is dumped by default in a directory called `build`, which is not tracked by Git.

### Configuring the high-precision C++ state scalar

Tudat's high-precision state interfaces use `tudat::HighPrecisionStateScalar`,
declared in the generated `tudat/config.hpp` header. It is Boost's
quad-precision binary float by default:

```sh
cmake -S . -B build
cmake --build build
```

The supported CMake values are `CPP_BIN_FLOAT_QUAD` (the default) and
`LONG_DOUBLE`. This selection affects the C++ state API and its binary interface,
so all libraries and C++ consumers must use the same configuration.
When TudatPy is built with `TUDATPY_BUILD_WITH_QUAD_PRECISION_EXPOSURE=ON`
(the default), Python can select the C++ state scalar before importing the
kernel:

```python
import tudatpy

if tudatpy.quad_precision_available():
    tudatpy.set_state_scalar_type("quad")

from tudatpy import kernel
```

The default Python mode is `"double"`. The selection is process-wide and is
locked by the first import of `tudatpy.kernel`; calling
`set_state_scalar_type` afterwards raises an error. The active selection is
reported by `tudatpy.get_state_scalar_type()`.

The Boost scalar is not exposed as a Python numeric type. Python scalar and
NumPy inputs and outputs remain binary64; selecting `"quad"` instantiates the
supported internal C++ state, propagation, observation, and estimation paths
with `HighPrecisionStateScalar`.

Eigen typedefs ending in `ld` always use literal `long double`; typedefs ending
in `hps` use the configured `HighPrecisionStateScalar`. Quad state arithmetic
does not make every input or external interface quad precision. In particular,
`Time` retains its split representation with a `long double` fractional field,
and data or APIs defined as `double` remain limited to that precision. The
effective epoch resolution therefore depends on the platform's `long double`
implementation (commonly extended precision on Linux/x86, but equivalent to
`double` with MSVC).

7. Install

```
python install.py -h                 # Show help and available flags
python install.py -e                 # Install in "editable mode"
```

> **Note**\
> This script installs Tudatpy in your active conda environment. If you install with the `-e` flag, you will not have to re-install every time you update the source code of the library.
> And that's it! The next step shows you what to do if you want to uninstall the libraries.

8. Uninstall

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
Use `pytest --no-remote-data` to skip tests that require a network connection.

Desired result:
````
=========================================== 6 passed in 1.78s ============================================
````
### Running `tudat` tests

Note that `tudat` tests are only built when using the `--tests` flag with `build.py`, for example `python build.py --tests -j4`.

2. Enter the `tudatpy/build` directory and run the tests using `ctest`
````
cd build
ctest -j <number-of-cores>
````

Desired result:
````
..
100% tests passed, 0 tests failed out of 224
Total Test time (real) = 490.77 sec
````

Note that when running tests in parallel with `-j`, CTest may execute tests in a non-sequential order to minimize total execution time.
