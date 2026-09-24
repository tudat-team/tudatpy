# Tudatpy

The **TU Delft Astrodynamics Toolbox in Python**, or **Tudatpy**, is a library that primarily exposes a powerful set of [C++  
libraries](https://tudat.tudelft.nl/) aiming at accelerating the implementation of simulations, real-data processing and analysis, and quality education in the field of Astrodynamics.
See the [documentation](https://tudat-space.readthedocs.io) for more.

For nominal usage, the use of our distributed **conda package** is recommended. For more details on the project, please refer to the [project website](https://docs.tudat.space/en/latest/) and the [project's Github page](https://github.com/tudat-team).

## Structure of the `Tudatpy` Repository

The `Tudatpy` repository contains both the source code and the binding code, together with the respective documentation and examples folders.
The next steps outline how to get to a working version of Tudatpy. First we list some prerequisites, and then we show how to set it up.

## Prerequisites

- ``conda``: You must have ``conda`` installed on your system to obtain all required dependencies. See our [user guide](https://docs.tudat.space/en/latest/getting-started/use-of-tools/conda.html) for an introduction.
- **Windows Users**: While local builds on Windows are possible, note that support from the core developer team is limited.
  - Install [Visual Studio 2022, Version 17](https://learn.microsoft.com/en-us/visualstudio/releases/2022/release-history#evergreen-bootstrappers)
  - As an alternative, Windows Subsystem for Linux ([WSL](https://docs.microsoft.com/en-us/windows/wsl/install), see also our [user guide](https://docs.tudat.space/en/latest/getting-started/use-of-tools/windows-subsystem-for-linux.html)) can be installed for a Linux environment inside Windows.

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

7. Install

```
python install.py -h                 # Show help and available flags
python install.py -e                 # Editable development installation
python install.py                    # Frozen installation of the current build
```

If you are using Windows, you might have to run the commands from a shell with admin privileges, since it modifies files inside your conda environment.

> **Note**\
> This script installs Tudatpy in your active conda environment. Editable mode links the Python files in the environment directly to this source checkout. Source edits and branch switches therefore affect the installed package immediately, while the compiled kernel remains the one in the selected build directory. Use editable mode for active development and keep the checkout on a compatible revision.
>
> Run `python install.py` without `-e` when the environment must remain independent of later source edits or branch switches. This installs a snapshot of the current build through CMake. Unlike editable installations, regular installations are not tracked by `uninstall.py` and must currently be removed manually.

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
