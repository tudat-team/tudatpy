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
