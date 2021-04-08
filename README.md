# tudat-bundle

This repository facillitates parallel development between the `tudat` (C++) and the
`tudatpy` (Python) library.

## Prerequisites

- [**Windows Users**] Windows Subsystem for Linux ([WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10))
  - All procedures, including the following prerequisite, assume the use of WSL. Power users who wish to do otherwise,
    must do so at their own risk, with reduced support from the team.
- Anaconda/Miniconda installation ([Installing Anaconda](https://tudat-space.readthedocs.io/en/latest/_src_first_steps/tudat_py.html#installing-anaconda))

## Setup

1. Clone the repository and enter directory

````
git clone --single-branch --branch minimal https://github.com/tudat-team/tudat-bundle
cd tudat-bundle
````

2. Clone the `tudat` & `tudatpy` submodules

````
git submodule update --init --recursive
````

3. [Optional] Switch `tudat` & `tudatpy` to their desired branches using

````
cd <tudat/tudatpy>
git checkout <branch-name>
````

4. Install the contained `environment.yaml` file to satisfy dependencies

````
conda env create -f environment.yaml
````

5. Activate the environment installed in step 4

````
conda activate tudat-bundle
````

6. Determine your `CONDA_PREFIX` path

````
echo $CONDA_PREFIX
````

7. Set the following CMake build configuration (See Notes below)

````
-DCMAKE_PREFIX_PATH=<CONDA_PREFIX>
-DCMAKE_CXX_STANDARD=14
-DBoost_NO_BOOST_CMAKE=ON
````

Alternatively, add the following to the `CMakeLists.txt` (extra vigilance required when committing changes):

````
set(CMAKE_PREFIX_PATH <CONDA_PREFIX>)
set(CMAKE_CXX_STANDARD 14)
set(Boost_NO_BOOST_CMAKE ON)
````

## Notes
- [**Windows ∩ CLion Users**] In CLion, be sure to set WSL as your Toolchain
  in `File>Settings>Build, Execution, Deployment>Toolchains`.

- [**CLion Users**] In CLion, the convention to set CMake arguments
  is to add them to `File>Settings>Build, Execution, Deployment>CMake Options`.

