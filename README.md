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

There are two directions you can go from here. CLion or the command line.

### CLion Build

6. Create a build profile in `File > Settings > Build, Execution, Deployment > CMake`. 
   - Note that the CMake configuration option `CMAKE_BUILD_TYPE` will be determined by the the build profile's `Build type` entry. A `Release` configuration
   will suppress a significant amount of harmless warnings during compilation. (*Currently: with the move to a 
   later version of boost, some warnings have cropped up that have either not been fixed in the source code, or 
   have not been suppressed via `tudat/cmake_modules/compiler.cmake`)

7. Add the CMake configuration to the `File > Settings > Build, Execution, Deployment > CMake > CMake options` text box:
   
```
-DCMAKE_PREFIX_PATH=<CONDA_PREFIX>
-DCMAKE_CXX_STANDARD=14
-DBoost_NO_BOOST_CMAKE=ON
```

> **Note** \
> The `CONDA_PREFIX` may be determined with by activating the environment installed in step 4 and printing its value:
> ````
> conda activate tudat-bundle && echo $CONDA_PREFIX
> ````

8. In the source tree on the left, right click the top level `CMakeLists.txt` then `Load/Reload CMake Project`.
   
9. `Build > Build Project`

### Command Line Build

6. Activate the environment installed in step 4

````
conda activate tudat-bundle
````

7. Run the `build.sh` script.

````
bash build.sh
````

## Verify your build

### Running `tudat` tests

1. Enter the `tudat` build directory
````
cd <build_directory>/tudat
````

2. Run the tests using `ctest` (packaged with CMake)
````
ctest
````

Desired result:
````
.. 
100% tests passed, 0 tests failed out of 224
Total Test time (real) = 490.77 sec
````

### Running `tudatpy` tests

1. Enter the `tudatpy` build directory
````
cd <build_directory>/tudatpy
````

2. Run the tests using `pytest`
````
pytest
````

Desired results:
````
=========================================== 6 passed in 1.78s ============================================
````

## Notes
- [**Windows Users ∩ CLion Users**] In CLion, be sure to set WSL as your Toolchain
  in `File>Settings>Build, Execution, Deployment>Toolchains`.

- [**CLion Users**] In CLion, the convention to set CMake arguments
  is to add them to `File>Settings>Build, Execution, Deployment>CMake Options`.
  
- [**All Users**] You can increase the number of cores used to compile `tudat` & `tudatpy` using the `-j<n>` 
  build argument, but **be aware** that the current complexity of the libraries can often result in your PC freezing indefinitely.

