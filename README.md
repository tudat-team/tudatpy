# tudat-bundle

This repository facilitates parallel development between the `tudat` (C++) and the
`tudatpy` (Python) library. 
Specific indications for documenting `tudat` or `tudapy` are reported in the `tudat-multidoc/README.md` file.


## Structure of the `tudat-bundle`

The `tudat-bundle` comprises the following repositories:

- `tudat`, where the tudat source code is located (this is a separate git repository);
- `tudatpy`, where the tudatpy binding code is located (this is a separate git repository);
- `tudat-multidoc`, where the documentation and the system to build the API is located (this is a separate git repository);
- `cli`, where the Python Command Line Interface scripts to build the API are located;

In addition, once the project is built, all the build output will be dumped in the `cmake-build-debug` directory, which
is not tracked by Git. If the API is also built, more untracked directories will appear, but this is explained in the
`tudat-multidoc/README.md` file.

## Prerequisites

- [**Windows Users**] Windows Subsystem for Linux ([WSL](https://docs.microsoft.com/en-us/windows/wsl/install))
  - All procedures, including the following prerequisite, assume the use of WSL. Power users who wish to do otherwise,
    must do so at their own risk, with reduced support from the team.
  - Note that WSL is a, partially separated, Ubuntu terminal environment for Windows. Anaconda/Miniconda, Python and any other dependencies you require while **executing code** from the `tudat-bundle`, must be installed in its Linux version via the Ubuntu terminal. This does not apply to PyCharm/CLion however, which can be configured to compile and/or run Python code through the WSL.
  - Note that, to access files and folders of WSL directly in Windows explorer, one can type `\\wsl$` or `Linux` in the Windows explorer access bar, then press enter.
  - At the opposite, please follow [this guide](https://docs.microsoft.com/en-us/windows/wsl/wsl2-mount-disk) to access Windows file trough WSL.
  - [This guide from Microsoft](https://docs.microsoft.com/en-us/windows/wsl/setup/environment) contains more information on the possibilities given trough WSL.
- Anaconda/Miniconda installation ([Installing Anaconda](https://tudat-space.readthedocs.io/en/latest/_src_first_steps/tudat_py.html#installing-anaconda))

## Setup

1. Clone the repository and enter directory

````
git clone https://github.com/tudat-team/tudat-bundle
cd tudat-bundle
````

> **Note** \
> The `tudat-bundle` repository uses git submodules, which "allow you to keep a Git repository as a subdirectory of 
> another Git repository" (from [the Git guide](https://git-scm.com/book/en/v2/Git-Tools-Submodules)). In particular,
> in the `tudat-bundle` there are four different subdirectories that are separate repositories: `tudat`, `tudatpy`,
> `tudat-multidoc` and `tudat-multidoc/multidoc`. Each repository has its own branches and functions separately from 
> the others. This is the reason why the following two steps are needed.

2. Clone the `tudat` & `tudatpy` submodules

````
git submodule update --init --recursive
````

3. [Optional] Switch `tudat` & `tudatpy` to their desired branches using

````
cd <tudat/tudatpy>
git checkout <branch-name>
````
Be advised that the branch from with the Conda packages are built, and that is being maintained the most, is `develop`.
See [here](https://github.com/tudat-team/tudatpy/tree/develop) for tudatpy develop branch, and [here](https://github.com/tudat-team/tudat/tree/develop) for tudat develop branch.

It is then recommended to switch to the `develop` branch using the commands above.

4. Install the contained `environment.yaml` file to satisfy dependencies

````
conda env create -f environment.yaml
````

> **Note**
>  Are you a Mac user with a M1 processor (a osx-arm64 system)? Then, edit the channels of the `environment.yaml` file as follows, before creating the environment:
>  ````
>  channels:
>  - https://conda.anaconda.org/conda-forge/osx-arm64
>  - https://conda.anaconda.org/tudat-team/osx-arm64
>  ````

There are two directions you can go from here. CLion or the command line.

### Build: CLion
> **Note**
> - [**Windows Users ∩ CLion Users**] In CLion, be sure to set WSL as your Toolchain
>  in `File>Settings>Build, Execution, Deployment>Toolchains`.
>
> - [**CLion Users**] In CLion, the convention to set CMake arguments
>  is to add them to `File>Settings>Build, Execution, Deployment>CMake Options`.

5. Open CLion, create a new project from `File > New Project` and select the directory that has been cloned under bullet
   point 1 (named `tudat-bundle`).
> **Note** \
> To avoid issues with CLion, the directory of the project should correspond exactly to the cloned directory named  `tudat-bundle`.


6. Create a build profile in `File > Settings > Build, Execution, Deployment > CMake`.
> **Note** \
> The CMake configuration option `CMAKE_BUILD_TYPE` will be determined by the the build profile's `Build type` entry.
> A `Release` configuration will suppress a significant amount of harmless warnings during compilation. Currently,
> with the move to a later version of boost, some warnings have cropped up that have either not been fixed in the
> source code, or have not been suppressed via `tudat/cmake_modules/compiler.cmake`.

7. Add the CMake configuration to the `File > Settings > Build, Execution, Deployment > CMake > CMake options` text box:

```
-DCMAKE_PREFIX_PATH=<CONDA_PREFIX>
-DCMAKE_CXX_STANDARD=14
-DBoost_NO_BOOST_CMAKE=ON
```

[**Optional**] Also add the following line to the `File > Settings > Build, Execution, Deployment > CMake > CMake options` text box to to build tudatpy with the tests.
````
-DTUDAT_BUILD_TESTS="${build_tests:-1}"
````

The following line can also be edited if you wish to build tudatpy with its debug info (switching from `Release` to `RelWithDebInfo`; note that `Debug` is also available):
````
-DCMAKE_BUILD_TYPE=RelWithDebInfo
````


> **Note** \
> The `CONDA_PREFIX` may be determined with by activating the environment installed in step 4 and printing its value:
> ````
> conda activate tudat-bundle && echo $CONDA_PREFIX
> ````

[**Optional**] Add `-j<n>` to `File > Settings > Build, Execution, Deployment > CMake > Build options` to use multiple
processors. It is likely that if you use all of your processors, your build will freeze your PC indefinitely. It is
recommended to start at `-j2` and work your way up with further builds, ensuring **no unsaved work** in the background.

8. In the source tree on the left, right click the top level `CMakeLists.txt` then `Load/Reload CMake Project`.

9. `Build > Build Project`

### Build: Command line

5. Activate the environment installed in step 4

````
conda activate tudat-bundle
````

[**Optional**] Edit the [build.sh](build.sh) script to build tudatpy with the tests by changing the `BUILD_TESTS` variable:
````
BUILD_TESTS="${build_tests:-1}"
````

The following line can also be edited if you wish to build tudatpy with its debug info (switching from `Release` to `RelWithDebInfo`; note that `Debug` is also available):
````
-DCMAKE_BUILD_TYPE=RelWithDebInfo
````

As building can take a while, you can build using multiple processors by appending `-j4` to the `cmake --build .` command as seen below, where 4 can be any number that you think your machine can handle. Do note that this can take up a few GB of RAM per processor used, so be aware of the other processes on your machine to avoid freezing or crashing. 
````
cmake --build . -j4
````

6. Run the [build.sh](build.sh) script.

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

Desired result:
````
=========================================== 6 passed in 1.78s ============================================
````

## Use your build
The path of the TudatPy kernel that has been manually compiled needs to be added before importing any `tudatpy.kernel` module.
This can be done with the following two lines, with `<kernel_path>` being similar to `<tudat-bundle_path>/build/tudatpy`:
```
import sys
sys.path.insert(0, <kernel_path>)
```

## Notes

- [**All Users**] You can increase the number of cores used to compile `tudat` & `tudatpy` using the `-j<n>`
  build argument, but **be aware** that the current complexity of the libraries can often result in your PC freezing indefinitely.

