# tudat-bundle

This repository facilitates parallel development between the `tudat` (C++) and the
`tudatpy` (Python) library.
Specific indications for documenting `tudat` or `tudapy` are reported in the `tudat-multidoc/README.md` file.

Please note that for normal usage, we recommend the use of our **conda packages**. For more details on the project, we refer to the [project website](https://docs.tudat.space/en/latest/) and our [project Github page](https://github.com/tudat-team).

## Structure of the `tudat-bundle`

The `tudat-bundle` comprises the following repositories:

- `tudat`, where the tudat source code is located (this is a separate git repository);
- `tudatpy`, where the tudatpy binding code is located (this is a separate git repository);
<!-- - `cli`, where the Python Command Line Interface scripts to build the API are located; -->

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
  - In the Ubuntu terminal environment under WSL, run the command `sudo apt-get install build-essential` to install the necessary compilation tools
- Anaconda/Miniconda installation ([Installing Anaconda](https://tudat-space.readthedocs.io/en/latest/_src_first_steps/tudat_py.html#installing-anaconda))
- CMake installation
  - Inside the Ubuntu terminal, install CMake by calling `sudo apt install cmake`.

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

3. Switch `tudat` & `tudatpy` to their desired branches using

````
cd tudat
git checkout develop
cd ..
cd tudatpy
git checkout develop
````
We recommend working with the `develop` branches ([tudatpy/develop](https://github.com/tudat-team/tudatpy/tree/develop), [tudat/develop](https://github.com/tudat-team/tudat/tree/develop)), as they receive frequent updates and are the ones used to build the Conda packages.

4. Install the contained `environment.yaml` file to satisfy dependencies

````
conda env create -f environment.yaml
````

It is possible that the creation of the environment will 'time out'. A likely reason for this is that the packages required cannot be found by the current channel, `conda-forge`. It is then advisable to add the channel `anaconda` to ensure a proper creation of the environment.

There are two directions you can go from here: the command line or CLion.

### Build & Install: Command line

5. Activate the conda environment

````
conda activate tudat-bundle
````

6. Build Tudat and TudatPy

```
python build.py -h                   # Show help and available flags
python build.py -j<number-of-cores>  # Compile tudat and tudatpy
```
This script compiles tudat and tudatpy. It will take some time to execute, but you can speed up the process by increasing the number of cores used with the `-j` flag.

7. Install

```
python install.py -h                 # Show help and available flags
python install.py -e                 # Install in "editable mode"
```
This script installs tudat and tudatpy in your active conda environment. If you install with the `-e` flag, you will not have to re-install every time you update the source code of the libraries.

And that's it! The next step shows you what to do if you want to uninstall the libraries.

8. Uninstall

```
python uninstall.py -h                # Show help and available flags
python uninstall.py                   # Uninstall tudat and tudatpy
```
This script will remove tudat and tudatpy from your Conda environment, but it will not delete the build directory.

### Alternative build: CLion
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

The `CONDA_PREFIX` may be determined by activating the environment installed in step 4 and printing its value:
````
conda activate tudat-bundle && echo $CONDA_PREFIX
````

The following line can also be edited if you wish to build tudatpy with its debug info (switching from `Release` to `RelWithDebInfo`; note that `Debug` is also available):
````
-DCMAKE_BUILD_TYPE=RelWithDebInfo
````

[**Optional**] Add `-j<n>` to `File > Settings > Build, Execution, Deployment > CMake > Build options` to use multiple
processors. It is likely that if you use all of your processors, your build will freeze your PC indefinitely. It is
recommended to start at `-j2` and work your way up with further builds, ensuring **no unsaved work** in the background.

8. In the source tree on the left, right click the top level `CMakeLists.txt` then `Load/Reload CMake Project`.

9. `Build > Build Project`

### Alternative build: VSCode

This section explains how to configure VSCode with CMake presets to build and manage your `tudat` build. The configuration supports parallel builds, switching between Debug and Release modes, enabling/disabling tests, and cleaning build directories.

#### Requirements

1. **VSCode** with the following extensions:
   - [CMake Tools](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools)
   - [C++ Extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools)

2. **CMake** installed on your system (minimum version 3.19).

3. **Conda** environment for managing dependencies.

4. **Make/Ninja** installed based on your operating system for build generation.

#### Steps to Configure VSCode

1. **Create the `CMakePresets.json` file** in the root directory of your project. Use the following template with necessary adjustments:

```json
{
  "version": 3,
  "cmakeMinimumRequired": {
    "major": 3,
    "minor": 19,
    "patch": 0
  },
  "configurePresets": [
    {
      "name": "default",
      "hidden": false,
      "generator": "Unix Makefiles", // Set this based on the OS: "Unix Makefiles" for Linux/macOS, "Ninja" or "Visual Studio" for Windows
      "binaryDir": "${sourceDir}/build", // Build directory
      "cacheVariables": {
        "CMAKE_PREFIX_PATH": "/path/to/your/conda/environment", // <--- CHANGE THIS: Set this to your conda environment or other toolchain path
        "CMAKE_INSTALL_PREFIX": "/path/to/your/conda/environment", // <--- CHANGE THIS
        "CMAKE_CXX_STANDARD": "14", // The C++ standard to use
        "Boost_NO_BOOST_CMAKE": "ON", // Ensures that the system-installed Boost CMake files are not used
        "CMAKE_BUILD_TYPE": "Release", // Default to Release build
        "TUDAT_BUILD_TESTS": "ON" // Set to ON to build tests by default
      }
    },
    {
      "name": "debug",
      "inherits": "default",
      "description": "Debug build configuration",
      "cacheVariables": {
        "CMAKE_BUILD_TYPE": "Debug" // Use Debug mode for this preset
      }
    },
    {
      "name": "no-tests",
      "inherits": "default",
      "description": "Build without tests",
      "cacheVariables": {
        "TUDAT_BUILD_TESTS": "OFF" // Disable test building in this preset
      }
    }
  ],
  "buildPresets": [
    {
      "name": "default",
      "hidden": false,
      "configurePreset": "default", // Use the default configuration preset
      "jobs": 8, //<--- CHANGE THIS: Parallel build with N cores
      "targets": [
        "all" // Build the default target (usually all)
      ]
    },
    {
      "name": "debug-build",
      "configurePreset": "debug", // Use the debug configuration preset
      "jobs": 8, //<--- CHANGE THIS: Parallel build with N cores
      "targets": [
        "all"
      ]
    },
    {
      "name": "build-no-tests",
      "configurePreset": "no-tests", // Use the no-tests configuration preset
      "jobs": 8, //<--- CHANGE THIS: Parallel build with N cores
      "targets": [
        "all"
      ]
    },
    {
      "name": "clean",
      "configurePreset": "default",
      "jobs": 8, //<--- CHANGE THIS: Parallel build with N cores
      "targets": [
        "clean" // Run the clean target
      ]
    }
  ]
}
```

##### Comments on the Configuration

- **CMAKE_PREFIX_PATH**: Set the `CMAKE_PREFIX_PATH` and `CMAKE_INSTALL_PREFIX` to your `tudat-bundle` environment. This is critical to ensure CMake can find the necessary libraries and dependencies.

- **Generator**: Depending on your OS, the generator should be:
  - `"Unix Makefiles"` for **macOS** and **Linux**.
  - `"Ninja"` for **Windows** (if Ninja is installed).
  - `"Visual Studio"` for **Windows** if using Visual Studio.

- **Parallel Jobs**: The `jobs` key is used to specify the number of jobs for parallel builds. Set the number of cores to be used for the build.

- **Presets**:
  - `default`: Builds the project in Release mode with tests.
  - `debug`: Builds the project in Debug mode.
  - `no-tests`: Builds the project without tests.
  - `clean`: Cleans the build directory.

#### How to Use Presets in VSCode

1. **Open VSCode** in your project directory.

2. Open the Command Palette (`F1` or `Ctrl+Shift+P`) and run the command `CMake: Select Configure Preset`.

3. Choose one of the configure presets:
   - `default`: for Release builds.
   - `debug`: for Debug builds.
   - `no-tests`: to skip building tests.

4. Once configured, run the command `CMake: Select Build Preset` from the Command Palette, and choose one of the build presets:
   - `default`: to build the project.
   - `debug-build`: to build in Debug mode.
   - `build-no-tests`: to build without tests.
   - `clean`: to clean the build directory.

5. To start the build, run `CMake: Build` from the Command Palette or use the build button in the status bar.

#### Troubleshooting

- **CMake Errors**: If CMake cannot find dependencies, check that the `CMAKE_PREFIX_PATH` and `CMAKE_INSTALL_PREFIX` are set correctly in the `CMakePresets.json` file.
- **Wrong Generator**: If the build fails due to an unsupported generator, change the `"generator"` field in the preset to match your OS and toolchain.

<!--
The following line can also be edited if you wish to build tudatpy with its debug info (switching from `Release` to `RelWithDebInfo`; note that `Debug` is also available):
````
-DCMAKE_BUILD_TYPE=RelWithDebInfo
````

As building can take a while, you can build using multiple processors by appending by modifying your [build.sh](build.sh) script. For instance, you can modify the existing line defining the ``NUMBER_OF_PROCESSORS`` to the following, to use 2 threads for compilation. Note that a single thread may use up to 4 GB of RAM, and using too many parallel threads will make the compilation run out of RAM and terminate.
````
NUMBER_OF_PROCESSORS=${number_of_processors:-2}
````

6. Run the [build.sh](build.sh) script.

````
bash build.sh
```` -->

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

## C++ Code Formatting with Clang-Format
`tudat-bundle` uses `clang-format` to enforce a consistent code style. The configuration file is located in the root of the repository as `.clang-format`.

### Prerequisites

Ensure `clang-format` is installed. You can check the installation by running:

```bash
clang-format --version
```

If it's not installed, install `clang-format` using a package manager:

- **On Ubuntu/Debian**: `sudo apt install clang-format`
- **On macOS**: `brew install clang-format`
- **On Windows**: Download and install from [LLVM's official site](https://releases.llvm.org/).

### Configuration File

This `.clang-format` file configures formatting settings, including indentation width, brace wrapping, spacing rules, and alignment settings. These settings will ensure your code remains organized, readable, and consistent with `tudat` code style.

### Usage

#### 1. Running Clang-Format on a Single File

To format a single file, run:

```bash
clang-format -i path/to/your/file.cpp
```

The `-i` flag formats the file in place.

#### 2. Running Clang-Format on Multiple Files

To format all `.cpp` and `.h` files in your project directory:

```bash
find path/to/your/project -name '*.cpp' -o -name '*.h' | xargs clang-format -i
```

Alternatively, if you're using an IDE like Visual Studio Code, you can configure it to automatically format on save.

#### 3. Setting Up Automatic Formatting in VSCode

If you’re using Visual Studio Code, you can automate code formatting:

1. Install the **Clang-Format** extension from the marketplace.
2. Open your project’s **Settings** (File > Preferences > Settings).
3. Search for `Clang-Format` and set the path to your `.clang-format` file.
4. Enable **Format on Save** by setting `"editor.formatOnSave": true` in your `settings.json`.

#### 4. Checking Format Without Applying

To preview formatting changes without modifying files, run:

```bash
clang-format -output-replacements-xml path/to/your/file.cpp
```

This will show XML output if there are any format issues. You can then decide to apply formatting manually.

<!-- ## Use your build
The path of the TudatPy kernel that has been manually compiled needs to be added before importing any `tudatpy.kernel` module.
This can be done with the following two lines, with `<kernel_path>` being similar to `<tudat-bundle_path>/build/tudatpy`:
```
import sys
sys.path.insert(0, <kernel_path>)
``` -->

<!-- ## Notes

- [**All Users**] You can increase the number of cores used to compile `tudat` & `tudatpy` using the `-j<n>`
  build argument, but **be aware** that the current complexity of the libraries can often result in your PC freezing indefinitely. -->
