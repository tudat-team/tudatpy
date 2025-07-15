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
conda activate tudat-bundle
````

> **Note**\
It is possible that the creation of the environment will 'time out'. A likely reason for this is that the packages required cannot be found by the current channel, `conda-forge`. It is then advisable to add the channel `anaconda` to ensure a proper creation of the environment.
>

5. Build TudatPy

```
python build.py -h                   # Show help and available flags
python build.py -j <number-of-cores>  # Compile Tudatpy
```
This script compiles Tudatpy. It will take some time to execute, but you can speed up the process by increasing the number of cores used with the `-j` flag.

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
`Tudatpy` uses `clang-format` to enforce a consistent code style. The configuration file is located in the root of the repository as `.clang-format`.

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
