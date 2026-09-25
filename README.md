![TU Delft Astrodynamics Toolbox (Tudat)](docs/_static/cover_black_small.png)

The **TU Delft Astrodynamics Toolbox (Tudat)** is a powerful set of libraries that support astrodynamics and space research.
It can be used for a wide variety of purposes, ranging from the simulation studies of reentry dynamics to the processing of real tracking data of interplanetary missions.
The core functionality of Tudat is implemented in C++ and exposed to Python in the ``tudatpy`` package.
For a comprehensive overview of functionality and example applications, see our [user guide](https://docs.tudat.space/en/latest/) and the [API documentation](https://py.api.tudat.space/en/latest/).

To get started with Tudat, we recommend the installation of the ``tudatpy`` conda package, described on the [Installation page](https://docs.tudat.space/en/latest/getting-started/installation.html) of our user guide.
After that, take a look at our [quickstart guide](https://docs.tudat.space/en/latest/getting-started/quickstart.html) to set up your first orbit simulation with Tudat.
For complete example applications see our [examples page](https://docs.tudat.space/en/latest/index-examples.html) or the corresponding [tudatpy-examples](https://github.com/tudat-team/tudatpy-examples) repository.

> [!TIP]
> If you run into any issues while setting up your Tudat installation or simulations, please do not hesitate to contact us through our [Github discussions](https://github.com/orgs/tudat-team/discussions?discussions_q=) forum.

## Table of Contents

- [Structure of the `tudatpy` Repository](#structure-of-the-tudatpy-repository)
- [Compile `tudatpy` from Source](#compile-tudatpy-from-source)
  - [Prerequisites](#prerequisites)
  - [Compilation Steps](#compilation-steps)
  - [Verify your Build](#verify-your-build)
- [Contributing](#contributing)
- [Citation](#citation)

## Structure of the `tudatpy` Repository

This repository contains both the source code and the binding code, together with the respective documentation and examples folders.
The `src/tudatpy` directory mirrors the `tudatpy` module structure, which holds the C++ binding code and Python source code, alongside the docstrings of all package functionality.
The C++ source code is stored in the `include/tudat` and `src/tudat` directories.
The `.rst` files to build the API documentation are located in `docs/tudatpy`.

## Compile `tudatpy` from Source

### Prerequisites

- ``conda``: You must have ``conda`` installed on your system to obtain all required dependencies. See our [user guide](https://docs.tudat.space/en/latest/getting-started/use-of-tools/conda.html) for an introduction.
- **Windows Users**: While local builds on Windows are possible, note that support from the core developer team is limited.
  - Install [Visual Studio 2022, Version 17](https://learn.microsoft.com/en-us/visualstudio/releases/2022/release-history#evergreen-bootstrappers)
  - As an alternative, Windows Subsystem for Linux ([WSL](https://docs.microsoft.com/en-us/windows/wsl/install), see also our [user guide](https://docs.tudat.space/en/latest/getting-started/use-of-tools/windows-subsystem-for-linux.html)) can be installed for a Linux environment inside Windows.

### Compilation Steps

1. Clone the repository and enter directory

````
git clone https://github.com/tudat-team/tudatpy
cd tudatpy
````

2. Clone the `examples/tudatpy` submodule:

````
git submodule update --init --recursive
````

3. (optionally): Switch `tudatpy` to a new or an already existing branch using:

````
git checkout <branch-name>
````

We use the `develop` branch as default branch, which is updated frequently and used to build the weekly conda development version packages.
To build the latest stable version, use the `master` branch.

4. Retrieve dependencies

Create a conda environment from the provided `environment.yaml` file to retrieve all required dependencies, then activate it:

```
conda env create -f environment.yaml
conda activate tudatpy-dev
```

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
python build.py -h                    # Show help and available flags
python build.py -j <number-of-cores>  # Compile Tudatpy
python build.py --docs                # Compile Tudatpy and build API docs
```

This script compiles `tudatpy`. It will take some time to execute, but you can speed up the process by increasing the number of cores used with the `-j` flag.
Once the project is built, all the build output is dumped by default in a directory called `build`, which is not tracked by Git.
You can add the `--docs` flag to build the API documentation, which will be saved inside the build directory.

7. Install

The `install.py` script installs `tudatpy` in your active conda environment.
When using editable mode (appending the `-e` flag) the Python files in the environment are linked the tudatpy source using symbolic links.
File modifications therefore affect the installed package immediately, while the compiled kernel remains the one in the selected build directory.

> [!WARNING]
> When using editable mode the Python files in your environment track the current state of your `tudatpy` source repository.
> The Python environment is therefore not guaranteed to be stable.
> Use editable mode only for active development and keep the checkout on a compatible revision.

Run `python install.py` without `-e` when the environment must remain independent of later source edits or branch switches. This installs a snapshot of the current build through CMake. Unlike editable installations, regular installations are not tracked by `uninstall.py` and must currently be removed manually.

> [!NOTE]
> If you are using Windows, you might have to run the following commands from a shell with admin privileges, since it modifies files inside your conda environment.

```
python install.py -h                 # Show help and available flags
python install.py -e                 # Editable development installation
python install.py                    # Frozen installation of the current build
```


8. Uninstall

This script will remove tudatpy from your conda environment, but it will not delete the build directory.

```
python uninstall.py -h                # Show help and available flags
python uninstall.py                   # Uninstall Tudatpy
```

### Verify your Build

#### Running `tudatpy` tests

Within the `tudatpy` directory, run `pytest`  (packaged with CMake)

````
pytest
````
Use `pytest --no-remote-data` to skip tests that require a network connection.

Desired result:
````
=========================================== 6 passed in 1.78s ============================================
````

#### Running `tudat` tests

Note that `tudat` tests are only built when using the `--tests` flag with `build.py`, for example `python build.py --tests -j4`.

Enter the `tudatpy/build` directory and run the tests using `ctest`

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

## Contributing

We are open to and appreciate external contributions!
Please take a look at our [contribution guidelines](https://github.com/tudat-team/tudatpy?tab=contributing-ov-file) for more information.


## Citation

If you use Tudat in your research, we appreciate a citation!
While we are preparing a journal publication, please use the following paper to cite the project in general:

```bibtex
@inproceedings{dirkx2022OpensourceAstrodynamicsTudatpy,
  title = {The Open-Source Astrodynamics {{Tudatpy}} Software – Overview for Planetary Mission Design and Science Analysis},
  booktitle = {{{EPSC}}},
  author = {Dirkx, Dominic and Fayolle, Marie and Garrett, Geoffrey and Avillez, Miguel and Cowan, Kevin and Cowan, Sean and Encarnacao, Joao and Fortuny Lombrana, Carlos and Gaffarel, Jérémie and Hener, Jonas and Hu, Xuanyu and Van Nistelrooij, Maarten and Oggionni, Filippo and Plumaris, Michael},
  date = {2022},
  doi = {10.5194/epsc2022-253},
  url = {https://meetingorganizer.copernicus.org/EPSC2022/EPSC2022-253.html},
}
```

To reference the orbit estimation capabilities of Tudat, use the following paper:

```bibtex
@inproceedings{gisolfi2025OpenSourceHighFidelityOrbit,
  title = {Open-{{Source High-Fidelity Orbit Estimation}} for {{Planetary Science}} and {{Space Situational Awareness Using}} the {{Tudat Software}}},
  booktitle = {{{IAF Astrodynamics Symposium}}},
  author = {Gisolfi, Luigi and Dirkx, Dominic and Avillez, Miguel and Dijkstra, Tristan and Fayolle, Sam and Filice, Valerio and Hener, Jonas and Minervino Amodio, Andrea and Sanchez Rodriguez, Alfonso and López Rivera, Antonio and Dahmani, Fabien and Garrett, Geoffrey and Cowan, Kevin and Hinüber, Lars and Kimon Plumaris, Michael and Reichel, Markus and Van Hulle, Simon and Alkahal, Riva and Cimò, Giuseppe and Encarnacao, Joao and Jeanjean, Marceau and Langbroek, Marco and Søndergaard, Martin and Molera Calves, Guifre and Root, Bart and Gehly, Steve and Hu, Xuanyu and Witte, Daan and Verdoes Kleijn, Gijs and Williams, Rees and Stiller, Dominik},
  date = {2025},
  pages = {1136--1155},
  publisher = {International Astronautical Federation (IAF)},
  location = {Sydney, Australia},
  doi = {10.52202/083087-0100},
  url = {http://www.proceedings.com/083087-0100.html},
  eventtitle = {{{IAC}}},
  isbn = {979-8-3313-2935-8}
}
```
