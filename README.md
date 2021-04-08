# tudat-bundle

This repository facillitates parallel development between the `tudat` (C++) and the
`tudatpy` (Python) library.

## Prerequisites

- [**Windows Users**] Windows Subsystem for Linux ([WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10))
  - All procedures, including the following prerequisite, assume the use of WSL. Power users who wish to do otherwise,
    must do so at their own risk, with reduced support from the team.
- Anaconda/Miniconda installation ([Installing Anaconda](https://tudat-space.readthedocs.io/en/latest/_src_first_steps/tudat_py.html#installing-anaconda))

## Setup

1. Clone the repository 

````
git clone https://github.com/tudat-team/tudat-bundle.git
````

2. (Optional) Switch to your desired branch

````
git checkout <branch>
````

3. Install the contained `environment.yaml` file to satisfy dependencies

````
conda env create -f environment.yaml
````

4. Activate the environment installed in step 1

````
conda activate tudat-bundle
````

5. Determine your `CONDA_PREFIX` path

````
echo $CONDA_PREFIX
````

6. Set the following CMake build configuration

````
-DCMAKE_PREFIX_PATH=<CONDA_PREFIX>
-DCMAKE_CXX_STANDARD=14
-DBoost_NO_BOOST_CMAKE=ON
````

## Notes

- [**CLion Users**] In CLion, the convention to set CMake arguments
  is to add them to `File>Settings>Build, Execution, Deployment>CMake Options`
