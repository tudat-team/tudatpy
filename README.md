
# TudatPy

TU Delft Astrodynamics Toolbox in Python, or **TudatPy**, is a library that primarily exposes the powerful set of C++ 
libraries, [Tudat](https://tudat.tudelft.nl/). TudatPy aims at accelerating the implementation of Tudat simulations,
providing and interface between Tudat and popular machine learning frameworks and establishing a platform to provide 
quality education in the field of astrodynamics. See the [documentation](https://ggarrett13.github.io/tudatpy/) for more.

## Installation

1. Install the Debian Python 3 development package (required for Python header inclusion for Pybind11).
   
        sudo apt-get install python3-dev
           
2. Install the package installer for Python 3 (required for final TudatPy installation step)

        sudo apt install python3-pip
           
3. Install the NumPy package (for Eigen matrix conversion to Python NumPy arrays)

        pip3 install numpy
           
**All together** (ignore if 1-3 followed)

    sudo apt-get install python3-dev &&
    sudo apt install python3-pip &&
    pip3 install numpy

### Clean TudatBundle Installation

4. Clone my [forked tudatBundle installation](https://github.com/ggarrett13/tudatBundle)

        git clone https://github.com/ggarrett13/tudatBundle.git
        
5. Enter the new directory

        cd tudatBundle
        
6. Checkout all the submodules

        git submodule update --init --recursive
        
7. Make a new build and directory and enter

         mkdir build && cd build
         
8. Initiate CMake for the project

        cmake ..
        
9. Build the project (If an error occurs, check below)

        make
        
10. Enter the TudatPy build folder

        cd build/tudatpy
        
11. Install the TudatPy shared library and Python files to your Python library

        pip install .
    
### Existing TudatBundle Installation

Since the project is in early development, I don't recommend this step. However if you choose to take this path,
I assume there's a reason for it and you know your way around cmake. (I don't guarantee this option being hiccup free)
 
1. Clone [Pybind11](https://github.com/pybind/pybind11) to the root directory of your tudatBundle.

        git clone https://github.com/pybind/pybind11

2. Clone [TudatPy](https://github.com/ggarrett13/tudatpy) to the root directory of your tudatBundle.

        git clone https://github.com/ggarrett13/tudatpy

3. Add the following lines to the tudatBundle ``CMakeList.txt`` to compile position independent code (near the top)

        set(CMAKE_POSITION_INDEPENDENT_CODE ON)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        
4. Add the following lines before the ``if(USE_PAGMO)`` in the tudatBundle ``CMakeList.txt`` (line 85 currently).

        set(pybind11_INCLUDE_DIR "${PROJECTROOT}pybind11/include")
        add_subdirectory("${PROJECTROOT}/pybind11/")
        list(APPEND BoostComponents python3)
        find_package(pybind11 2.0 REQUIRED)

5. Add the following line to the end of the tudatBundle ``CMakeList.txt``.

        add_subdirectory("${PROJECTROOT}/tudatpy")

6. Cross your fingers and follow steps 1-3 and 8-11.

## Possible Errors

### Missing Python header ``#include<Python.h>``.

Ensure steps 1-3 are followed with the Python interpreter found during the cmake configuration (step 8). The
configuration process will print out the Python ``interpreter``, ``interpreter version`` and ``include dir``. If the
detected Python installation corresponds to the one used to carry out steps 1-3, see 'Creating an issue on GitHub' below.

## Creating an issue on GitHub

post an issue on the [TudatPy](https://github.com/ggarrett13/tudatpy)
GitHub page using the issue template ``<tudatpy>/docs/issue_template.md``.