***************************************
# TudatPy
***************************************

TU Delft Astrodynamics Toolbox in Python, or **TudatPy**, is a library that primarily exposes the powerful set of C++ 
libraries, [Tudat](https://tudat.tudelft.nl/). TudatPy aims at accelerating the implementation of Tudat simulations,
providing and interface between Tudat and popular machine learning frameworks and establishing a platform to provide 
quality education in the field of astrodynamics. See the [documentation](https://ggarrett13.github.io/tudatpy/) for more.

## Installation

1. Follow my [forked tudatBundle installation](https://github.com/ggarrett13/tudatBundle) replacing the first step with

        git clone https://github.com/tudat/tudatBundle.git
        
2. After successfully building the bundle, enter the build directory for tudatpy after step 6 (build/tudatpy)

        cd tudatpy
        
3. Ensure the following Python dependencies are met

        sudo apt-get install python3-dev &&
        sudo apt install python3-pip &&
        pip3 install numpy
        
4. Install the compiled Python bindings to your Python library using the ``setup.py`` script

        pip install .