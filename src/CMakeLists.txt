#    Copyright (c) 2010-2018, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# TODO: [@DominicDirkx] Setup the CMakelist.txt such that all .so objects are compiled under tudatpy/ and NOT tudatpy/src.
# TODO: Setup the CMakelist.txt such that tudatPy is installed to the detected dist-packages of python library.
# TODO: Setup the tests using the pyTest framework, or an alternative that works with "make test".
FILE(COPY __init__.py DESTINATION ..)
PYTHON_ADD_MODULE(core core.cpp)

# Constant value definition tests.
FILE(COPY test/test_constants.py DESTINATION .)
add_test(NAME constants COMMAND ${PYTHON_EXECUTABLE} -m pytest -rsv test_constants.py)

## Simulation setup tests.
#FILE(COPY test/simulation_setup.py DESTINATION .)
#ADD_TEST(NAME SimulationSetup COMMAND ${PYTHON_EXECUTABLE}/simulation_setup.py)
