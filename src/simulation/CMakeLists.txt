#    Copyright (c) 2010-2018, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

PYTHON_ADD_MODULE(simulation_setup SimulationSetup.cpp)
FILE(COPY simulation_setup.py DESTINATION .)
ADD_TEST(NAME src COMMAND ${PYTHON_EXECUTABLE} simulation_setup.py)
