#    Copyright (c) 2010-2018, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#
PYTHON_ADD_MODULE(test_module test_module.cpp)
FILE(COPY __init__.py DESTINATION .)
FILE(COPY test_constants.py DESTINATION .)
FILE(COPY test_environment_setup.py DESTINATION .)
FILE(COPY test_to_python_conversion.py DESTINATION .)
