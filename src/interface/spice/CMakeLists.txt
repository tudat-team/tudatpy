#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

#Set the source files.
set(spice_interface_SOURCES
        "spiceEphemeris.cpp"
        "spiceRotationalEphemeris.cpp"
        "spiceInterface.cpp"
        )

# Set the header files.
set(spice_interface_HEADERS
        "spiceEphemeris.h"
        "spiceRotationalEphemeris.h"
        "spiceInterface.h"
        )

# Add library.
TUDAT_ADD_LIBRARY("spice_interface"
        "${spice_interface_SOURCES}"
        "${spice_interface_HEADERS}"
        PRIVATE_LINKS "${Boost_LIBRARIES}" "${CSpice_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" "${CSpice_INCLUDE_DIRS}"
        )

# Add tests.
TUDAT_ADD_TEST_CASE(SpiceInterface
        PRIVATE_LINKS
        tudat_ephemerides
        tudat_basic_mathematics
        tudat_spice_interface
        tudat_basic_astrodynamics
        )
