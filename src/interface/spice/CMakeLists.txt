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
        PRIVATE_LINKS "${CSpice_LIBRARIES}"
        )
