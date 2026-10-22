#    Copyright (c) 2010-2026, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Set the source files.
set(close_encounters_SOURCES
        "bPlane.cpp"
        )

# Set the header files.
set(close_encounters_HEADERS
        "bPlane.h"
        )

# Add library.
TUDAT_ADD_LIBRARY("close_encounters"
        "${close_encounters_SOURCES}"
        "${close_encounters_HEADERS}"
        )
