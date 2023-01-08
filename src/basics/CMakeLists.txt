#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Add source files.
set(basics_SOURCES
        "utilities.cpp"
        "deprecationWarnings.cpp"
        )

# Add header files.
set(basics_HEADERS
        "utilities.h"
        "testMacros.h"
        "utilityMacros.h"
        "timeType.h"
        "basicTypedefs.h"
        "identityElements.h"
        "tudatTypeTraits.h"
        "deprecationWarnings.h"
        )

# Add library.
TUDAT_ADD_LIBRARY("basics"
        "${basics_SOURCES}"
        "${basics_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
