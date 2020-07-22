#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Add header files.
set(root_finders_HEADERS
        "halleyRootFinder.h"
        "newtonRaphson.h"
        "rootFinder.h"
        "secantRootFinder.h"
        "bisection.h"
        "terminationConditions.h"
        "createRootFinder.h"
        )

# Add source files.
set(root_finders_SOURCES
        "createRootFinder.cpp"
        )

# Add library.
TUDAT_ADD_LIBRARY("root_finders"
        "${root_finders_SOURCES}"
        "${root_finders_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
