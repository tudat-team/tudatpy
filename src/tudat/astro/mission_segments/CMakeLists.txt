#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Set the source files.
set(mission_segments_SOURCES
        "escapeAndCapture.cpp"
        "gravityAssist.cpp"
        "improvedInversePolynomialWall.cpp"
        "lambertTargeter.cpp"
        "lambertTargeterIzzo.cpp"
        "lambertTargeterGooding.cpp"
        "lambertRoutines.cpp"
        "multiRevolutionLambertTargeterIzzo.cpp"
        "oscillatingFunctionNovak.cpp"
        "zeroRevolutionLambertTargeterIzzo.cpp"
        "transferNode.cpp"
        "transferLeg.cpp"
        "transferTrajectory.cpp"
        "createTransferTrajectory.cpp"
        )

# Set the header files.
set(mission_segments_HEADERS
        "escapeAndCapture.h"
        "gravityAssist.h"
        "improvedInversePolynomialWall.h"
        "lambertTargeter.h"
        "lambertTargeterIzzo.h"
        "lambertTargeterGooding.h"
        "lambertRoutines.h"
        "multiRevolutionLambertTargeterIzzo.h"
        "oscillatingFunctionNovak.h"
        "zeroRevolutionLambertTargeterIzzo.h"
        "transferNode.h"
        "transferLeg.h"
        "transferTrajectory.h"
        "createTransferTrajectory.h"
        )

# Add library.
TUDAT_ADD_LIBRARY("mission_segments"
        "${mission_segments_SOURCES}"
        "${mission_segments_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
