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
set(trajectory_design_SOURCES
        "captureLeg.cpp"
        "departureLegMga.cpp"
        "departureLegMga1DsmPosition.cpp"
        "departureLegMga1DsmVelocity.cpp"
        "exportTrajectory.cpp"
        "planetTrajectory.cpp"
        "swingbyLegMga.cpp"
        "swingbyLegMga1DsmPosition.cpp"
        "swingbyLegMga1DsmVelocity.cpp"
        "trajectory.cpp"
        )

# Set the header files.
set(trajectory_design_HEADERS
        "captureLeg.h"
        "departureLeg.h"
        "departureLegMga.h"
        "departureLegMga1DsmPosition.h"
        "departureLegMga1DsmVelocity.h"
        "exportTrajectory.h"
        "missionLeg.h"
        "planetTrajectory.h"
        "spaceLeg.h"
        "swingbyLeg.h"
        "swingbyLegMga.h"
        "swingbyLegMga1DsmPosition.h"
        "swingbyLegMga1DsmVelocity.h"
        "trajectory.h"
        )

TUDAT_ADD_LIBRARY("trajectory_design"
        "${trajectory_design_SOURCES}"
        "${trajectory_design_HEADERS}"
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
