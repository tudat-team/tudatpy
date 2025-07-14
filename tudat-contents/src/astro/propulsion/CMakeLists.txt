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
set(propulsion_SOURCES
        "thrustAccelerationModel.cpp"
        "thrustGuidance.cpp"
        "thrustFunctions.cpp"
        "costateBasedThrustGuidance.cpp"
        "thrustMagnitudeWrapper.cpp"
        )

# Set the header files.
set(propulsion_HEADERS
        "thrustAccelerationModel.h"
        "thrustGuidance.h"
        "thrustFunctions.h"
        "thrustMagnitudeWrapper.h"
        "massRateFromThrust.h"
        "costateBasedThrustGuidance.h"
        )

TUDAT_ADD_LIBRARY("propulsion"
        "${propulsion_SOURCES}"
        "${propulsion_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
