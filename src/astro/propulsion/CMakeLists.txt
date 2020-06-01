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
        "thrustGuidance.cpp"
        "thrustFunctions.cpp"
        "costateBasedThrustGuidance.cpp"
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
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )


# Add static libraries.
#add_library(tudat_propulsion STATIC ${PROPULSION_SOURCES} ${PROPULSION_HEADERS})
#setup_tudat_library_target(tudat_propulsion "${SRCROOT}${PROPULSIONDIR}")
#
#add_executable(test_ThrustAcceleration "${SRCROOT}${PROPULSIONDIR}/tests/unitTestThrustAcceleration.cpp")
#setup_custom_test_program(test_ThrustAcceleration "${SRCROOT}${PROPULSIONDIR}")
#target_link_libraries(test_ThrustAcceleration ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
