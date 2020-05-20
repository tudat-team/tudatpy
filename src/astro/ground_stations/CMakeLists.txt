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
set(ground_stations_SOURCES
        "groundStation.cpp"
        "groundStationState.cpp"
        "pointingAnglesCalculator.cpp"
        )

# Set the header files.
set(ground_stations_HEADERS
        "groundStation.h"
        "groundStationState.h"
        "pointingAnglesCalculator.h"
        )

TUDAT_ADD_LIBRARY("ground_stations"
        "${ground_stations_SOURCES}"
        "${ground_stations_HEADERS}"
        )

# Add static libraries.
#add_library(tudat_ground_stations STATIC ${GROUND_STATIONS_SOURCES} ${GROUND_STATIONS_HEADERS})
#setup_tudat_library_target(tudat_ground_stations "${SRCROOT}{GROUNDSTATIONSDIR}")
#
#
#add_executable(test_GroundStationState "${SRCROOT}${GROUNDSTATIONSDIR}/tests/unitTestGroundStationState.cpp")
#setup_custom_test_program(test_GroundStationState "${SRCROOT}${GROUNDSTATIONSDIR}")
#target_link_libraries(test_GroundStationState ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_PointingAnglesCalculator "${SRCROOT}${GROUNDSTATIONSDIR}/tests/unitTestPointingAnglesCalculator.cpp")
#setup_custom_test_program(test_PointingAnglesCalculator "${SRCROOT}${GROUNDSTATIONSDIR}")
#target_link_libraries(test_PointingAnglesCalculator ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
