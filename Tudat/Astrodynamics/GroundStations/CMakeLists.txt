 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Set the source files.
set(GROUND_STATIONS_SOURCES
  "${SRCROOT}${GROUNDSTATIONSDIR}/groundStation.cpp"
  "${SRCROOT}${GROUNDSTATIONSDIR}/groundStationState.cpp"
  "${SRCROOT}${GROUNDSTATIONSDIR}/pointingAnglesCalculator.cpp"
)

# Set the header files.
set(GROUND_STATIONS_HEADERS
  "${SRCROOT}${GROUNDSTATIONSDIR}/groundStation.h"
  "${SRCROOT}${GROUNDSTATIONSDIR}/groundStationState.h"
  "${SRCROOT}${GROUNDSTATIONSDIR}/pointingAnglesCalculator.h"
)

# Add static libraries.
add_library(tudat_ground_stations STATIC ${GROUND_STATIONS_SOURCES} ${GROUND_STATIONS_HEADERS})
setup_tudat_library_target(tudat_ground_stations "${SRCROOT}{GROUNDSTATIONSDIR}")


add_executable(test_GroundStationState "${SRCROOT}${GROUNDSTATIONSDIR}/UnitTests/unitTestGroundStationState.cpp")
setup_custom_test_program(test_GroundStationState "${SRCROOT}${GROUNDSTATIONSDIR}")
target_link_libraries(test_GroundStationState ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_PointingAnglesCalculator "${SRCROOT}${GROUNDSTATIONSDIR}/UnitTests/unitTestPointingAnglesCalculator.cpp")
setup_custom_test_program(test_PointingAnglesCalculator "${SRCROOT}${GROUNDSTATIONSDIR}")
target_link_libraries(test_PointingAnglesCalculator ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
