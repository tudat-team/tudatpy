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
        )


# Add static libraries, second line only if to be used later on outside this application.
#add_library(tudat_trajectory_design STATIC ${SPACETRAJECTORIES_SOURCES} ${SPACETRAJECTORIES_HEADERS})
#setup_tudat_library_target(tudat_trajectory_design "${SRCROOT}${TRAJECTORYDIR}")
#
## Add unit tests.
#add_executable(test_CaptureLeg "${SRCROOT}${TRAJECTORYDIR}/tests/unitTestCaptureLeg.cpp")
#setup_unit_test_executable_target(test_CaptureLeg "${SRCROOT}${TRAJECTORYDIR}")
#target_link_libraries(test_CaptureLeg tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})
#
## Add unit tests.
#add_executable(test_DepartureLegMga "${SRCROOT}${TRAJECTORYDIR}/tests/unitTestDepartureLegMga.cpp")
#setup_unit_test_executable_target(test_DepartureLegMga "${SRCROOT}${TRAJECTORYDIR}")
#target_link_libraries(test_DepartureLegMga tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})
#
## Add unit tests.
#add_executable(test_DepartureLegMga1DsmPosition "${SRCROOT}${TRAJECTORYDIR}/tests/unitTestDepartureLegMga1DsmPosition.cpp")
#setup_unit_test_executable_target(test_DepartureLegMga1DsmPosition "${SRCROOT}${TRAJECTORYDIR}")
#target_link_libraries(test_DepartureLegMga1DsmPosition tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})
#
## Add unit tests.
#add_executable(test_DepartureLegMga1DsmVelocity "${SRCROOT}${TRAJECTORYDIR}/tests/unitTestDepartureLegMga1DsmVelocity.cpp")
#setup_unit_test_executable_target(test_DepartureLegMga1DsmVelocity "${SRCROOT}${TRAJECTORYDIR}")
#target_link_libraries(test_DepartureLegMga1DsmVelocity tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})
#
## Add unit tests.
#add_executable(test_SwingbyLegMga "${SRCROOT}${TRAJECTORYDIR}/tests/unitTestSwingbyLegMga.cpp")
#setup_unit_test_executable_target(test_SwingbyLegMga "${SRCROOT}${TRAJECTORYDIR}")
#target_link_libraries(test_SwingbyLegMga tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})
#
## Add unit tests.
#add_executable(test_SwingbyLegMga1DsmPosition "${SRCROOT}${TRAJECTORYDIR}/tests/unitTestSwingbyLegMga1DsmPosition.cpp")
#setup_unit_test_executable_target(test_SwingbyLegMga1DsmPosition "${SRCROOT}${TRAJECTORYDIR}")
#target_link_libraries(test_SwingbyLegMga1DsmPosition tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})
#
## Add unit tests.
#add_executable(test_SwingbyLegMga1DsmVelocity "${SRCROOT}${TRAJECTORYDIR}/tests/unitTestSwingbyLegMga1DsmVelocity.cpp")
#setup_unit_test_executable_target(test_SwingbyLegMga1DsmVelocity "${SRCROOT}${TRAJECTORYDIR}")
#target_link_libraries(test_SwingbyLegMga1DsmVelocity tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})
#
## Add unit tests.
#add_executable(test_Trajectory "${SRCROOT}${TRAJECTORYDIR}/tests/unitTestTrajectory.cpp")
#setup_unit_test_executable_target(test_Trajectory "${SRCROOT}${TRAJECTORYDIR}")
#target_link_libraries(test_Trajectory tudat_trajectory_design tudat_mission_segments tudat_ephemerides tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
