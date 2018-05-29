 #    Copyright (c) 2010-2012 Delft University of Technology.
 #
 #    This software is protected by national and international copyright.
 #    Any unauthorized use, reproduction or modification is unlawful and
 #    will be prosecuted. Commercial and non-private application of the
 #    software in any form is strictly prohibited unless otherwise granted
 #    by the authors.
 #
 #    The code is provided without any warranty; without even the implied
 #    warranty of merchantibility or fitness for a particular purpose.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      120607    P. Musegaas       Creation of file. 
 #

# Set the global macros for setting up targets. Do not modify!
# Different sub-directories are created for applications and unit tests.
macro(setup_application_executable_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}/applications")
endmacro(setup_application_executable_target)

macro(setup_library_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
  install(TARGETS ${target_name} LIBRARY DESTINATION "${BINROOT}" ARCHIVE DESTINATION "${LIBROOT}")
endmacro(setup_library_target)

macro(setup_unit_test_executable_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/unit_tests")
  get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
  add_test("${target_name}" "${BINROOT}/unit_tests/${target_name}")
endmacro(setup_unit_test_executable_target)

# Set the source files.
set(SPACETRAJECTORIES_SOURCES
  "${SRCROOT}${TRAJECTORYDIR}/captureLeg.cpp"
  "${SRCROOT}${TRAJECTORYDIR}/departureLegMga.cpp"
  "${SRCROOT}${TRAJECTORYDIR}/departureLegMga1DsmPosition.cpp"
  "${SRCROOT}${TRAJECTORYDIR}/departureLegMga1DsmVelocity.cpp"
  "${SRCROOT}${TRAJECTORYDIR}/exportTrajectory.cpp"
  "${SRCROOT}${TRAJECTORYDIR}/planetTrajectory.cpp"
  "${SRCROOT}${TRAJECTORYDIR}/swingbyLegMga.cpp"
  "${SRCROOT}${TRAJECTORYDIR}/swingbyLegMga1DsmPosition.cpp"
  "${SRCROOT}${TRAJECTORYDIR}/swingbyLegMga1DsmVelocity.cpp"
  "${SRCROOT}${TRAJECTORYDIR}/trajectory.cpp"
)

# Set the header files.
set(SPACETRAJECTORIES_HEADERS 
  "${SRCROOT}${TRAJECTORYDIR}/captureLeg.h"
  "${SRCROOT}${TRAJECTORYDIR}/departureLeg.h"
  "${SRCROOT}${TRAJECTORYDIR}/departureLegMga.h"
  "${SRCROOT}${TRAJECTORYDIR}/departureLegMga1DsmPosition.h"
  "${SRCROOT}${TRAJECTORYDIR}/departureLegMga1DsmVelocity.h"
  "${SRCROOT}${TRAJECTORYDIR}/exportTrajectory.h"
  "${SRCROOT}${TRAJECTORYDIR}/missionLeg.h"
  "${SRCROOT}${TRAJECTORYDIR}/planetTrajectory.h"
  "${SRCROOT}${TRAJECTORYDIR}/spaceLeg.h"
  "${SRCROOT}${TRAJECTORYDIR}/swingbyLeg.h"
  "${SRCROOT}${TRAJECTORYDIR}/swingbyLegMga.h"
  "${SRCROOT}${TRAJECTORYDIR}/swingbyLegMga1DsmPosition.h"
  "${SRCROOT}${TRAJECTORYDIR}/swingbyLegMga1DsmVelocity.h"
  "${SRCROOT}${TRAJECTORYDIR}/trajectory.h"
)

# Add static libraries, second line only if to be used later on outside this application.
add_library(tudat_trajectory_design STATIC ${SPACETRAJECTORIES_SOURCES} ${SPACETRAJECTORIES_HEADERS})
setup_tudat_library_target(tudat_trajectory_design "${SRCROOT}${TRAJECTORYDIR}")

# Add unit tests.
add_executable(test_CaptureLeg "${SRCROOT}${TRAJECTORYDIR}/UnitTests/unitTestCaptureLeg.cpp")
setup_unit_test_executable_target(test_CaptureLeg "${SRCROOT}${TRAJECTORYDIR}")
target_link_libraries(test_CaptureLeg tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})

# Add unit tests.
add_executable(test_DepartureLegMga "${SRCROOT}${TRAJECTORYDIR}/UnitTests/unitTestDepartureLegMga.cpp")
setup_unit_test_executable_target(test_DepartureLegMga "${SRCROOT}${TRAJECTORYDIR}")
target_link_libraries(test_DepartureLegMga tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})

# Add unit tests.
add_executable(test_DepartureLegMga1DsmPosition "${SRCROOT}${TRAJECTORYDIR}/UnitTests/unitTestDepartureLegMga1DsmPosition.cpp")
setup_unit_test_executable_target(test_DepartureLegMga1DsmPosition "${SRCROOT}${TRAJECTORYDIR}")
target_link_libraries(test_DepartureLegMga1DsmPosition tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})

# Add unit tests.
add_executable(test_DepartureLegMga1DsmVelocity "${SRCROOT}${TRAJECTORYDIR}/UnitTests/unitTestDepartureLegMga1DsmVelocity.cpp")
setup_unit_test_executable_target(test_DepartureLegMga1DsmVelocity "${SRCROOT}${TRAJECTORYDIR}")
target_link_libraries(test_DepartureLegMga1DsmVelocity tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})

# Add unit tests.
add_executable(test_SwingbyLegMga "${SRCROOT}${TRAJECTORYDIR}/UnitTests/unitTestSwingbyLegMga.cpp")
setup_unit_test_executable_target(test_SwingbyLegMga "${SRCROOT}${TRAJECTORYDIR}")
target_link_libraries(test_SwingbyLegMga tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})

# Add unit tests.
add_executable(test_SwingbyLegMga1DsmPosition "${SRCROOT}${TRAJECTORYDIR}/UnitTests/unitTestSwingbyLegMga1DsmPosition.cpp")
setup_unit_test_executable_target(test_SwingbyLegMga1DsmPosition "${SRCROOT}${TRAJECTORYDIR}")
target_link_libraries(test_SwingbyLegMga1DsmPosition tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})

# Add unit tests.
add_executable(test_SwingbyLegMga1DsmVelocity "${SRCROOT}${TRAJECTORYDIR}/UnitTests/unitTestSwingbyLegMga1DsmVelocity.cpp")
setup_unit_test_executable_target(test_SwingbyLegMga1DsmVelocity "${SRCROOT}${TRAJECTORYDIR}")
target_link_libraries(test_SwingbyLegMga1DsmVelocity tudat_trajectory_design tudat_mission_segments tudat_basic_mathematics ${Boost_LIBRARIES})

# Add unit tests.
add_executable(test_Trajectory "${SRCROOT}${TRAJECTORYDIR}/UnitTests/unitTestTrajectory.cpp")
setup_unit_test_executable_target(test_Trajectory "${SRCROOT}${TRAJECTORYDIR}")
target_link_libraries(test_Trajectory tudat_trajectory_design tudat_mission_segments tudat_ephemerides tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
