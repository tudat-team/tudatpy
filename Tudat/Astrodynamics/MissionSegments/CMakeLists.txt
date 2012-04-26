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
 #      110820    S.M. Persson      File created.
 #      111025    K. Kumar          Adapted file to work with Revision 194.
 #      111026    K. Kumar          Adapted file so all headers show in project tree in Qt Creator.
 #

# Set the source files.
set(MISSIONSEGMENTS_SOURCES
  "${SRCROOT}${MISSIONSEGMENTSDIR}/escapeAndCapture.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/gravityAssist.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeter.cpp"
)

# Set the header files.
set(MISSIONSEGMENTS_HEADERS 
  "${SRCROOT}${MISSIONSEGMENTSDIR}/capturePhase.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/deepSpaceManeuver.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/escapeAndCapture.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/escapePhase.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/gravityAssist.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeter.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/trajectoryDesignMethod.h"
)

# Add static libraries.
add_library(tudat_mission_segments STATIC ${MISSIONSEGMENTS_SOURCES} ${MISSIONSEGMENTS_HEADERS})
setup_tudat_library_target(tudat_mission_segments "${SRCROOT}${MISSIONSEGMENTSDIR}")
setup_headers("${MISSIONSEGMENTS_HEADERS}" "${MISSIONSEGMENTSDIR}")

# Add unit tests.
add_executable(test_DeepSpaceManeuver "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestDeepSpaceManeuver.cpp")
setup_custom_test_program(test_DeepSpaceManeuver "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_DeepSpaceManeuver tudat_states tudat_mission_segments ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_EscapeAndCapture "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestEscapeAndCapture.cpp")
setup_custom_test_program(test_EscapeAndCapture "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_EscapeAndCapture tudat_bodies tudat_mission_segments tudat_gravitation tudat_basic_astrodynamics tudat_root_finding_methods tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_GravityAssist "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestGravityAssist.cpp")
setup_custom_test_program(test_GravityAssist "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_GravityAssist tudat_bodies tudat_gravitation tudat_mission_segments tudat_basic_astrodynamics tudat_root_finding_methods tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_LambertTargeter "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestLambertTargeter.cpp")
setup_custom_test_program(test_LambertTargeter "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_LambertTargeter tudat_bodies tudat_gravitation tudat_mission_segments tudat_root_finding_methods tudat_basic_astrodynamics tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
