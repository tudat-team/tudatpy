# \file CMakeLists.txt
 #    This is the CMake file for the MissionSegments sub-directory.
 #
 #    Path              : /Astrodynamics/MissionSegments/
 #    Version           : 1
 #    Check status      : Checked
 #
 #    Author            : S.M. Persson
 #    Affiliation       : McGill University, Department of Mechanical Engineering
 #    E-mail address    : mikael.s.persson@gmail.com
 #
 #    Author            : K. Kumar
 #    Affiliation       : Delft University of Technology
 #    E-mail address    : K.Kumar@tudelft.nl
 #
 #    Author            : S.M. Persson
 #    Affiliation       : McGill University, Department of Mechanical Engineering
 #    E-mail address    : mikael.s.persson@gmail.com
 #
 #    Date created      : 20 August, 2011
 #    Last modified     : 25 October, 2011
 #
 #    References
 #
 #    Notes
 #
 #    Copyright (c) 2010-2011 Delft University of Technology.
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
 #/

# Set the source files.
set(MISSIONSEGMENTS_SOURCES
  "${SRCROOT}${MISSIONSEGMENTSDIR}/escapeAndCapture.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/gravityAssist.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeter.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/librationPoint.cpp"
)

# Set the header files.
set(MISSIONSEGMENTS_HEADERS 
  "${MISSIONSEGMENTSDIR}/capturePhase.h"
  "${MISSIONSEGMENTSDIR}/deepSpaceManeuver.h"
  "${MISSIONSEGMENTSDIR}/escapeAndCapture.h"
  "${MISSIONSEGMENTSDIR}/escapePhase.h"
  "${MISSIONSEGMENTSDIR}/gravityAssist.h"
  "${MISSIONSEGMENTSDIR}/lambertTargeter.h"
  "${MISSIONSEGMENTSDIR}/librationPoint.h"
  "${MISSIONSEGMENTSDIR}/trajectoryDesignMethod.h"
)

# Add static libraries.
add_library(tudat_mission_segments STATIC ${MISSIONSEGMENTS_SOURCES})
setup_tudat_library_target(tudat_mission_segments "${SRCROOT}${MISSIONSEGMENTSDIR}")
setup_headers("${MISSIONSEGMENTS_HEADERS}" "${MISSIONSEGMENTSDIR}")

# Add unit tests.
add_executable(test_DeepSpaceManeuver "${SRCROOT}${MISSIONSEGMENTSDIR}/unitTestDeepSpaceManeuver.cpp")
setup_custom_test_program(test_DeepSpaceManeuver "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_DeepSpaceManeuver tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics)

add_executable(test_EscapeAndCapture "${SRCROOT}${MISSIONSEGMENTSDIR}/unitTestEscapeAndCapture.cpp")
setup_custom_test_program(test_EscapeAndCapture "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_EscapeAndCapture tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics)

add_executable(test_GravityAssist "${SRCROOT}${MISSIONSEGMENTSDIR}/unitTestGravityAssist.cpp")
setup_custom_test_program(test_GravityAssist "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_GravityAssist tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics)

add_executable(test_LambertTargeter "${SRCROOT}${MISSIONSEGMENTSDIR}/unitTestLambertTargeter.cpp")
setup_custom_test_program(test_LambertTargeter "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_LambertTargeter tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics)

add_executable(test_LibrationPoints "${SRCROOT}${MISSIONSEGMENTSDIR}/unitTestLibrationPoints.cpp")
setup_custom_test_program(test_LibrationPoints "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_LibrationPoints tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics)

# End of file.
