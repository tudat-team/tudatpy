# \file CMakeLists.txt
 #    This is the CMake file for the States sub-directory.
 #
 #    Path              : /Astrodynamics/States/
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

# Add source files.
set(STATES_SOURCES
  "${SRCROOT}${STATESDIR}/approximatePlanetPositions.cpp"
  "${SRCROOT}${STATESDIR}/approximatePlanetPositionsBase.cpp"
  "${SRCROOT}${STATESDIR}/approximatePlanetPositionsCircularCoplanar.cpp"
  "${SRCROOT}${STATESDIR}/convertMeanAnomalyToEccentricAnomaly.cpp"
  "${SRCROOT}${STATESDIR}/convertMeanAnomalyToHyperbolicEccentricAnomaly.cpp"
  "${SRCROOT}${STATESDIR}/orbitalElementConversions.cpp"

  "${SRCROOT}${STATESDIR}/approximatePlanetPositionsDataContainer.h"
  "${SRCROOT}${STATESDIR}/cartesianElements.h"
  "${SRCROOT}${STATESDIR}/cartesianPositionElements.h"
  "${SRCROOT}${STATESDIR}/cartesianVelocityElements.h"
  "${SRCROOT}${STATESDIR}/convertMeanAnomalyBase.h"
  "${SRCROOT}${STATESDIR}/ephemeris.h"
  "${SRCROOT}${STATESDIR}/keplerianElements.h"
  "${SRCROOT}${STATESDIR}/state.h"
)

# Add header files.
set(STATES_HEADERS 
  "${STATESDIR}/approximatePlanetPositions.h"
  "${STATESDIR}/approximatePlanetPositionsBase.h"
  "${STATESDIR}/approximatePlanetPositionsCircularCoplanar.h"
  "${STATESDIR}/convertMeanAnomalyToEccentricAnomaly.h"
  "${STATESDIR}/convertMeanAnomalyToHyperbolicEccentricAnomaly.h"
  "${STATESDIR}/orbitalElementConversions.h"
)

# Add static libraries.
add_library(tudat_states STATIC ${STATES_SOURCES})
setup_tudat_library_target(tudat_states "${SRCROOT}${STATESDIR}")
setup_headers("${STATES_HEADERS}" "${STATESDIR}")

# Add unit tests.
add_executable(test_ApproximatePlanetPositions "${SRCROOT}${STATESDIR}/unitTestApproximatePlanetPositions.cpp")
setup_custom_test_program(test_ApproximatePlanetPositions "${SRCROOT}${STATESDIR}")
target_link_libraries(test_ApproximatePlanetPositions tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics)

add_executable(test_CartesianElements "${SRCROOT}${STATESDIR}/unitTestCartesianElements.cpp")
setup_custom_test_program(test_CartesianElements "${SRCROOT}${STATESDIR}")
target_link_libraries(test_CartesianElements tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics)

add_executable(test_KeplerianElements "${SRCROOT}${STATESDIR}/unitTestKeplerianElements.cpp")
setup_custom_test_program(test_KeplerianElements "${SRCROOT}${STATESDIR}")
target_link_libraries(test_KeplerianElements tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics)

add_executable(test_OrbitalElementConversions "${SRCROOT}${STATESDIR}/unitTestOrbitalElementConversions.cpp")
setup_custom_test_program(test_OrbitalElementConversions "${SRCROOT}${STATESDIR}")
target_link_libraries(test_OrbitalElementConversions tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics)

# End of file.
