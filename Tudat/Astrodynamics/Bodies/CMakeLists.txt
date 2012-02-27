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
 #      120209    B. Tong Minh      Merge Ephemeris package
 #

# Set Ephemeris path.
set(EPHEMERISDIR "${BODIESDIR}/Ephemeris")

# Set the source files.
set(BODIES_SOURCES
  "${SRCROOT}${BODIESDIR}/planet.cpp"
  "${SRCROOT}${EPHEMERISDIR}/approximatePlanetPositions.cpp"
  "${SRCROOT}${EPHEMERISDIR}/approximatePlanetPositionsBase.cpp"
  "${SRCROOT}${EPHEMERISDIR}/approximatePlanetPositionsCircularCoplanar.cpp"
)

# Set the header files.
set(BODIES_HEADERS 
  "${SRCROOT}${BODIESDIR}/body.h"
  "${SRCROOT}${BODIESDIR}/celestialBody.h"
  "${SRCROOT}${BODIESDIR}/planet.h"
  "${SRCROOT}${BODIESDIR}/vehicle.h"
  "${SRCROOT}${BODIESDIR}/vehicleExternalModel.h"
  "${SRCROOT}${EPHEMERISDIR}/approximatePlanetPositions.h"
  "${SRCROOT}${EPHEMERISDIR}/approximatePlanetPositionsBase.h"
  "${SRCROOT}${EPHEMERISDIR}/approximatePlanetPositionsCircularCoplanar.h"
  "${SRCROOT}${EPHEMERISDIR}/approximatePlanetPositionsDataContainer.h"
  "${SRCROOT}${EPHEMERISDIR}/ephemeris.h"
)

# Add static libraries.
add_library(tudat_bodies STATIC ${BODIES_SOURCES} ${BODIES_HEADERS})
setup_tudat_library_target(tudat_bodies "${SRCROOT}${BODIESDIR}")
setup_headers("${BODIES_HEADERS}" "${BODIESDIR}")

# Add unit tests.
add_executable(test_ApproximatePlanetPositions "${SRCROOT}${EPHEMERISDIR}/UnitTests/unitTestApproximatePlanetPositions.cpp")
setup_custom_test_program(test_ApproximatePlanetPositions "${SRCROOT}${EPHEMERISDIR}")
target_link_libraries(test_ApproximatePlanetPositions tudat_states tudat_bodies tudat_gravitation tudat_root_finding_methods tudat_basic_astrodynamics tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
