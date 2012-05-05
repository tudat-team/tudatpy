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
 #	120203    K. Kumar	    Adapted for new BasicAstrodynamics package.
 #      120216    K. Kumar          Removed propagator files.
 #
 #    References
 #

# Set the source files.
set(BASICASTRODYNAMICS_SOURCES
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/keplerPropagator.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/convertMeanAnomalyToEccentricAnomaly.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/convertMeanAnomalyToHyperbolicEccentricAnomaly.cpp"
)

# Set the header files.
set(BASICASTRODYNAMICS_HEADERS 
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/keplerPropagator.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/convertMeanAnomalyToEccentricAnomaly.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/convertMeanAnomalyToHyperbolicEccentricAnomaly.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/pureMomentModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/momentDueToForceModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/forceModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/generalizedForceModel.h"
)

# Add static libraries.
add_library(tudat_basic_astrodynamics STATIC ${BASICASTRODYNAMICS_SOURCES} ${BASICASTRODYNAMICS_HEADERS})
setup_tudat_library_target(tudat_basic_astrodynamics "${SRCROOT}${BASICASTRODYNAMICSDIR}")

# Add unit tests.
add_executable(test_MeanToEccenticAnomalyConversion "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestConvertMeanAnomalyToEccentricAnomaly.cpp")
setup_custom_test_program(test_MeanToEccenticAnomalyConversion "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_MeanToEccenticAnomalyConversion tudat_basic_astrodynamics tudat_root_finding_methods ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_KeplerPropagator "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestKeplerPropagator.cpp")
setup_custom_test_program(test_KeplerPropagator "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_KeplerPropagator tudat_input_output tudat_bodies tudat_gravitation tudat_basic_astrodynamics tudat_states tudat_root_finding_methods ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_MomentDueToForce "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestMomentDueToForce.cpp")
setup_custom_test_program(test_MomentDueToForce "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_MomentDueToForce tudat_basic_astrodynamics)
