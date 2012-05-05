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
set(AERODYNAMICS_SOURCES
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicCoefficientGenerator.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicMoment.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicForce.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamics.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/exponentialAtmosphere.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/hypersonicLocalInclinationAnalysis.cpp"
  "${SRCROOT}${AERODYNAMICSDIR}/tabulatedAtmosphere.cpp"
)

# Set the header files.
set(AERODYNAMICS_HEADERS 
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicCoefficientGenerator.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicCoefficientInterface.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicMoment.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicAcceleration.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicForce.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicRotationalAcceleration.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamics.h"
  "${SRCROOT}${AERODYNAMICSDIR}/exponentialAtmosphere.h"
  "${SRCROOT}${AERODYNAMICSDIR}/hypersonicLocalInclinationAnalysis.h"
  "${SRCROOT}${AERODYNAMICSDIR}/tabulatedAtmosphere.h"
  "${SRCROOT}${AERODYNAMICSDIR}/aerodynamicForce.h"
  "${SRCROOT}${AERODYNAMICSDIR}/atmosphereModel.h"
  "${SRCROOT}${AERODYNAMICSDIR}/standardAtmosphere.h"
)

# Add static libraries.
add_library(tudat_aerodynamics STATIC ${AERODYNAMICS_SOURCES} ${AERODYNAMICS_HEADERS})
setup_tudat_library_target(tudat_aerodynamics "${SRCROOT}{AERODYNAMICSDIR}")

# Add unit tests.
add_executable(test_AerodynamicMomentAndAerodynamicForce "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestAerodynamicMomentAndAerodynamicForce.cpp")
setup_custom_test_program(test_AerodynamicMomentAndAerodynamicForce "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_AerodynamicMomentAndAerodynamicForce tudat_aerodynamics ${Boost_LIBRARIES})

add_executable(test_AerodynamicsNamespace "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestAerodynamicsNamespace.cpp")
setup_custom_test_program(test_AerodynamicsNamespace "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_AerodynamicsNamespace tudat_aerodynamics)

add_executable(test_CoefficientGenerator "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestCoefficientGenerator.cpp")
setup_custom_test_program(test_CoefficientGenerator "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_CoefficientGenerator tudat_aerodynamics tudat_geometric_shapes tudat_basic_mathematics ${TUDAT_CORE_LIBRARIES})

add_executable(test_ExponentialAtmosphere "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestExponentialAtmosphere.cpp")
setup_custom_test_program(test_ExponentialAtmosphere "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_ExponentialAtmosphere tudat_aerodynamics ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_TabulatedAtmosphere "${SRCROOT}${AERODYNAMICSDIR}/UnitTests/unitTestTabulatedAtmosphere.cpp")
setup_custom_test_program(test_TabulatedAtmosphere "${SRCROOT}${AERODYNAMICSDIR}")
target_link_libraries(test_TabulatedAtmosphere tudat_aerodynamics tudat_interpolators tudat_basic_mathematics tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
