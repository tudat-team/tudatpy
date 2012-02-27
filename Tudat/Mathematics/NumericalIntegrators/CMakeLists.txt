 #    Copyright (c) 2012 Delft University of Technology.
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
 #      120127    B. Tong Minh      File created.
 #

# Add source files.
set(NUMERICALINTEGRATORS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/rungeKuttaCoefficients.cpp"
)

# Add header files.
set(NUMERICALINTEGRATORS_HEADERS 
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/euler.h"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/rungeKuttaCoefficients.h"
)

# Add static libraries.
add_library(tudat_numerical_integrators STATIC ${NUMERICALINTEGRATORS_SOURCES} ${NUMERICALINTEGRATORS_HEADERS})
setup_tudat_library_target(tudat_numerical_integrators "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
setup_headers("${NUMERICALINTEGRATORS_HEADERS}" "${MATHEMATICSDIR}/NumericalIntegrators")

# Add unit tests.
add_executable(test_EulerIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestEulerIntegrator.cpp")
setup_custom_test_program(test_EulerIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_EulerIntegrator tudat_numerical_integrators ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_RungeKuttaVariableStepsizeIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaVariableStepSizeIntegrator.cpp")
setup_custom_test_program(test_RungeKuttaVariableStepsizeIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaVariableStepsizeIntegrator tudat_numerical_integrators ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
