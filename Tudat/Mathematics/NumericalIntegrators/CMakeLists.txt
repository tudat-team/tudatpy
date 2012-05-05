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
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/burdenAndFairesNumericalIntegratorTest.h"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/matlabNumericalIntegratorTest.h"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/numericalIntegratorTestFunctionSuite.h"
)

# Add static libraries.
add_library(tudat_numerical_integrators STATIC ${NUMERICALINTEGRATORS_SOURCES} ${NUMERICALINTEGRATORS_HEADERS})
setup_tudat_library_target(tudat_numerical_integrators "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")

# Add unit tests.
add_executable(test_EulerIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestEulerIntegrator.cpp")
setup_custom_test_program(test_EulerIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_EulerIntegrator tudat_numerical_integrators ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_RungeKuttaVariableStepSizeIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaVariableStepSizeIntegrator.cpp")
setup_custom_test_program(test_RungeKuttaVariableStepSizeIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaVariableStepSizeIntegrator tudat_numerical_integrators ${Boost_LIBRARIES})

add_executable(test_RungeKuttaFehlberg45Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaFehlberg45Integrator.cpp")
setup_custom_test_program(test_RungeKuttaFehlberg45Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaFehlberg45Integrator tudat_numerical_integrators tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_RungeKuttaFehlberg56Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaFehlberg56Integrator.cpp")
setup_custom_test_program(test_RungeKuttaFehlberg56Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaFehlberg56Integrator tudat_numerical_integrators tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_RungeKuttaFehlberg78Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaFehlberg78Integrator.cpp")
setup_custom_test_program(test_RungeKuttaFehlberg78Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaFehlberg78Integrator tudat_numerical_integrators tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_RungeKutta87DormandPrinceIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKutta87DormandPrinceIntegrator.cpp")
setup_custom_test_program(test_RungeKutta87DormandPrinceIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKutta87DormandPrinceIntegrator tudat_numerical_integrators tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_RungeKuttaCoefficients "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaCoefficients.cpp")
setup_custom_test_program(test_RungeKuttaCoefficients "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaCoefficients tudat_numerical_integrators ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
