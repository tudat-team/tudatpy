 #    Copyright (c) 2010-2017, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Add source files.
set(NUMERICALINTEGRATORS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/rungeKuttaCoefficients.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/bulirschStoerVariableStepsizeIntegrator.cpp"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/numericalIntegratorTests.cpp"
)

# Add header files.
set(NUMERICALINTEGRATORS_HEADERS 
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/rungeKuttaCoefficients.h"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/burdenAndFairesNumericalIntegratorTest.h"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/bulirschStoerVariableStepsizeIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/euler.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/numericalIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/reinitializableNumericalIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/rungeKutta4Integrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/createNumericalIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/numericalIntegratorTests.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/numericalIntegratorTestFunctions.h"
)

# Add static libraries.
add_library(tudat_numerical_integrators STATIC ${NUMERICALINTEGRATORS_SOURCES} ${NUMERICALINTEGRATORS_HEADERS})
setup_tudat_library_target(tudat_numerical_integrators "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")

# Add unit tests.
add_executable(test_EulerIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestEulerIntegrator.cpp")
setup_custom_test_program(test_EulerIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_EulerIntegrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_NumericalIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestNumericalIntegrator.cpp")
setup_custom_test_program(test_NumericalIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_NumericalIntegrator tudat_numerical_integrators ${Boost_LIBRARIES})

add_executable(test_RungeKutta4Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKutta4Integrator.cpp")
setup_custom_test_program(test_RungeKutta4Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKutta4Integrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_RungeKuttaVariableStepSizeIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaVariableStepSizeIntegrator.cpp")
setup_custom_test_program(test_RungeKuttaVariableStepSizeIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaVariableStepSizeIntegrator tudat_numerical_integrators ${Boost_LIBRARIES})

add_executable(test_RungeKuttaCoefficients "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaCoefficients.cpp")
setup_custom_test_program(test_RungeKuttaCoefficients "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaCoefficients tudat_numerical_integrators ${Boost_LIBRARIES})

add_executable(test_RungeKuttaFehlberg45Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaFehlberg45Integrator.cpp")
setup_custom_test_program(test_RungeKuttaFehlberg45Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaFehlberg45Integrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_RungeKuttaFehlberg56Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaFehlberg56Integrator.cpp")
setup_custom_test_program(test_RungeKuttaFehlberg56Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaFehlberg56Integrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_RungeKuttaFehlberg78Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaFehlberg78Integrator.cpp")
setup_custom_test_program(test_RungeKuttaFehlberg78Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaFehlberg78Integrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_RungeKutta87DormandPrinceIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKutta87DormandPrinceIntegrator.cpp")
setup_custom_test_program(test_RungeKutta87DormandPrinceIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKutta87DormandPrinceIntegrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_BulirschStoerVariableStepSizeIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestBulirschStoerVariableStepSizeIntegrator.cpp")
 setup_custom_test_program(test_BulirschStoerVariableStepSizeIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_BulirschStoerVariableStepSizeIntegrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

