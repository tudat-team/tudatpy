 #    Copyright (c) 2010-2018, Delft University of Technology
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
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/rungeKuttaCoefficients.cpp"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/bulirschStoerVariableStepsizeIntegrator.cpp"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/numericalIntegratorTests.cpp"
)

# Add header files.
set(NUMERICALINTEGRATORS_HEADERS 
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/adamsBashforthMoultonIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/createNumericalIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/bulirschStoerVariableStepsizeIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/euler.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/numericalIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/reinitializableNumericalIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/rungeKutta4Integrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/rungeKuttaCoefficients.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/rungeKuttaVariableStepSizeIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/burdenAndFairesNumericalIntegratorTest.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/numericalIntegratorTests.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/numericalIntegratorTestFunctions.h"
)

# Add static libraries.
add_library(tudat_numerical_integrators STATIC ${NUMERICALINTEGRATORS_SOURCES} ${NUMERICALINTEGRATORS_HEADERS})
setup_tudat_library_target(tudat_numerical_integrators "${SRCROOT}${NUMERICALINTEGRATORSDIR}")

# Add unit tests.

add_executable(test_AdamsBashforthMoultonIntegrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/unitTestAdamsBashforthMoultonIntegrator.cpp")
setup_custom_test_program(test_AdamsBashforthMoultonIntegrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}")
target_link_libraries(test_AdamsBashforthMoultonIntegrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_EulerIntegrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/unitTestEulerIntegrator.cpp")
setup_custom_test_program(test_EulerIntegrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}")
target_link_libraries(test_EulerIntegrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_NumericalIntegrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/unitTestNumericalIntegrator.cpp")
setup_custom_test_program(test_NumericalIntegrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}")
target_link_libraries(test_NumericalIntegrator tudat_numerical_integrators ${Boost_LIBRARIES})

add_executable(test_RungeKutta4Integrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/unitTestRungeKutta4Integrator.cpp")
setup_custom_test_program(test_RungeKutta4Integrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}")
target_link_libraries(test_RungeKutta4Integrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_RungeKuttaVariableStepSizeIntegrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/unitTestRungeKuttaVariableStepSizeIntegrator.cpp")
setup_custom_test_program(test_RungeKuttaVariableStepSizeIntegrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}")
target_link_libraries(test_RungeKuttaVariableStepSizeIntegrator tudat_numerical_integrators ${Boost_LIBRARIES})

add_executable(test_RungeKuttaCoefficients "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/unitTestRungeKuttaCoefficients.cpp")
setup_custom_test_program(test_RungeKuttaCoefficients "${SRCROOT}${NUMERICALINTEGRATORSDIR}")
target_link_libraries(test_RungeKuttaCoefficients tudat_numerical_integrators ${Boost_LIBRARIES})

add_executable(test_RungeKuttaFehlberg45Integrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/unitTestRungeKuttaFehlberg45Integrator.cpp")
setup_custom_test_program(test_RungeKuttaFehlberg45Integrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}")
target_link_libraries(test_RungeKuttaFehlberg45Integrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_RungeKuttaFehlberg56Integrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/unitTestRungeKuttaFehlberg56Integrator.cpp")
setup_custom_test_program(test_RungeKuttaFehlberg56Integrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}")
target_link_libraries(test_RungeKuttaFehlberg56Integrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_RungeKuttaFehlberg78Integrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/unitTestRungeKuttaFehlberg78Integrator.cpp")
setup_custom_test_program(test_RungeKuttaFehlberg78Integrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}")
target_link_libraries(test_RungeKuttaFehlberg78Integrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_RungeKutta87DormandPrinceIntegrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/unitTestRungeKutta87DormandPrinceIntegrator.cpp")
setup_custom_test_program(test_RungeKutta87DormandPrinceIntegrator "${SRCROOT}${NUMERICALINTEGRATORSDIR}")
target_link_libraries(test_RungeKutta87DormandPrinceIntegrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_BulirschStoerVariableStepSizeIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestBulirschStoerVariableStepSizeIntegrator.cpp")
setup_custom_test_program(test_BulirschStoerVariableStepSizeIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_BulirschStoerVariableStepSizeIntegrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

