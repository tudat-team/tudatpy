#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

SET(CMAKE_INSTALL_RPATH "/home/ggarrett/Repositories/new/tudat_bundle/install_test/lib")

# Add source files.
set(numerical_integrators_SOURCES
        "rungeKuttaCoefficients.cpp"
        "createNumericalIntegrator.cpp"
        "bulirschStoerVariableStepsizeIntegrator.cpp"
#        "tests/numericalIntegratorTests.cpp"
        "numericalIntegrator.cpp"
        "reinitializableNumericalIntegrator.cpp"
        "rungeKutta4Integrator.cpp"
        "rungeKuttaVariableStepSizeIntegrator.cpp"
        "adamsBashforthMoultonIntegrator.cpp"
        "bulirschStoerVariableStepsizeIntegrator.cpp"
        )

# Add header files.
set(numerical_integrators_HEADERS
        "adamsBashforthMoultonIntegrator.h"
        "createNumericalIntegrator.h"
        "bulirschStoerVariableStepsizeIntegrator.h"
        "euler.h"
        "numericalIntegrator.h"
        "reinitializableNumericalIntegrator.h"
        "rungeKutta4Integrator.h"
        "rungeKuttaCoefficients.h"
        "rungeKuttaVariableStepSizeIntegrator.h"
#        "tests/burdenAndFairesNumericalIntegratorTest.h"
#        "tests/numericalIntegratorTests.h"
#        "tests/numericalIntegratorTestFunctions.h"
        )

# Add library.
TUDAT_ADD_LIBRARY("numerical_integrators"
        "${numerical_integrators_SOURCES}"
        "${numerical_integrators_HEADERS}"
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

# Add tests.
TUDAT_ADD_TEST_CASE(AdamsBashforthMoultonIntegrator
        PRIVATE_LINKS tudat_numerical_integrators tudat_input_output
        )

TUDAT_ADD_TEST_CASE(EulerIntegrator
        PRIVATE_LINKS tudat_numerical_integrators tudat_input_output
        SOURCES "tests/numericalIntegratorTests.cpp"
        )

TUDAT_ADD_TEST_CASE(NumericalIntegrator
        PRIVATE_LINKS tudat_numerical_integrators
        SOURCES "tests/numericalIntegratorTests.cpp"
        )

TUDAT_ADD_TEST_CASE(RungeKutta4Integrator
        PRIVATE_LINKS tudat_numerical_integrators tudat_input_output
        SOURCES "tests/numericalIntegratorTests.cpp"
        )

TUDAT_ADD_TEST_CASE(RungeKuttaVariableStepSizeIntegrator
        PRIVATE_LINKS tudat_numerical_integrators
        SOURCES "tests/numericalIntegratorTests.cpp"
        )

TUDAT_ADD_TEST_CASE(RungeKuttaCoefficients
        PRIVATE_LINKS tudat_numerical_integrators
        SOURCES "tests/numericalIntegratorTests.cpp"
        )

TUDAT_ADD_TEST_CASE(RungeKuttaFehlberg45Integrator
        PRIVATE_LINKS tudat_numerical_integrators tudat_input_output
        SOURCES "tests/numericalIntegratorTests.cpp"
        )

TUDAT_ADD_TEST_CASE(RungeKuttaFehlberg78Integrator
        PRIVATE_LINKS tudat_numerical_integrators tudat_input_output
        SOURCES "tests/numericalIntegratorTests.cpp"
        )

TUDAT_ADD_TEST_CASE(RungeKutta87DormandPrinceIntegrator
        PRIVATE_LINKS tudat_numerical_integrators tudat_input_output
        SOURCES "tests/numericalIntegratorTests.cpp"
        )

TUDAT_ADD_TEST_CASE(BulirschStoerVariableStepSizeIntegrator
        PRIVATE_LINKS tudat_numerical_integrators tudat_input_output
        SOURCES "tests/numericalIntegratorTests.cpp"
        )
