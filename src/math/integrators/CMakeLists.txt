#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Add source files.
set(numerical_integrators_SOURCES
        "rungeKuttaCoefficients.cpp"
        "createNumericalIntegrator.cpp"
        "bulirschStoerVariableStepsizeIntegrator.cpp"
        "tests/numericalIntegratorTests.cpp"
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
        "tests/burdenAndFairesNumericalIntegratorTest.h"
        "tests/numericalIntegratorTests.h"
        "tests/numericalIntegratorTestFunctions.h"
        )

# Add library.
TUDAT_ADD_LIBRARY("numerical_integrators"
        "${numerical_integrators_SOURCES}"
        "${numerical_integrators_HEADERS}"
        )

# Add tests.
TUDAT_ADD_TESTCASE(AdamsBashforthMoultonIntegrator tudat_numerical_integrators tudat_input_output)
TUDAT_ADD_TESTCASE(EulerIntegrator tudat_numerical_integrators tudat_input_output)
TUDAT_ADD_TESTCASE(NumericalIntegrator tudat_numerical_integrators)
TUDAT_ADD_TESTCASE(RungeKutta4Integrator tudat_numerical_integrators tudat_input_output)
TUDAT_ADD_TESTCASE(RungeKuttaVariableStepSizeIntegrator tudat_numerical_integrators)
TUDAT_ADD_TESTCASE(RungeKuttaCoefficients tudat_numerical_integrators)
TUDAT_ADD_TESTCASE(RungeKuttaFehlberg45Integrator tudat_numerical_integrators tudat_input_output)
TUDAT_ADD_TESTCASE(RungeKuttaFehlberg78Integrator tudat_numerical_integrators tudat_input_output)
TUDAT_ADD_TESTCASE(RungeKutta87DormandPrinceIntegrator tudat_numerical_integrators tudat_input_output)
TUDAT_ADD_TESTCASE(BulirschStoerVariableStepSizeIntegrator tudat_numerical_integrators tudat_input_output)
