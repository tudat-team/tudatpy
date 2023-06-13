#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

TUDAT_ADD_TEST_CASE(AdamsBashforthMoultonIntegrator PRIVATE_LINKS
        tudat_numerical_integrators
        tudat_input_output)

TUDAT_ADD_TEST_CASE(EulerIntegrator PRIVATE_LINKS
        tudat_test_support
        tudat_numerical_integrators
        tudat_input_output)

TUDAT_ADD_TEST_CASE(NumericalIntegrator
        PRIVATE_LINKS tudat_test_support tudat_numerical_integrators)

TUDAT_ADD_TEST_CASE(RungeKutta4Integrator
        PRIVATE_LINKS tudat_test_support tudat_numerical_integrators tudat_input_output)

TUDAT_ADD_TEST_CASE(RungeKuttaVariableStepSizeIntegrator
        PRIVATE_LINKS tudat_test_support tudat_numerical_integrators)

TUDAT_ADD_TEST_CASE(RungeKuttaCoefficients
        PRIVATE_LINKS tudat_test_support tudat_numerical_integrators)

TUDAT_ADD_TEST_CASE(RungeKuttaFehlberg45Integrator
        PRIVATE_LINKS tudat_test_support tudat_numerical_integrators tudat_input_output)

TUDAT_ADD_TEST_CASE(RungeKuttaFehlberg78Integrator
        PRIVATE_LINKS tudat_test_support tudat_numerical_integrators tudat_input_output)

TUDAT_ADD_TEST_CASE(RungeKutta87DormandPrinceIntegrator
        PRIVATE_LINKS tudat_test_support tudat_numerical_integrators tudat_input_output)

TUDAT_ADD_TEST_CASE(BulirschStoerVariableStepSizeIntegrator
        PRIVATE_LINKS tudat_test_support tudat_numerical_integrators tudat_input_output)

TUDAT_ADD_TEST_CASE(PerBlockStepSizeControl
        PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})
