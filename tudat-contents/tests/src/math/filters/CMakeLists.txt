#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

TUDAT_ADD_TEST_CASE(LinearKalmanFilter PRIVATE_LINKS
        tudat_filters
        tudat_numerical_integrators
        tudat_statistics
        tudat_basics
        tudat_basic_mathematics
        tudat_input_output)

TUDAT_ADD_TEST_CASE(ExtendedKalmanFilter PRIVATE_LINKS
        tudat_filters
        tudat_numerical_integrators
        tudat_statistics
        tudat_basics
        tudat_basic_mathematics
        tudat_input_output)

TUDAT_ADD_TEST_CASE(UnscentedKalmanFilter PRIVATE_LINKS
        tudat_filters
        tudat_numerical_integrators
        tudat_statistics
        tudat_basics
        tudat_basic_mathematics
        tudat_input_output)
