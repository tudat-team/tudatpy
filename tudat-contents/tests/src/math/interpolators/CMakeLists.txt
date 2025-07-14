#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

TUDAT_ADD_TEST_CASE(CubicSplineInterpolator
        PRIVATE_LINKS
        tudat_input_output
        tudat_interpolators
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(HermiteCubicSplineInterpolator
        PRIVATE_LINKS
        tudat_input_output
        tudat_interpolators
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(LinearInterpolator
        PRIVATE_LINKS
        tudat_input_output
        tudat_interpolators
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(MultiLinearInterpolator
        PRIVATE_LINKS
        tudat_input_output
        tudat_interpolators
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(LagrangeInterpolators
        PRIVATE_LINKS
        tudat_input_output
        tudat_interpolators
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(InterpolatorVectorConversion
        PRIVATE_LINKS
        tudat_spice_interface
        tudat_input_output
        tudat_interpolators
        tudat_basic_mathematics
        )

