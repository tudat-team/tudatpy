#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

TUDAT_ADD_TEST_CASE(MathematicalConstants)

TUDAT_ADD_TEST_CASE(LinearAlgebra PRIVATE_LINKS tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(CoordinateConversions PRIVATE_LINKS tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(NearestNeighbourSearch PRIVATE_LINKS tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(NumericalDerivative PRIVATE_LINKS tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(LegendrePolynomials PRIVATE_LINKS tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(SphericalHarmonics PRIVATE_LINKS tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(RotationAboutArbitraryAxis PRIVATE_LINKS tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(RotationPartials PRIVATE_LINKS tudat_basic_mathematics tudat_reference_frames)

