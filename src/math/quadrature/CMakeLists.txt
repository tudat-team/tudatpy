#   Copyright (c) 2010-2019, Delft University of Technology
#   All rigths reserved
#
#   This file is part of the Tudat. Redistribution and use in source and
#   binary forms, with or without modification, are permitted exclusively
#   under the terms of the Modified BSD license. You should have received
#   a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Add header files.
set(numerical_quadrature_HEADERS
        "numericalQuadrature.h"
        "trapezoidQuadrature.h"
        "gaussianQuadrature.h"
        "createNumericalQuadrature.h"
        )

# Add source files.
set(numerical_quadrature_SOURCES
        "dummyNumericalQuadrature.cpp"
        "gaussianQuadrature.cpp"
        "createNumericalQuadrature.cpp"
        )

# Add library.
TUDAT_ADD_LIBRARY("numerical_quadrature"
        "${numerical_quadrature_SOURCES}"
        "${numerical_quadrature_HEADERS}"
        )

# Add tests.
TUDAT_ADD_TESTCASE(TrapezoidQuadrature tudat_numerical_quadrature)
TUDAT_ADD_TESTCASE(GaussianQuadrature tudat_input_output tudat_numerical_quadrature)
