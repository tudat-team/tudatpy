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
set(interpolators_SOURCES
        "cubicSplineInterpolator.cpp"
        "linearInterpolator.cpp"
        "lagrangeInterpolator.cpp"
        "interpolator.cpp"
        "multiLinearInterpolator.cpp"
        )

# Add header files.
set(interpolators_HEADERS
        "cubicSplineInterpolator.h"
        "hermiteCubicSplineInterpolator.h"
        "linearInterpolator.h"
        "lagrangeInterpolator.h"
        "interpolator.h"
        "lookupScheme.h"
        "multiDimensionalInterpolator.h"
        "oneDimensionalInterpolator.h"
        "multiLinearInterpolator.h"
        "piecewiseConstantInterpolator.h"
        "jumpDataLinearInterpolator.h"
        "createInterpolator.h"
        )

# Add library.
TUDAT_ADD_LIBRARY("interpolators"
        "${interpolators_SOURCES}"
        "${interpolators_HEADERS}"
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
