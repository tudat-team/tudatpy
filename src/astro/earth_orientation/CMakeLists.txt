#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(earth_orientation_SOURCES
        "earthOrientationCalculator.cpp"
        "terrestrialTimeScaleConverter.cpp"
        "eopReader.cpp"
        "polarMotionCalculator.cpp"
        "precessionNutationCalculator.cpp"
        "readAmplitudeAndArgumentMultipliers.cpp"
#        "tests/sofaEarthOrientationCookbookExamples.cpp"
        "shortPeriodEarthOrientationCorrectionCalculator.cpp"
        )

# Set the header files.
set(earth_orientation_HEADERS
        "earthOrientationCalculator.h"
        "terrestrialTimeScaleConverter.h"
        "eopReader.h"
        "polarMotionCalculator.h"
        "precessionNutationCalculator.h"
        "readAmplitudeAndArgumentMultipliers.h"
        "shortPeriodEarthOrientationCorrectionCalculator.h"
#        "tests/sofaEarthOrientationCookbookExamples.h"
        )

TUDAT_ADD_LIBRARY("earth_orientation"
        "${earth_orientation_SOURCES}"
        "${earth_orientation_HEADERS}"
        PRIVATE_LINKS "${Sofa_LIBRARIES}" "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" "${Sofa_INCLUDE_DIRS}"
        )
