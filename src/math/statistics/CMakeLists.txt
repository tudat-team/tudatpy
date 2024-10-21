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
set(statistics_SOURCES
        "basicStatistics.cpp"
        "simpleLinearRegression.cpp"
        "multiVariateGaussianProbabilityDistributions.cpp"
        "continuousProbabilityDistributions.cpp"
        "boostProbabilityDistributions.cpp"
        "kernelDensityDistribution.cpp"
        "randomSampling.cpp"
        "randomVariableGenerator.cpp"
        "allanVariance.cpp"
        )

# Add header files.
set(statistics_HEADERS
        "basicStatistics.h"
        "simpleLinearRegression.h"
        "multiVariateGaussianProbabilityDistributions.h"
        "continuousProbabilityDistributions.h"
        "boostProbabilityDistributions.h"
        "kernelDensityDistribution.h"
        "randomSampling.h"
        "randomVariableGenerator.h"
        "allanVariance.h"
        )

# Add library.
if(TUDAT_BUILD_WITH_FFTW3)

    set(statistics_HEADERS "${statistics_HEADERS};fastFourierTransform.h;powerLawNoiseGeneration.h")
    set(statistics_SOURCES "${statistics_SOURCES};fastFourierTransform.cpp;powerLawNoiseGeneration.cpp")

    TUDAT_ADD_LIBRARY("statistics"
            "${statistics_SOURCES}"
            "${statistics_HEADERS}"
            "${FFTW3_LIBRARIES}"
            )
else()
    TUDAT_ADD_LIBRARY("statistics"
            "${statistics_SOURCES}"
            "${statistics_HEADERS}"
            "${FFTW3_LIBRARIES}"
            )
endif()
