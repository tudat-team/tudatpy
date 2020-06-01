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
        )

# Add library.
TUDAT_ADD_LIBRARY("statistics"
        "${statistics_SOURCES}"
        "${statistics_HEADERS}"
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

# Add test library interface.
add_library(statistics_test_library INTERFACE)
target_link_libraries(statistics_test_library
        INTERFACE
        tudat_statistics
        tudat_basics
        )

# Add tests.
TUDAT_ADD_TEST_CASE(SimpleLinearRegression PRIVATE_LINKS statistics_test_library)
TUDAT_ADD_TEST_CASE(BasicStatistics PRIVATE_LINKS statistics_test_library)
TUDAT_ADD_TEST_CASE(MultiVariateGaussianProbabilityDistributions PRIVATE_LINKS statistics_test_library)
TUDAT_ADD_TEST_CASE(BoostRandomVariables PRIVATE_LINKS statistics_test_library)
TUDAT_ADD_TEST_CASE(KernelDensityDistribution PRIVATE_LINKS statistics_test_library)
TUDAT_ADD_TEST_CASE(RandomSampling PRIVATE_LINKS statistics_test_library)
