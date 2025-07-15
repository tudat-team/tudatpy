#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

TUDAT_ADD_TEST_CASE(SimpleLinearRegression PRIVATE_LINKS
        tudat_statistics
        tudat_basics
        )

TUDAT_ADD_TEST_CASE(BasicStatistics PRIVATE_LINKS
        tudat_statistics
        tudat_basics
        )

TUDAT_ADD_TEST_CASE(MultiVariateGaussianProbabilityDistributions PRIVATE_LINKS
        tudat_statistics
        tudat_basics
        )

TUDAT_ADD_TEST_CASE(BoostRandomVariables PRIVATE_LINKS
        tudat_statistics
        tudat_basics
        )

TUDAT_ADD_TEST_CASE(KernelDensityDistribution PRIVATE_LINKS
        tudat_statistics
        tudat_basics
        )

TUDAT_ADD_TEST_CASE(RandomSampling PRIVATE_LINKS
        tudat_statistics
        tudat_basics
        )

if(TUDAT_BUILD_WITH_FFTW3)
    TUDAT_ADD_TEST_CASE(PowerLawNoise PRIVATE_LINKS
            tudat_statistics
            tudat_basic_mathematics
            tudat_basics
            ${FFTW3_LIBRARIES}
            )

    TUDAT_ADD_TEST_CASE(Fft PRIVATE_LINKS
            tudat_statistics tudat_basic_mathematics ${FFTW3_LIBRARIES})
endif()