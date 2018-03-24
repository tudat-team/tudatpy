/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#define BOOST_TEST_MAIN

#include <vector>
#include <limits>

#include <Eigen/Core>

#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/Statistics/randomSampling.h"
#include "Tudat/Mathematics/Statistics/basicStatistics.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_Random_Sampler )

BOOST_AUTO_TEST_CASE( test_randomVectorUniform )
{
    int dimension = 3;
    int numberOfSamples = 1E6;
    int seed = 511;

    Eigen::VectorXd lower( dimension );
    Eigen::VectorXd upper( dimension );
    lower << 0.0, 1.0, -2.0;
    upper << 1.0, 3.0, 4.0;

    Eigen::VectorXd width = upper - lower;
    Eigen::VectorXd average = (upper + lower) / 2.0;

    Eigen::VectorXd standardDeviation = std::sqrt( 1.0 / 12.0 ) * width;
    {
        std::vector< Eigen::VectorXd > samples =
                tudat::statistics::generateUniformRandomSample( seed, numberOfSamples, lower, upper );

        // Compute sample mean and standard deviation
        Eigen::VectorXd sampleMean = statistics::computeSampleMean( samples );
        Eigen::ArrayXd sampleStandardDeviations = statistics::computeSampleVariance( samples ).array( ).sqrt( );

        BOOST_CHECK_SMALL( std::fabs( average( 0 ) - sampleMean( 0 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( average( 1 ) - sampleMean( 1 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( average( 2 ) - sampleMean( 2 ) ), 5.0E-3 );

        BOOST_CHECK_SMALL( std::fabs( standardDeviation( 0 ) - sampleStandardDeviations( 0 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( standardDeviation( 1 ) - sampleStandardDeviations( 1 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( standardDeviation( 2 ) - sampleStandardDeviations( 2 ) ), 5.0E-3 );
    }

    {
        std::vector< Eigen::VectorXd > samples = tudat::statistics::generateUniformRandomSample( seed, numberOfSamples, dimension );

        // Compute sample mean and standard deviation
        Eigen::VectorXd sampleMean = statistics::computeSampleMean( samples );
        Eigen::ArrayXd sampleStandardDeviations = statistics::computeSampleVariance( samples ).array( ).sqrt( );

        BOOST_CHECK_SMALL( std::fabs( 0.5 - sampleMean( 0 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( 0.5 - sampleMean( 1 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( 0.5 - sampleMean( 2 ) ), 5.0E-3 );

        BOOST_CHECK_SMALL( std::fabs( std::sqrt( 1.0 / 12.0 ) - sampleStandardDeviations( 0 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( std::sqrt( 1.0 / 12.0 ) - sampleStandardDeviations( 1 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( std::sqrt( 1.0 / 12.0 ) - sampleStandardDeviations( 2 ) ), 5.0E-3 );
    }
}

BOOST_AUTO_TEST_CASE( test_randomVectorGaussian )
{
    int dimension = 3;
    int numberOfSamples = 1E6;
    int seed = 511;

    Eigen::VectorXd mean( dimension );
    Eigen::VectorXd standardDeviation( dimension );
    mean << 0.0, 1.0, -2.0;
    standardDeviation << 1.0, 3.0, 4.0;

    {
        std::vector< Eigen::VectorXd > samples =
                tudat::statistics::generateGaussianRandomSample( seed, numberOfSamples, mean, standardDeviation );

        // Compute sample mean and standard deviation
        Eigen::VectorXd sampleMean = statistics::computeSampleMean( samples );
        Eigen::ArrayXd sampleStandardDeviations = statistics::computeSampleVariance( samples ).array( ).sqrt( );

        BOOST_CHECK_SMALL( std::fabs( mean( 0 ) - sampleMean( 0 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( mean( 1 ) - sampleMean( 1 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( mean( 2 ) - sampleMean( 2 ) ), 5.0E-3 );

        BOOST_CHECK_SMALL( std::fabs( standardDeviation( 0 ) - sampleStandardDeviations( 0 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( standardDeviation( 1 ) - sampleStandardDeviations( 1 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( standardDeviation( 2 ) - sampleStandardDeviations( 2 ) ), 5.0E-3 );
    }

    {
        std::vector< Eigen::VectorXd > samples =
                tudat::statistics::generateGaussianRandomSample( seed, numberOfSamples, 3 );

        // Compute sample mean and standard deviation
        Eigen::VectorXd sampleMean = statistics::computeSampleMean( samples );
        Eigen::ArrayXd sampleStandardDeviations = statistics::computeSampleVariance( samples ).array( ).sqrt( );

        BOOST_CHECK_SMALL( std::fabs( 0.0 - sampleMean( 0 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( 0.0 - sampleMean( 1 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( 0.0 - sampleMean( 2 ) ), 5.0E-3 );

        BOOST_CHECK_SMALL( std::fabs( 1.0 - sampleStandardDeviations( 0 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( 1.0 - sampleStandardDeviations( 1 ) ), 5.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( 1.0 - sampleStandardDeviations( 2 ) ), 5.0E-3 );
    }
}


#if USE_GSL

//! Test if Sobol sampler interface is working correctly. Note that this test is somewhat minimal, but the core of the
//! Sobol sampling is tested in the GSL unit tests.
BOOST_AUTO_TEST_CASE( test_Sobol_Sampler )
{
    // Define settings and generate samples.
    int dimension = 3;
    int numberOfSamples = 1E6;

    Eigen::VectorXd lower( dimension );
    Eigen::VectorXd upper( dimension );
    lower << 0.0, 1.0, -2.0;
    upper << 1.0, 3.0, 4.0;

    Eigen::VectorXd width = upper - lower;
    Eigen::VectorXd average = (upper + lower) / 2.0;

    std::vector< Eigen::VectorXd > sobolSamples = tudat::statistics::generateVectorSobolSample(
                numberOfSamples, lower, upper );

    // Compute sample average
    Eigen::VectorXd sampleMean = statistics::computeSampleMean( sobolSamples );

    BOOST_CHECK_SMALL( std::fabs( average( 0 ) - sampleMean( 0 ) ), 2.0E-6 );
    BOOST_CHECK_SMALL( std::fabs( average( 1 ) - sampleMean( 1 ) ), 2.0E-6 );
    BOOST_CHECK_SMALL( std::fabs( average( 2 ) - sampleMean( 2 ) ), 2.0E-6 );

    sobolSamples = tudat::statistics::generateVectorSobolSample( dimension, numberOfSamples );

    // Compute sample average
    sampleMean = statistics::computeSampleMean( sobolSamples );

    BOOST_CHECK_SMALL( std::fabs( 0.5 - sampleMean( 0 ) ), 2.0E-6 );
    BOOST_CHECK_SMALL( std::fabs( 0.5 - sampleMean( 1 ) ), 2.0E-6 );
    BOOST_CHECK_SMALL( std::fabs( 0.5 - sampleMean( 2 ) ), 2.0E-6 );
}
#endif

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
