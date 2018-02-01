/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/Statistics/basicStatistics.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_basic_statistics )

//! Test if sample mean is computed correctly.
BOOST_AUTO_TEST_CASE( testSampleMean )
{
    // Test computation of sample mean on finite population using unbiased estimators.
    // The expected values are computed using the Microsoft Excel the AVERAGE( ) function.

    // Declare vector of sample data.
    std::vector< double > sampleData;

    // Populate vector with sample data.
    sampleData.push_back( 2.5 );
    sampleData.push_back( 6.4 );
    sampleData.push_back( 8.9 );
    sampleData.push_back( 12.7 );
    sampleData.push_back( 15.0 );

    // Set expected sample mean.
    double expectedSampleMean = 9.1;

    // Compute sample mean.
    double computedSampleMean = statistics::computeSampleMean( sampleData );

    // Check if computed sample mean matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedSampleMean, expectedSampleMean,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test if sample variance is computed correctly.
BOOST_AUTO_TEST_CASE( testSampleVariance )
{
    // Test computation of sample variance on finite population using unbiased estimators.
    // The expected values are computed using the Microsoft Excel the VAR( ) function.

    // Declare vector of sample data.
    std::vector< double > sampleData;

    // Populate vector with sample data.
    sampleData.push_back( 2.5 );
    sampleData.push_back( 6.4 );
    sampleData.push_back( 8.9 );
    sampleData.push_back( 12.7 );
    sampleData.push_back( 15.0 );

    // Declare expected sample variance.
    double expectedSampleVariance = 24.665;

    // Compute sample variance.
    double computedSampleVariance = statistics::computeSampleVariance( sampleData );

    // Check if computed sample variance matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedSampleVariance, expectedSampleVariance,
                                std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
