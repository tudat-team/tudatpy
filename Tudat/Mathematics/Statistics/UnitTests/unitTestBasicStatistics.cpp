/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110207    B. Romgens        File created.
 *      110215    K. Kumar          Minor modifications to layout, comments
 *                                  and variable-naming.
 *      110411    K. Kumar          Added unit test for
 *                                  convertCartesianToSpherical( ) function.
 *      110701    K. Kumar          Updated failing tests with relative errors.
 *      110708    K. Kumar          Added unit tests for computeSampleMean( )
 *                                  and computeSampleVariance( ) functions.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      111111    K. Kumar          Strange error with convertCylindricalToCartesian function;
 *                                  achieved precision of results is less than machine precision,
 *                                  fixed by using slightly larger precision tolerance.
 *      120202    K. Kumar          Moved unit tests from unitTestBasicMathematicsFunctions.h/.cpp;
 *                                  rewrote unit tests using Boost unit test framework.
 *
 *    References
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
    double computedSampleMean = tudat::mathematics::statistics::computeSampleMean( sampleData );

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
    double computedSampleVariance = tudat::mathematics::statistics::
            computeSampleVariance( sampleData );

    // Check if computed sample variance matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedSampleVariance, expectedSampleVariance,
                                std::numeric_limits< double >::epsilon( ) );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
