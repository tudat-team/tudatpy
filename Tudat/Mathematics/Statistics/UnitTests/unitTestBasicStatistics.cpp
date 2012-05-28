/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
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
