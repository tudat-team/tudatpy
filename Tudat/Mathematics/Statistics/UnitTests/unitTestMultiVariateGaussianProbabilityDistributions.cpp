/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <iostream>
#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include <boost/math/distributions/lognormal.hpp>

#include "Tudat/Mathematics/Statistics/multiVariateGaussianProbabilityDistributions.h"

namespace tudat
{
namespace unit_tests
{

double gaussianCopulaProbabilityDensity( const Eigen::VectorXd& independentVariables , int dimension_ , Eigen::MatrixXd correlationMatrix_ )
{
    Eigen::MatrixXd inverseCorrelationMatrix_ = correlationMatrix_.inverse() ;
    double determinant_ = correlationMatrix_.determinant() ;

    double probabilityDensity = 0.0 ;
    Eigen::VectorXd y( dimension_ ) ;
    boost::math::normal distribution( 0.0 , 1.0 );

    for(int i = 0 ; i < dimension_ ; i++)
    {
        y(i) = boost::math::quantile( distribution , independentVariables(i) ); // Inverse cdf
    }

    // Calculate probability density
    Eigen::MatrixXd location = - 0.5 * ( y.transpose() *
                                         (inverseCorrelationMatrix_ - Eigen::MatrixXd::Identity( dimension_ , dimension_ ) ) * y ) ;

    probabilityDensity = ( ( 1.0/( std::sqrt( determinant_ ) ) )*std::exp( location(0,0) ) ) ;
    return probabilityDensity;
}


// Exact solution of 2D Gaussian distribution
// Montgomery, D. C. & Runger, G. C. Applied Statistics and Probability for engineers Wiley, 2014
double computeBiGaussianPdf(
        const double standardDeviationX, const double standardDeviationY, const double correlation, const double meanX, const double meanY, const Eigen::VectorXd independentVariables )
{
    using tudat::mathematical_constants::PI;
    using std::pow;
    using std::exp;

    return ( 1.0 / ( 2.0 * PI*standardDeviationX*standardDeviationY*sqrt( 1.0 - pow( correlation, 2.0 ) ) ) ) *
            exp( ( -1.0 / ( 2.0 * ( 1.0 - pow( correlation, 2.0) ) ) ) * (
                     ( pow( independentVariables( 0 ) - meanX, 2.0 ) ) /( pow( standardDeviationX, 2.0 ) )
                     - ( 2.0 * correlation * ( independentVariables( 0 ) - meanX ) * ( independentVariables( 1 )- meanY ) ) / ( standardDeviationX * standardDeviationY )
                     + ( pow( independentVariables( 1 ) - meanY, 2.0 ) ) / ( pow( standardDeviationY, 2.0) ) ) );
}


BOOST_AUTO_TEST_SUITE( test_probability_distributions )

//! Test if ND Gaussian distribution class works correctly.
BOOST_AUTO_TEST_CASE( testMultiDimensionalGaussianDistribution )
{
    using namespace tudat::statistics;


    Eigen::VectorXd mean(2);
    mean << 0.0 , 1.0 ;
    Eigen::MatrixXd covariance(2,2);
    covariance <<   3.0 , -1.0 ,
            -1.0 , 3.0 ;
    GaussianDistributionXd distribution( mean, covariance );

    Eigen::VectorXd independentVariables( 2 );
    independentVariables << 0.0 , 1.0 ;

    {
        double standardDeviationX = std::sqrt( covariance( 0, 0 ) );
        double standardDeviationY = std::sqrt( covariance( 1, 1 ) );
        double meanX = mean( 0 );
        double meanY = mean( 1 );
        double correlation = covariance( 0 , 1 ) / ( standardDeviationX * standardDeviationY );

        double exactSolution = computeBiGaussianPdf( standardDeviationX, standardDeviationY, correlation, meanX, meanY, independentVariables );


        BOOST_CHECK_SMALL( std::fabs( distribution.evaluatePdf(independentVariables) - exactSolution ) ,
                           std::numeric_limits<double>::epsilon( ) );
    }

    // Second test
    mean << -1.0 , 2.0 ;
    covariance <<   4.0 , 1.5 ,
            1.5 , 4.0 ;
    independentVariables << 2.5 , 0.0 ;

    {
        double standardDeviationX = std::sqrt( covariance( 0, 0 ) );
        double standardDeviationY = std::sqrt( covariance( 1, 1 ) );
        double meanX = mean( 0 );
        double meanY = mean( 1 );
        double correlation = covariance( 0 , 1 ) / ( standardDeviationX * standardDeviationY );

        GaussianDistributionXd distribution2(mean,covariance);

        double exactSolution = computeBiGaussianPdf( standardDeviationX, standardDeviationY, correlation, meanX, meanY, independentVariables );

        BOOST_CHECK_SMALL( std::fabs( distribution2.evaluatePdf(independentVariables) - exactSolution ) ,
                           std::numeric_limits<double>::epsilon() ) ;
    }
}

//! Test if Gaussian Copula distribution class works correctly.
BOOST_AUTO_TEST_CASE( testGaussianCopula )
{
    using namespace tudat::statistics;
    using tudat::mathematical_constants::PI;

    int dimension = 2 ;
    Eigen::MatrixXd correlationMatrix( dimension , dimension ) ;
    correlationMatrix << 1.0 , 0.3 ,
            0.3 , 1.0 ;
    GaussianCopulaDistributionXd distribution( dimension , correlationMatrix );

    Eigen::VectorXd independentVariables(2);
    independentVariables << 0.5 , 0.3 ;

    BOOST_CHECK_SMALL( std::fabs( gaussianCopulaProbabilityDensity( independentVariables , dimension , correlationMatrix )
                                  - distribution.evaluatePdf( independentVariables ) ), 1E-15 );

    // Test out of bounds
    independentVariables << 1.1 , 0.5 ;
    BOOST_CHECK_SMALL( std::fabs( distribution.evaluatePdf( independentVariables ) - 0.0 ) , 1E-15 );

    independentVariables << 0.1 , -0.5 ;
    BOOST_CHECK_SMALL( std::fabs( distribution.evaluatePdf( independentVariables ) - 0.0 ) , 1E-15 );

    independentVariables << -0.1 , -0.5 ;
    BOOST_CHECK_SMALL( std::fabs( distribution.evaluatePdf( independentVariables ) - 0.0 ) , 1E-15 );

    independentVariables << 0.1 , 1.5 ;
    BOOST_CHECK_SMALL( std::fabs( distribution.evaluatePdf( independentVariables ) - 0.0 ) , 1E-15 );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
