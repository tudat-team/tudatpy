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

// Test function for Gaussian cupola distribution (computation of pdf).
// NOTE: This manner of testing the Gaussian cupola is not ideal, but no usable test data has been identified.
double gaussianCopulaProbabilityDensity(
        const Eigen::VectorXd& independentVariables , int dimension_ , Eigen::MatrixXd correlationMatrix_ )
{
    Eigen::MatrixXd inverseCorrelationMatrix_ = correlationMatrix_.inverse() ;
    double determinant_ = correlationMatrix_.determinant( );

    double probabilityDensity = 0.0 ;
    Eigen::VectorXd y( dimension_ ) ;
    boost::math::normal distribution( 0.0 , 1.0 );

    for( int i = 0 ; i < dimension_; i++ )
    {
        y( i ) = boost::math::quantile( distribution , independentVariables( i ) ); // Inverse cdf
    }

    // Calculate probability density
    Eigen::MatrixXd location = -0.5 *
            ( y.transpose( )  * ( inverseCorrelationMatrix_ - Eigen::MatrixXd::Identity( dimension_ , dimension_ ) ) * y ) ;

    probabilityDensity = ( ( 1.0 / ( std::sqrt( determinant_ ) ) ) * std::exp( location( 0, 0 ) ) ) ;
    return probabilityDensity;
}


// Exact solution of 2D Gaussian distribution
// Montgomery, D. C. & Runger, G. C. Applied Statistics and Probability for engineers Wiley, 2014
double computeBiGaussianPdf(
        const double standardDeviationX, const double standardDeviationY, const double correlation,
        const double meanX, const double meanY, const Eigen::VectorXd independentVariables )
{
    using tudat::mathematical_constants::PI;
    using std::pow;
    using std::exp;

    return ( 1.0 / ( 2.0 * PI*standardDeviationX*standardDeviationY*sqrt( 1.0 - pow( correlation, 2.0 ) ) ) ) *
            exp( ( -1.0 / ( 2.0 * ( 1.0 - pow( correlation, 2.0) ) ) ) * (
                     ( pow( independentVariables( 0 ) - meanX, 2.0 ) ) /( pow( standardDeviationX, 2.0 ) )
                     - ( 2.0 * correlation * ( independentVariables( 0 ) - meanX ) *
                         ( independentVariables( 1 )- meanY ) ) / ( standardDeviationX * standardDeviationY )
                     + ( pow( independentVariables( 1 ) - meanY, 2.0 ) ) / ( pow( standardDeviationY, 2.0 ) ) ) );
}


BOOST_AUTO_TEST_SUITE( test_probability_distributions )

//! Test if Multi-dimensional Gaussian distribution class works correctly for 2 dimensions, using analytical solution for
//! bigaussian distribution.
BOOST_AUTO_TEST_CASE( testMultiDimensionalGaussianDistribution )
{
    using namespace tudat::statistics;

    // Defined mean and covariance.
    Eigen::VectorXd mean( 2 );
    mean << 0.0 , 1.0;

    Eigen::MatrixXd covariance( 2, 2 );
    covariance << 3.0, -1.0, -1.0, 3.0;

    // Create distribution
    GaussianDistributionXd distribution( mean, covariance );

    // Defined new mean and covariance.
    Eigen::VectorXd mean2( 2 );
    mean2 << -1.0 , 2.0;

    Eigen::MatrixXd covariance2( 2, 2 );
    covariance2 << 4.0, 1.5, 1.5, 4.0;

    // Create second distribution
    GaussianDistributionXd distribution2( mean2,covariance2 );

    // Define test independent variables.
    Eigen::VectorXd independentVariables( 2 );

    // Test distributions for range of independent variables.
    for( unsigned int i = 0; i < 81; i++ )
    {
        for( unsigned int j = 0; j < 81; j++ )
        {
            independentVariables( 0 ) = -4.0 + 0.1 * static_cast< double >( i );
            independentVariables( 1 ) = -4.0 + 0.1 * static_cast< double >( j );

            {
                // Compute exact solution of pdf
                double standardDeviationX = std::sqrt( covariance( 0, 0 ) );
                double standardDeviationY = std::sqrt( covariance( 1, 1 ) );
                double meanX = mean( 0 );
                double meanY = mean( 1 );
                double correlation = covariance( 0 , 1 ) / ( standardDeviationX * standardDeviationY );

                double exactSolution = computeBiGaussianPdf(
                            standardDeviationX, standardDeviationY, correlation, meanX, meanY, independentVariables );

                // Test pdf value
                BOOST_CHECK_SMALL( std::fabs( distribution.evaluatePdf( independentVariables) - exactSolution ),
                                   std::numeric_limits< double >::epsilon( ) );
            }

            {
                // Compute exact solution of pdf
                double standardDeviationX = std::sqrt( covariance2( 0, 0 ) );
                double standardDeviationY = std::sqrt( covariance2( 1, 1 ) );
                double meanX = mean2( 0 );
                double meanY = mean2( 1 );
                double correlation = covariance2( 0 , 1 ) / ( standardDeviationX * standardDeviationY );

                double exactSolution = computeBiGaussianPdf(
                            standardDeviationX, standardDeviationY, correlation, meanX, meanY, independentVariables );

                // Test pdf value
                BOOST_CHECK_SMALL( std::fabs( distribution2.evaluatePdf( independentVariables ) - exactSolution ) ,
                                   std::numeric_limits<double>::epsilon() ) ;
            }
        }
    }
}

//! Test if Gaussian Copula distribution class works correctly.
BOOST_AUTO_TEST_CASE( testGaussianCopula )
{
    using namespace tudat::statistics;
    using tudat::mathematical_constants::PI;

    // Define properties of distribution.
    int dimension = 2 ;
    Eigen::MatrixXd correlationMatrix( dimension , dimension ) ;
    correlationMatrix << 1.0, 0.3, 0.3, 1.0 ;

    // Create distribution
    GaussianCopulaDistributionXd distribution( correlationMatrix );

    Eigen::VectorXd independentVariables( 2 );
    independentVariables << 0.5 , 0.3 ;

    for( unsigned int i = 1; i < 40; i++ )
    {
        for( unsigned int j = 1; j < 40; j++ )
        {
            independentVariables( 0 ) = 0.025 * static_cast< double >( i );
            independentVariables( 1 ) = 0.025 * static_cast< double >( j );

            // Test Gaussian cupola at current independentVariables/
            // NOTE: This manner of testing the Gaussian cupola (similar computation in source and test file) is not ideal,
            // but no usable test data has been identified.
            BOOST_CHECK_SMALL( std::fabs( gaussianCopulaProbabilityDensity(
                                              independentVariables , dimension , correlationMatrix )
                                          - distribution.evaluatePdf( independentVariables ) ), 1E-15 );

        }
    }

    // Test out of bounds for distribution (pdf equals 0).
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
