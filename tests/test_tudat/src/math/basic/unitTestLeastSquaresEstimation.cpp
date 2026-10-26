/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <cmath>
#include <random>

#include <boost/test/included/unit_test.hpp>

#include <Eigen/Sparse>

#include "tudat/math/basic/leastSquaresEstimation.h"

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_least_squares_estimation )

//! Test that the consider parameter covariance contribution for a sparse design matrix equals the dense version
BOOST_AUTO_TEST_CASE( testSparseConsiderParametersCovarianceContribution )
{
    const int numberOfObservations = 400;
    const int numberOfParameters = 15;
    const int numberOfConsiderParameters = 3;

    std::mt19937 generator( 7 );
    std::normal_distribution< double > normalDistribution( 0.0, 1.0 );
    std::uniform_real_distribution< double > uniformDistribution( 0.0, 1.0 );

    // Design matrix with about 10% non-zero entries
    Eigen::MatrixXd designMatrix = Eigen::MatrixXd::Zero( numberOfObservations, numberOfParameters );
    for( int i = 0; i < numberOfObservations; i++ )
    {
        for( int j = 0; j < numberOfParameters; j++ )
        {
            if( uniformDistribution( generator ) < 0.1 )
            {
                designMatrix( i, j ) = normalDistribution( generator );
            }
        }
    }
    Eigen::SparseMatrix< double > sparseDesignMatrix = designMatrix.sparseView( );

    // Dense consider partials, positive weights and symmetric positive (semi-)definite covariances
    Eigen::MatrixXd considerDesignMatrix = Eigen::MatrixXd::Random( numberOfObservations, numberOfConsiderParameters );
    Eigen::VectorXd weightsDiagonal = ( Eigen::VectorXd::Random( numberOfObservations ).array( ) + 2.0 ).matrix( );
    Eigen::MatrixXd covarianceFactor = Eigen::MatrixXd::Random( numberOfParameters, numberOfParameters );
    Eigen::MatrixXd normalisedCovariance =
            covarianceFactor * covarianceFactor.transpose( ) + Eigen::MatrixXd::Identity( numberOfParameters, numberOfParameters );
    Eigen::MatrixXd considerCovarianceFactor = Eigen::MatrixXd::Random( numberOfConsiderParameters, numberOfConsiderParameters );
    Eigen::MatrixXd considerCovariance = considerCovarianceFactor * considerCovarianceFactor.transpose( );

    Eigen::MatrixXd denseContribution = linear_algebra::calculateConsiderParametersCovarianceContribution< Eigen::MatrixXd >(
            normalisedCovariance, designMatrix, weightsDiagonal, considerDesignMatrix, considerCovariance );
    Eigen::MatrixXd sparseContribution = linear_algebra::calculateConsiderParametersCovarianceContribution(
            normalisedCovariance, sparseDesignMatrix, weightsDiagonal, considerDesignMatrix, considerCovariance );

    BOOST_CHECK_EQUAL( sparseContribution.rows( ), numberOfParameters );
    BOOST_CHECK_EQUAL( sparseContribution.cols( ), numberOfParameters );

    // The two versions differ only in the order of the products
    const double maximumEntry = denseContribution.cwiseAbs( ).maxCoeff( );
    for( int i = 0; i < numberOfParameters; i++ )
    {
        for( int j = 0; j < numberOfParameters; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( sparseContribution( i, j ) - denseContribution( i, j ) ), 1.0E-12 * maximumEntry );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
