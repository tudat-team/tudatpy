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

#include <cmath>
#include <numeric>

#include <Eigen/LU>

#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/lognormal.hpp>

#include "tudat/math/statistics/multiVariateGaussianProbabilityDistributions.h"
namespace tudat
{
namespace statistics
{

GaussianDistributionXd::GaussianDistributionXd( const Eigen::VectorXd& mean, const Eigen::MatrixXd& covarianceMatrix ):
    mean_( mean ), covarianceMatrix_( covarianceMatrix )
{
    if( covarianceMatrix.rows( ) != covarianceMatrix.cols( ) )
    {
        throw std::runtime_error( "Error, covarianceMatrix input to GaussianDistributionXd is not square" );
    }

    dimension_ = static_cast< double >( mean_.rows( ) );
    determinant_ = covarianceMatrix_.determinant( );
    inverseCovarianceMatrix_ = covarianceMatrix_.inverse( );
}

double GaussianDistributionXd::evaluatePdf( const Eigen::VectorXd& independentVariables )
{
    const Eigen::VectorXd distanceFromMean = ( independentVariables - mean_ );
    const Eigen::MatrixXd location = -0.5 * ( distanceFromMean.transpose( ) * inverseCovarianceMatrix_ * distanceFromMean );

    return std::exp( location( 0, 0 ) ) / ( std::pow( 2.0 * mathematical_constants::PI, dimension_ / 2.0 ) * std::sqrt( determinant_ ) );
}

GaussianCopulaDistributionXd::GaussianCopulaDistributionXd( const Eigen::MatrixXd& correlationMatrix ):
    correlationMatrix_( correlationMatrix )
{
    if( correlationMatrix.rows( ) != correlationMatrix.cols( ) )
    {
        throw std::runtime_error( "Error, correlationMatrix input to GaussianCopulaDistributionXd is not square" );
    }

    dimension_ = correlationMatrix.rows( );
    inverseCorrelationMatrix_ = correlationMatrix_.inverse( );
    determinant_ = correlationMatrix_.determinant( );
}

//! Function to evaluate pdf of Gaussian cupola distribution
double GaussianCopulaDistributionXd::evaluatePdf( const Eigen::VectorXd& independentVariables )
{
    double probabilityDensity = 0.0;

    // Check if vector independentVariables is inside [0,1]
    int inBound = 0;
    for( int i = 0; i < dimension_; i++ )
    {
        if( independentVariables( i ) > 0.0 && independentVariables( i ) < 1.0 )
        {
            inBound++;
        }
    }

    // If data is in bounds
    if( inBound == dimension_ )
    {
        // Convert U[0,1] to N[0,1] using inverse CDF of standard normal distribution
        Eigen::VectorXd gaussianQuantiles( dimension_ );
        boost::math::normal distribution( 0.0, 1.0 );

        for( int i = 0; i < dimension_; i++ )
        {
            gaussianQuantiles( i ) = boost::math::quantile( distribution, independentVariables( i ) );  // Inverse cdf
        }

        // Calculate probability density
        Eigen::MatrixXd location = -0.5 *
                ( gaussianQuantiles.transpose( ) * ( inverseCorrelationMatrix_ - Eigen::MatrixXd::Identity( dimension_, dimension_ ) ) *
                  gaussianQuantiles );

        probabilityDensity = ( ( 1.0 / ( std::sqrt( determinant_ ) ) ) * std::exp( location( 0, 0 ) ) );
    }

    return probabilityDensity;
}

}  // namespace statistics
}  // namespace tudat
