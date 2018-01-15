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

#include <cmath>
#include <numeric>

#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/lognormal.hpp>

#include "Tudat/Mathematics/Statistics/multiVariateGaussianProbabilityDistributions.h"
namespace tudat
{
namespace statistics
{


//! Function to evaluate pdf of Gaussian cupola distribution
double GaussianCopulaDistributionXd::evaluatePdf(
        const Eigen::VectorXd& independentVariables )
{
    double probabilityDensity = 0.0 ;

    // Check if vector independentVariables is inside [0,1]
    int inBound = 0 ;
    for( int i = 0 ; i < dimension_ ; i++ )
    {
        if( independentVariables(i) > 0.0 && independentVariables(i) < 1.0 )
        {
            inBound++;
        }
    }

    // If data is in bounds
    if( inBound == dimension_ )
    {
        // Convert U[0,1] to N[0,1] using inverse CDF of standard normal distribution
        Eigen::VectorXd gaussianQuantiles( dimension_ ) ;
        boost::math::normal distribution( 0.0 , 1.0 );

        for( int i = 0 ; i < dimension_ ; i++ )
        {
            gaussianQuantiles( i ) = boost::math::quantile( distribution , independentVariables( i ) ); // Inverse cdf
        }

        // Calculate probability density
        Eigen::MatrixXd location = - 0.5 * ( gaussianQuantiles.transpose() *
                       ( inverseCorrelationMatrix_ - Eigen::MatrixXd::Identity( dimension_ , dimension_ ) ) *
                                             gaussianQuantiles );

        probabilityDensity = ( ( 1.0 / ( std::sqrt( determinant_ ) ) ) * std::exp( location( 0, 0 ) ) ) ;
    }

    return probabilityDensity ;
}

} // namespace statistics
} // namespace tudat
