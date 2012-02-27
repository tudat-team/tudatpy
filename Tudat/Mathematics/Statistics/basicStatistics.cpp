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
 *      090807    J. Melman         First creation of code.
 *      100930    D. Dirkx          Modified to comply with Tudat standards
 *      100930    J. Melman         Implemented namespace, minor comment
 *                                  changes.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120202    K. Kumar          Moved functions from linearAlgebra.cpp.
 *
 *    References
 *
 */

#include <cmath>
#include <numeric>
#include "Tudat/Mathematics/Statistics/basicStatistics.h"

namespace tudat
{
namespace mathematics
{
namespace statistics
{

//! Compute average of the components of a vector.
double computeAverageOfVectorComponents( const Eigen::VectorXd& vectorOfData )
{ return vectorOfData.sum( ) / vectorOfData.rows( ); }

//! Compute standard deviation of the components of a vector.
double computeStandardDeviationOfVectorComponents( const Eigen::VectorXd& vectorOfData )
{
    // Compute average of components.
    double averageOfComponents = tudat::mathematics::statistics::
            computeAverageOfVectorComponents( vectorOfData );

    // Declare variance of components.
    double varianceOfComponents = 0.0;

    // Compute variance of components.
    for ( int i = 0 ; i < vectorOfData.rows( ) ; i++ )
    {
        varianceOfComponents += std::pow( ( vectorOfData( i ) - averageOfComponents ), 2.0 );
    }

    varianceOfComponents /= static_cast< double >( vectorOfData.rows( ) - 1 );

    // Return square root of variance ( = standard deviation ).
    return std::sqrt( varianceOfComponents );
}

//! Compute sample mean.
double computeSampleMean( const std::vector< double >& sampleData )
{
    // Return sample mean.
    return std::accumulate( sampleData.begin( ), sampleData.end( ), 0.0 )
            / static_cast< double >( sampleData.size( ) );
}

//! Compute sample variance.
double computeSampleVariance( const std::vector< double >& sampleData )
{
    // Declare local variables.
    // Declare and compute sample mean.
    double sampleMean_ = computeSampleMean( sampleData );

    // Declare and initialize sum of residuals squared.
    double sumOfResidualsSquared_ = 0.0;

    // Compute sum of residuals of sample data squared.
    for ( unsigned int i = 0; i < sampleData.size( ); i++ )
    {
        sumOfResidualsSquared_ += std::pow( sampleData.at( i ) - sampleMean_, 2.0 );
    }

    // Return sample variance.
    return 1.0 / ( static_cast< double >( sampleData.size( ) ) - 1.0 ) * sumOfResidualsSquared_;
}

} // namespace statistics
} // namespace mathematics
} // namespace tudat
