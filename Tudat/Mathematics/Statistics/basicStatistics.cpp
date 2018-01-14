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

#include "Tudat/Mathematics/Statistics/basicStatistics.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

namespace tudat
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
    double averageOfComponents = computeAverageOfVectorComponents( vectorOfData );

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

//! Compute sample mean for a sample of VectorXd
Eigen::VectorXd computeSampleMean( const std::vector< Eigen::VectorXd >& sampleData )
{
    // Return sample mean.
    Eigen::VectorXd initialValue = Eigen::VectorXd::Zero( sampleData.at( 0 ).rows( ) );
    return std::accumulate( sampleData.begin( ), sampleData.end( ), initialValue )
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

//! Compute Sample median
double computeSampleMedian( std::vector< double > sampleData )
{
    // Sort data
    std::sort( sampleData.begin() , sampleData.end());

    // Check if odd number of samples or even
    double numberOfSamples = static_cast< double >( sampleData.size() );
    int odd = tudat::basic_mathematics::computeModulo( numberOfSamples , 2.0 ) ;

    // Calculate sample median
    double sampleMedian;
    if( odd == 0 ) // even
    {
        // 0 .. 99 (100) ->
        int index = static_cast< int >( numberOfSamples / 2.0 - 0.5 ) ;
        sampleMedian = (sampleData[index] + sampleData[index+1])/2.0 ;
    }
    else // odd
    {
        // 0 .. 100 (101) -> 50
        int index = static_cast< int >( numberOfSamples / 2.0 - 0.5 ) ;
        sampleMedian = sampleData[index] ;
    }

    return sampleMedian;
}

//! Compute sample variance for a sample of VectorXd.
Eigen::VectorXd computeSampleVariance( const std::vector< Eigen::VectorXd >& sampleData )
{
    // Declare local variables.
    // Declare and compute sample mean.
    Eigen::VectorXd sampleMean_ = computeSampleMean( sampleData );

    // Declare and initialize sum of residuals squared.
    Eigen::VectorXd sumOfResidualsSquared_ = Eigen::VectorXd::Zero( sampleData.at( 0 ).rows( ) );

    // Compute sum of residuals of sample data squared.
    for ( unsigned int i = 0; i < sampleData.size( ); i++ )
    {
        sumOfResidualsSquared_ += ( sampleData.at( i ) - sampleMean_ ).cwiseProduct( sampleData.at( i ) - sampleMean_ );
    }

    // Return sample variance.
    return 1.0 / ( static_cast< double >( sampleData.size( ) ) - 1.0 ) * sumOfResidualsSquared_;
}

} // namespace statistics
} // namespace tudat
