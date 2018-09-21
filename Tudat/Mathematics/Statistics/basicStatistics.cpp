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
#include "Tudat/Basics/utilities.h"

namespace tudat
{

namespace statistics
{

//! Compute average of the components of a vector.
double computeAverageOfVectorComponents( const Eigen::VectorXd& vectorOfData )
{
    return vectorOfData.sum( ) / vectorOfData.rows( );
}

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
    double sampleMean = computeSampleMean( sampleData );

    // Declare and initialize sum of residuals squared.
    double sumOfResidualsSquared = 0.0;

    // Compute sum of residuals of sample data squared.
    for ( unsigned int i = 0; i < sampleData.size( ); i++ )
    {
        sumOfResidualsSquared += std::pow( sampleData.at( i ) - sampleMean, 2.0 );
    }

    // Return sample variance.
    return ( 1.0 / ( static_cast< double >( sampleData.size( ) ) - 1.0 ) * sumOfResidualsSquared );
}

//! Compute Sample median
double computeSampleMedian( std::vector< double > sampleData )
{
    // Sort data
    std::sort( sampleData.begin() , sampleData.end());

    // Check if odd number of samples or even
    double numberOfSamples = static_cast< double >( sampleData.size( ) );
    int odd = static_cast< int >( basic_mathematics::computeModulo( numberOfSamples, 2.0 ) );

    // Calculate sample median
    double sampleMedian;
    if( odd == 0 ) // even
    {
        // 0 .. 99 (100) ->
        unsigned int index = static_cast< unsigned int >( numberOfSamples / 2.0 - 0.5 );
        sampleMedian = ( sampleData[ index ] + sampleData[ index + 1 ] ) / 2.0;
    }
    else // odd
    {
        // 0 .. 100 (101) -> 50
        unsigned int index = static_cast< unsigned int >( numberOfSamples / 2.0 - 0.5 );
        sampleMedian = sampleData[ index ];
    }

    return sampleMedian;
}

//! Compute sample variance for a sample of VectorXd.
Eigen::VectorXd computeSampleVariance( const std::vector< Eigen::VectorXd >& sampleData )
{
    // Declare local variables.
    // Declare and compute sample mean.
    Eigen::VectorXd sampleMean = computeSampleMean( sampleData );

    // Declare and initialize sum of residuals squared.
    Eigen::VectorXd sumOfResidualsSquared = Eigen::VectorXd::Zero( sampleData.at( 0 ).rows( ) );

    // Compute sum of residuals of sample data squared.
    for ( unsigned int i = 0; i < sampleData.size( ); i++ )
    {
        sumOfResidualsSquared += ( sampleData.at( i ) - sampleMean ).cwiseProduct( sampleData.at( i ) - sampleMean );
    }

    // Return sample variance.
    return 1.0 / ( static_cast< double >( sampleData.size( ) ) - 1.0 ) * sumOfResidualsSquared;
}

//! Compute moving average of an Eigen vector.
Eigen::VectorXd computeMovingAverage( const Eigen::VectorXd& sampleData, const unsigned int numberOfAveragingPoints )
{
    // Declare moving average vector
    int numberOfSamplePoints = sampleData.rows( );
    Eigen::VectorXd movingAverage;
    movingAverage.resize( numberOfSamplePoints );

    // Check that number of samples to be used for averaging is odd
    if ( ( numberOfAveragingPoints % 2 ) == 0 )
    {
        throw std::runtime_error( "Error while computing moving average. The number of points used to average has to be odd." );
    }
    int numberOfPointsOnEachSide = ( numberOfAveragingPoints - 1 ) / 2;

    // Loop over each element and compute moving average
    unsigned int lowerIndex, upperIndex;
    for ( int i = 0; i < numberOfSamplePoints; i++ )
    {
        lowerIndex = ( ( i - numberOfPointsOnEachSide ) < 0 ) ? 0 : ( i - numberOfPointsOnEachSide );
        upperIndex = ( ( i + numberOfPointsOnEachSide ) > ( numberOfSamplePoints - 1 ) ) ?
                  ( numberOfSamplePoints - 1 ) : ( i + numberOfPointsOnEachSide );
        movingAverage[ i ] = computeAverageOfVectorComponents( sampleData.segment( lowerIndex, ( upperIndex - lowerIndex + 1 ) ) );
    }
    return movingAverage;
}

//! Compute moving average of a set of Eigen vectors in a STL vector.
std::vector< Eigen::Vector3d > computeMovingAverage(
        const std::vector< Eigen::Vector3d >& sampleData, const unsigned int numberOfAveragingPoints )
{
    // Convert map to Eigen matrix
    Eigen::MatrixXd matrixOfSampleData = utilities::convertStlVectorToEigenMatrix< double, 3 >( sampleData );

    // Declare moving average vector
    Eigen::MatrixXd movingAverage;
    movingAverage.resizeLike( matrixOfSampleData );

    // Loop over rows and compute sample mean
    for ( int i = 0; i < movingAverage.rows( ); i++ )
    {
        movingAverage.row( i ) = computeMovingAverage( matrixOfSampleData.row( i ).transpose( ), numberOfAveragingPoints ).transpose( );
    }

    // Output data as vector
    std::vector< Eigen::Vector3d > outputData;
    for ( unsigned int i = 0; i < movingAverage.cols( ); i++ )
    {
        outputData.push_back( movingAverage.col( i ) );
    }
    return outputData;
}

//! Compute moving average of a set of Eigen vectors in a map.
std::map< double, Eigen::VectorXd > computeMovingAverage(
        const std::map< double, Eigen::VectorXd >& sampleData, const unsigned int numberOfAveragingPoints )
{
    // Convert map to Eigen matrix
    Eigen::MatrixXd matrixOfSampleData =
            utilities::extractKeyAndValuesFromMap< double, double, Eigen::Dynamic >( sampleData ).second;

    // Declare moving average vector
    Eigen::MatrixXd movingAverage;
    movingAverage.resizeLike( matrixOfSampleData );

    // Loop over rows and compute sample mean
    for ( int i = 0; i < movingAverage.rows( ); i++ )
    {
        movingAverage.row( i ) = computeMovingAverage( matrixOfSampleData.row( i ).transpose( ), numberOfAveragingPoints ).transpose( );
    }

    // Output data as map
    unsigned int i = 0;
    std::map< double, Eigen::VectorXd > outputData;
    for ( typename std::map< double, Eigen::VectorXd >::const_iterator
          mapIterator = sampleData.begin( ); mapIterator != sampleData.end( ); mapIterator++, i++ )
    {
        outputData[ mapIterator->first ] = movingAverage.col( i );
    }
    return outputData;
}

} // namespace statistics

} // namespace tudat
