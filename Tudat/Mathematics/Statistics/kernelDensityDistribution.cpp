/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/Mathematics/Statistics/kernelDensityDistribution.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <Tudat/Mathematics/Statistics/basicStatistics.h>

#include <Tudat/InputOutput/matrixTextFileReader.h>

namespace tudat
{

namespace statistics
{

using tudat::mathematical_constants::PI;

//! Get probability density
double EpanechnikovKernelDistribution::evaluatePdf( const double& independentVariable )
{
    if( ( independentVariable - mean_ ) >= ( -bandWidth_ ) && ( independentVariable - mean_ ) <= ( bandWidth_ ) )
    {
        return ( 3.0 / ( 4.0 * bandWidth_ ) ) * ( 1.0 - std::pow( ( independentVariable - mean_ ) / bandWidth_, 2.0 ) );
    }
    else
    {
        return 0.0;
    }
}

//! Get probability mass
double EpanechnikovKernelDistribution::evaluateCdf( const double& independentVariable )
{
    if( ( independentVariable - mean_ ) >= ( -bandWidth_ ) && ( independentVariable - mean_ ) <= ( bandWidth_ ) )
    {
        return ( 3.0 / ( 4.0 * bandWidth_ ) ) *
                ( ( independentVariable - mean_ ) - ( std::pow( ( independentVariable - mean_ ), 3.0 )
                                                      / ( 3.0 * std::pow( bandWidth_, 2.0 ) ) ) ) + 0.5;
    }
    else if( ( independentVariable - mean_ ) < ( -bandWidth_ ) )
    {
        return 0.0;
    }
    else
    {
        return 1.0;
    }

}

//! Constructor
KernelDensityDistribution::KernelDensityDistribution(
        const std::vector< Eigen::VectorXd >& samples,
        const double bandWidthFactor,
        const KernelType kernelType,
        const Eigen::VectorXd& standardDeviation,
        const Eigen::VectorXd& manualBandwidth )
{
    // Load data
    dataSamples_ = samples;
    dimensions_ = samples[ 0 ].rows( );
    numberOfSamples_ = static_cast< int >( samples.size( ) );

    // Compute datasample properties
    computeSampleMean( );
    computeSampleVariance( );

    // Adjust samples with standard deviation and recompute datasample properties
    if( standardDeviation.rows( ) > 0 )
    {
        scaleSamplesWithVariance( standardDeviation );
        computeSampleMean( );
        computeSampleVariance( );
    }


    // Compute bandwidth
    if( manualBandwidth.rows( ) > 0 )
    {
        bandWidth_ = manualBandwidth * bandWidthFactor;
    }
    else
    {
        computeOptimalBandWidth( );
        bandWidth_ = optimalBandwidth_ * bandWidthFactor;
    }

    // Construct kernel matrix rows: samples, cols: dimensions_
    kernelType_ = kernelType;



    generateKernelPointerMatrix( );
}

//! Set kernel pointer matrix
void KernelDensityDistribution::generateKernelPointerMatrix( )
{
    std::vector< std::vector< boost::shared_ptr< ContinuousProbabilityDistribution< double > > > > matrix( 0 );
    kernelPointersMatrix_ = matrix;

    for( unsigned int i = 0; i < bandWidth_.rows( ); i++ )
    {
        if( bandWidth_( i ) < 10.0 * std::numeric_limits< double >::epsilon( ) )
        {
            std::string errorMessage = "Error in kernel density distribution, index " +
                    boost::lexical_cast< std::string >( i ) + "has bandwidth " +
                    boost::lexical_cast< std::string >( bandWidth_( i ) );
            throw std::runtime_error( errorMessage );
        }
    }
    // Fill kernel pointer matrix with distribution pointer objects
    std::vector< boost::shared_ptr< ContinuousProbabilityDistribution< double > > > vector( dataSamples_[ 0 ].rows( ) );
    for( unsigned int i = 0; i < dataSamples_.size( ); i++ )
    {
        for( unsigned int j = 0; j < vector.size( ); j++ )
        {
            vector[ j ] = constructKernel( dataSamples_[ i ]( j ), bandWidth_( j ) );
        }
        kernelPointersMatrix_.push_back( vector );
    }
}

void KernelDensityDistribution::computeSampleMean( )
{
    sampleMean_ = Eigen::VectorXd::Zero( dimensions_ );
    for( unsigned int i = 0; i < dataSamples_.size( ); i++ )
    {
        sampleMean_ +=  dataSamples_[ i ];
    }
    sampleMean_ = sampleMean_ / static_cast< double >( numberOfSamples_ );
}

void KernelDensityDistribution::computeSampleVariance( )
{
    sampleVariance_ = Eigen::VectorXd::Zero( dimensions_ );
    for( unsigned int i = 0; i < dataSamples_.size( ); i++ )
    {
        sampleVariance_ += ( dataSamples_[ i ] - sampleMean_ ).cwiseAbs2( );
    }
    sampleVariance_ = sampleVariance_ / ( static_cast< double >( numberOfSamples_ ) - 1.0 );
    sampleStandardDeviation_ = sampleVariance_.cwiseSqrt( );
}

//! Scale the samples such to the required variance
void KernelDensityDistribution::scaleSamplesWithVariance(
        const Eigen::VectorXd& standardDeviation )
{
    for( unsigned int i = 0; i < dataSamples_.size( ); i++ )
    {
        for( int j = 0; j < dataSamples_[ i ].rows( ); j++ )
        {
            dataSamples_[ i ]( j ) = dataSamples_[ i ]( j ) / sampleStandardDeviation_( j ) * standardDeviation( j );
        }
    }
}

//! Compute the optimal bandwidth
void KernelDensityDistribution::computeOptimalBandWidth( )
{
    int dims = dimensions_;

    optimalBandwidth_ = Eigen::VectorXd::Zero( dims, 1);

    // Calculate sigma (median absolute deviation estimator)
    Eigen::VectorXd medianOfSamples = getMedian( dataSamples_ );

    std::vector< Eigen::VectorXd > dataSamples2( dataSamples_.size( ) );
    for( unsigned int i = 0; i < dataSamples_.size( ); i++ )
    {
        dataSamples2[ i ] = dataSamples_[ i ] - medianOfSamples;
        dataSamples2[ i ] = dataSamples2[ i ].cwiseAbs( );
    }
    Eigen::VectorXd medianOfSamples2 = getMedian( dataSamples2 );
    Eigen::VectorXd sigma = medianOfSamples2 / 0.6745;

    // Calculate optimal Bandwidth
    optimalBandwidth_ = std::pow( 4.0 / ( ( static_cast< double >( dimensions_ ) + 2.0 ) *
                                         static_cast< double >( numberOfSamples_ ) ),
                                 1.0 / ( static_cast< double >( dimensions_ )  + 4.0 ) ) * sigma;
}

//! Get median from datasamples
Eigen::VectorXd KernelDensityDistribution::getMedian( std::vector< Eigen::VectorXd > samples )
{
    int dims = samples[ 0 ].rows( );
    std::vector< std::vector< double > > dataSamplesVectors( dims );

    std::vector< double > data( 0 );
    for( int i = 0; i < dims; i++ )
    {
        for( unsigned int j = 0; j < samples.size( ); j++ )
        {
            data.push_back( samples[ j ](i) );
        }
        dataSamplesVectors[ i ] = data;
        data.clear( );
    }

    Eigen::VectorXd medianOfSamples( dims );
    for( int i = 0; i < dims; i++ )
    {
        medianOfSamples[ i ] = tudat::statistics::computeSampleMedian( dataSamplesVectors[ i ] );
    }
    return medianOfSamples;
}

//! Get probability density of the kernel density distribution
double KernelDensityDistribution::evaluatePdf( const Eigen::VectorXd& independentVariables )
{
    double propbabilityDensity = 0.0;
    double value = 1.0;
    for( int i = 0; i < numberOfSamples_; i++ )
    {
        for( int dim = 0; dim < dimensions_; dim++ )
        {
            value = value * kernelPointersMatrix_[ i ][ dim ]->evaluatePdf( independentVariables( dim ) );
        }
        propbabilityDensity += value;
        value = 1.0;
    }
    return propbabilityDensity / static_cast< double >( numberOfSamples_ );
}

//! Get cumulative probability of the kernel density distribution
double KernelDensityDistribution::evaluateCdf( const Eigen::VectorXd& independentVariables )
{
    double cumulativeProbability = 0.0;
    double value = 1.0;
    for( int i = 0; i < numberOfSamples_; i++ )
    {
        for( int dim = 0; dim < dimensions_; dim++ )
        {
            value = value * kernelPointersMatrix_[ i ][ dim ]->evaluateCdf( independentVariables( dim ) );
        }
        cumulativeProbability += value;
        value = 1.0;
    }
    return cumulativeProbability / static_cast< double >( numberOfSamples_ );
}

//! Get cumulative probability of marginal distribution
double KernelDensityDistribution::getCumulativeMarginalProbability(int marginalDimension, double x )
{
    double cumulativeProbability = 0.0;
    for( int i = 0; i < numberOfSamples_; i++ )
    {
        cumulativeProbability += kernelPointersMatrix_[ i ][ marginalDimension ]->evaluateCdf( x );
    }
    return cumulativeProbability / static_cast< double >( numberOfSamples_ );
}

//! Get probability density of marginal distribution
double KernelDensityDistribution::getMarginalProbabilityDensity(
        const std::vector< int >& marginalDimensions, const Eigen::VectorXd& independentVariables )
{
    double probabilityDensity = 0.0;
    double value = 1.0;
    for( int i = 0; i < numberOfSamples_; i++ )
    {
        for( unsigned int j = 0; j < marginalDimensions.size( ); j++ )
        {
            value = value * kernelPointersMatrix_[ i ][ marginalDimensions[ j ] ]->evaluatePdf( independentVariables( j ) );
        }

        probabilityDensity += value;
        value = 1.0;
    }
    return probabilityDensity / static_cast< double >( numberOfSamples_ );
}

//! Get probability density of marginal distribution
double KernelDensityDistribution::getMarginalProbabilityDensityd(
        const int marginalDimension, const double independentVariable )
{
    double probabilityDensity = 0.0;
    for( int i = 0; i < numberOfSamples_; i++ )
    {
        probabilityDensity += kernelPointersMatrix_[ i ][ marginalDimension ]->evaluatePdf( independentVariable );
    }
    return probabilityDensity / static_cast< double >( numberOfSamples_ );
}

//! get cumulative probability of marginal distribution of conditional distribution
double KernelDensityDistribution::getCumulativeConditionalMarginalProbability(
        const std::vector< int >& conditionDimensions,
        const std::vector< double >& conditions,
        const int marginalDimension, const double independentVariable )
{
    double value = 1.0;
    double normalizationFactor = 0.0;
    double marginalValue = 0.0;

    // multiply all conditions * marginal, divide marginals of conditions
    for( int i = 0; i < numberOfSamples_; i++ )
    {
        value = 1.0;

        for( unsigned int j = 0; j < conditionDimensions.size( ); j++ )
        {
            value = value * kernelPointersMatrix_[ i ][ conditionDimensions[ j ] ]->evaluateCdf( conditions[ j ] );
        }
        normalizationFactor += value;

        value = value * kernelPointersMatrix_[ i ][ marginalDimension ]->evaluateCdf( independentVariable );
        marginalValue += value;
    }

    return marginalValue / normalizationFactor;
}

//! get probability density of marginal distribution of conditional distribution
double KernelDensityDistribution::getConditionalMarginalProbabilityDensity(
        const std::vector< int >& conditionDimensions,
        const std::vector< double >& conditions,
        const int marginalDimension, const double independentVariable )
{
    double value = 1.0;
    double normalizationFactor = 0.0;
    double marginalValue = 0.0;

    // multiply all conditions * marginal, divide marginals of conditions
    for( int i = 0; i < numberOfSamples_; i++ )
    {
        value = 1.0;

        for( unsigned int j = 0; j < conditionDimensions.size( ); j++ )
        {
            value = value * kernelPointersMatrix_[ i ][ conditionDimensions[ j ] ]->evaluatePdf( conditions[ j ] );
        }

        normalizationFactor += value;

        value = value * kernelPointersMatrix_[ i ][ marginalDimension ]->evaluatePdf( independentVariable );
        marginalValue += value;
    }

    return marginalValue / ( normalizationFactor );
}

} // Close Namespace statistics

} // Close Namespace tudat

