/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Mathematics/Statistics/kernelDensityDistribution.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/Statistics/basicStatistics.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"
namespace tudat
{

namespace statistics
{

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
    dimensions_ = dataSamples_.at( 0 ).rows( );
    numberOfSamples_ = static_cast< int >( samples.size( ) );

    // Check input consistency
    if( ( standardDeviation.rows( ) > 0 ) && ( standardDeviation.rows( ) != dimensions_  ) )
    {
        std::string errorMessage = "Error when creating KernelDensityDistribution, manual standard deviation size is inconsistent, should have size : " +
                std::to_string( dimensions_ ) + "but has size " +
                std::to_string( standardDeviation.rows( ) );
        throw std::runtime_error( errorMessage );

    }

    if( ( manualBandwidth.rows( ) > 0 ) && ( manualBandwidth.rows( ) != dimensions_  ) )
    {
        std::string errorMessage = "Error when creating KernelDensityDistribution, manual bandwidth size is inconsistent, should have size : " +
                std::to_string( dimensions_ ) + "but has size " +
                std::to_string( manualBandwidth.rows( ) );
        throw std::runtime_error( errorMessage );
    }

    for( unsigned int i = 0; i < samples.size( ); i++ )
    {
        if( dataSamples_.at( i ).rows( ) != dimensions_ )
        {
            std::string errorMessage = "Error when creating KernelDensityDistribution, samples size is inconsistent, should have size : " +
                    std::to_string( dimensions_ ) + "but entry " +
                    std::to_string( i ) + " has size " +
                    std::to_string( dataSamples_.at( i ).rows( ) );
            throw std::runtime_error( errorMessage );
        }
    }

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

//! Function that generates the kernel density distribution based on the samples and kernel type that is provided
void KernelDensityDistribution::generateKernelPointerMatrix( )
{
    // Clear existing kernels.
    kernelPointersMatrix_.clear( );

    // Check for numerical problems with bandwidths.
    for( int i = 0; i < bandWidth_.rows( ); i++ )
    {
        if( bandWidth_( i ) < 10.0 * std::numeric_limits< double >::epsilon( ) )
        {
            std::string errorMessage = "Error in kernel density distribution, index " +
                    std::to_string( i ) + "has bandwidth " +
                    std::to_string( bandWidth_( i ) );
            throw std::runtime_error( errorMessage );
        }
    }

    // Fill kernel pointer matrix with distribution pointer objects
    std::vector< boost::shared_ptr< ContinuousProbabilityDistribution< double > > > vector( dataSamples_[ 0 ].rows( ) );

    // Iterate over all samples; create kernel for each entry in each sample
    for( unsigned int i = 0; i < dataSamples_.size( ); i++ )
    {
        for( unsigned int j = 0; j < vector.size( ); j++ )
        {
            vector[ j ] = constructKernel( dataSamples_[ i ]( j ), bandWidth_( j ) );
        }
        kernelPointersMatrix_.push_back( vector );
    }
}

//! Function that computes and sets the sample mean.
void KernelDensityDistribution::computeSampleMean( )
{
    sampleMean_ = Eigen::VectorXd::Zero( dimensions_ );
    for( unsigned int i = 0; i < dataSamples_.size( ); i++ )
    {
        sampleMean_ +=  dataSamples_[ i ];
    }
    sampleMean_ = sampleMean_ / static_cast< double >( numberOfSamples_ );
}

//! Function that computes and sets the sample variance and standard deviation.
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

//! Function to scale the samples to the required standard deviation
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
    optimalBandwidth_ = Eigen::VectorXd::Zero( dimensions_, 1);

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

//! Function to retrieve the sample set median
Eigen::VectorXd KernelDensityDistribution::getMedian(
        const std::vector< Eigen::VectorXd >& samples )
{
    int sampleDimensions = samples[ 0 ].rows( );
    std::vector< std::vector< double > > dataSamplesVectors( sampleDimensions );

    // Create vector of j^th entry from each sample
    std::vector< double > data( 0 );
    for( int i = 0; i < sampleDimensions; i++ )
    {
        for( unsigned int j = 0; j < samples.size( ); j++ )
        {
            data.push_back( samples[ j ]( i ) );
        }
        dataSamplesVectors[ i ] = data;
        data.clear( );
    }

    // Compute sample median for each entry of samples.
    Eigen::VectorXd medianOfSamples( sampleDimensions );
    for( int i = 0; i < sampleDimensions; i++ )
    {
        medianOfSamples[ i ] = tudat::statistics::computeSampleMedian( dataSamplesVectors[ i ] );
    }
    return medianOfSamples;
}

//! Get probability density of the kernel density distribution
double KernelDensityDistribution::evaluatePdf( const Eigen::VectorXd& independentVariables )
{
    double propbabilityDensity = 0.0;
    double currentKernelPdf = 1.0;

    // Iterative over all kernels
    for( int i = 0; i < numberOfSamples_; i++ )
    {     
        currentKernelPdf = 1.0;

        // Compute pdf of current kernel
        for( int currentDimension = 0; currentDimension < dimensions_; currentDimension++ )
        {
            currentKernelPdf *= kernelPointersMatrix_[ i ][ currentDimension ]->evaluatePdf(
                        independentVariables( currentDimension ) );
        }
        propbabilityDensity += currentKernelPdf;
    }

    // Average over all kernels
    return propbabilityDensity / static_cast< double >( numberOfSamples_ );
}

//! Get cumulative probability of the kernel density distribution
double KernelDensityDistribution::evaluateCdf( const Eigen::VectorXd& independentVariables )
{
    double cumulativeProbability = 0.0;
    double currentKernelCdf = 1.0;

    // Iterative over all kernels
    for( int i = 0; i < numberOfSamples_; i++ )
    {
        currentKernelCdf = 1.0;

        // Compute cdf of current kernel
        for( int currentDimension = 0; currentDimension < dimensions_; currentDimension++ )
        {
            currentKernelCdf *= kernelPointersMatrix_[ i ][ currentDimension ]->evaluateCdf(
                        independentVariables( currentDimension ) );
        }
        cumulativeProbability += currentKernelCdf;
    }

    // Average over all kernels
    return cumulativeProbability / static_cast< double >( numberOfSamples_ );
}

//! Get cumulative probability of marginal distribution
double KernelDensityDistribution::evaluateCumulativeMarginalProbability(
        const int marginalDimension, const double independentVariable )
{
    // Compute cdf at independentVariable in marginalDimension, averaged over all samples
    double cumulativeProbability = 0.0;
    for( int i = 0; i < numberOfSamples_; i++ )
    {
        cumulativeProbability += kernelPointersMatrix_[ i ][ marginalDimension ]->evaluateCdf( independentVariable );
    }
    return cumulativeProbability / static_cast< double >( numberOfSamples_ );
}

//! Function to evaluate probability density of joint marginal distribution.
double KernelDensityDistribution::evaluateMarginalProbabilityDensity(
        const std::vector< int >& marginalDimensions, const Eigen::VectorXd& independentVariables )
{
    double probabilityDensity = 0.0;
    double marginalPdfOfCurrentKernel = 1.0;

    // Iterate over all kernels
    for( int i = 0; i < numberOfSamples_; i++ )
    {
        marginalPdfOfCurrentKernel = 1.0;

        // Compute marginal pdf for current kernel
        for( unsigned int j = 0; j < marginalDimensions.size( ); j++ )
        {
            marginalPdfOfCurrentKernel *= kernelPointersMatrix_[ i ][ marginalDimensions[ j ] ]->evaluatePdf( independentVariables( j ) );
        }

        probabilityDensity += marginalPdfOfCurrentKernel;
    }

    // Average marginal pdf over all samples
    return probabilityDensity / static_cast< double >( numberOfSamples_ );
}

//! Function to evaluate marginal distribution density at single dimension.
double KernelDensityDistribution::evaluateMarginalProbabilityDensity(
        const int marginalDimension, const double independentVariable )
{
    double probabilityDensity = 0.0;

    // Compute pdf at independentVariable in marginalDimension, averaged over all samples
    for( int i = 0; i < numberOfSamples_; i++ )
    {
        probabilityDensity += kernelPointersMatrix_[ i ][ marginalDimension ]->evaluatePdf( independentVariable );
    }
    return probabilityDensity / static_cast< double >( numberOfSamples_ );
}

//! Function to evaluate cumulative conditional probability of marginal distribution at single dimension
double KernelDensityDistribution::evaluateCumulativeConditionalMarginalProbability(
        const std::vector< int >& conditionDimensions,
        const std::vector< double >& conditions,
        const int marginalDimension, const double independentVariable )
{
    if( std::find( conditionDimensions.begin( ), conditionDimensions.end( ), marginalDimension ) !=
            conditionDimensions.end( ) )
    {
        throw std::runtime_error( "Error when evaluatiing cumulative conditional marginal kernel density probability, repeated indices found" );
    }

    double marginalConditionalCdfOfCurrentKernel = 1.0;
    double normalizationFactor = 0.0;
    double marginalValue = 0.0;

    // Iterate over all kernels
    for( int i = 0; i < numberOfSamples_; i++ )
    {
        marginalConditionalCdfOfCurrentKernel = 1.0;

        for( unsigned int j = 0; j < conditionDimensions.size( ); j++ )
        {
            marginalConditionalCdfOfCurrentKernel *= kernelPointersMatrix_[ i ][ conditionDimensions[ j ] ]->evaluateCdf( conditions[ j ] );
        }

        // Current value of marginalConditionalCdfOfCurrentKernel is the marginal cdf at the given conditional
        normalizationFactor += marginalConditionalCdfOfCurrentKernel;

        // Compute cdf for current kernel at marginal dimension
        marginalConditionalCdfOfCurrentKernel *= kernelPointersMatrix_[ i ][ marginalDimension ]->evaluateCdf( independentVariable );
        marginalValue += marginalConditionalCdfOfCurrentKernel;
    }

    return marginalValue / normalizationFactor;
}

//! Function to evaluate conditional probability density of marginal distribution at single dimension
double KernelDensityDistribution::evaluateConditionalMarginalProbabilityDensity(
        const std::vector< int >& conditionDimensions,
        const std::vector< double >& conditions,
        const int marginalDimension, const double independentVariable )
{
    if( std::find( conditionDimensions.begin( ), conditionDimensions.end( ), marginalDimension ) !=
            conditionDimensions.end( ) )
    {
        throw std::runtime_error( "Error when evaluatiing conditional marginal kernel density probability, repeated indices found" );
    }

    double marginalConditionalPdfOfCurrentKernel = 1.0;
    double normalizationFactor = 0.0;
    double marginalValue = 0.0;

    // Iterate over all kernels
    for( int i = 0; i < numberOfSamples_; i++ )
    {
        marginalConditionalPdfOfCurrentKernel = 1.0;
        for( unsigned int j = 0; j < conditionDimensions.size( ); j++ )
        {
            marginalConditionalPdfOfCurrentKernel *=
                    kernelPointersMatrix_[ i ][ conditionDimensions[ j ] ]->evaluatePdf( conditions[ j ] );
        }

        // Current value of marginalConditionalPdfOfCurrentKernel is the marginal pdf at the given conditional
        normalizationFactor += marginalConditionalPdfOfCurrentKernel;

        // Compute pdf for current kernel at marginal dimension
        marginalConditionalPdfOfCurrentKernel *=
                kernelPointersMatrix_[ i ][ marginalDimension ]->evaluatePdf( independentVariable );
        marginalValue += marginalConditionalPdfOfCurrentKernel;
    }

    return marginalValue / normalizationFactor;
}

} // namespace statistics

} // namespace tudat

