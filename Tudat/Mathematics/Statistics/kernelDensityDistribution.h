/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_KERNEL_DENSITY_DISTRIBUTION_H
#define TUDAT_KERNEL_DENSITY_DISTRIBUTION_H

#include <iostream>
#include <map>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Tudat/Mathematics/Statistics/continuousProbabilityDistributions.h>
#include <Tudat/Mathematics/Statistics/boostProbabilityDistributions.h>

namespace tudat
{

namespace statistics
{

//! Epanechnikov kernel distribution
/*!
 * Class that uses random samples to generate a distribution using Kernel Density Estimation.
 */
class EpanechnikovKernelDistribution: public tudat::statistics::ContinuousProbabilityDistribution< double >
{
public:

    //! Constructor
    EpanechnikovKernelDistribution( const double mean, const double bandWidth ):
        mean_( mean ), bandWidth_( bandWidth ){ }

    //! Get probability density
    double evaluatePdf( const double& independentVariable );

    //! Get cumulative probability
    double evaluateCdf( const double& independentVariable );

protected:

private:

    //! Mean of the kernel
    double mean_;

    //! Bandwidth of the kernel
    double bandWidth_;
};

//! Kernel type that is used in the Kernel Density Estimate
enum KernelType
{
    Gaussian = 0,
    Epanechnikov = 1
};

//! Kernel Density Estimate Distribution.
/*!
 * Class that uses random samples to generate a distribution using Kernel Density Estimation.
 */
class KernelDensityDistribution: public tudat::statistics::ContinuousProbabilityDistribution< Eigen::VectorXd >
{


public:

    //! Constructor
    KernelDensityDistribution(
            const std::vector< Eigen::VectorXd >& samples,
            const double bandWidthFactor = 1.0,
            const KernelType kernel_type = KernelType::Gaussian,
            const Eigen::VectorXd& manualStandardDeviation = Eigen::VectorXd::Zero( 0 ),
            const Eigen::VectorXd& manualBandwidth = Eigen::VectorXd::Zero( 0 ) );

    //! Get probability density
    double evaluatePdf( const Eigen::VectorXd& independentVariables );

    //! Get cumulative probability ( F(x) = P(x<X) )
    double evaluateCdf( const Eigen::VectorXd& independentVariables );

    //! Get probability density of marginal distribution
    double getMarginalProbabilityDensity(
            const std::vector< int >& marginalDimensions, const Eigen::VectorXd& x );

    //! Get probability density of marginal distribution
    double getMarginalProbabilityDensityd( const int marginalDimension, const double independentVariable );

    //! Get cumulative probability of marginal distribution
    double getCumulativeMarginalProbability( const int marginalDimension, const double independentVariable );

    //! get cumulative probability of marginal distribution of conditional distribution
    double getCumulativeConditionalMarginalProbability(
            const std::vector< int >& conditionDimensions,
            const std::vector< double >& conditions,
            const int marginalDimension, const double independentVariable );

    //! get probability density of marginal distribution of conditional distribution
    double getConditionalMarginalProbabilityDensity(
            const std::vector< int >& conditionDimensions,
            const std::vector< double >& conditions,
            const int marginalDimension, const double independentVariable );

    //! Compute optimal bandwidth
    void computeOptimalBandWidth( );

    //! Get optimal bandWidth_
    Eigen::VectorXd getOptimalBandWidth( )
    {
        return optimalBandwidth_;
    }

    //! Get bandWidth_
    Eigen::VectorXd getBandWidth( )
    {
        return bandWidth_;
    }

    //! Get samples
    std::vector< Eigen::VectorXd > getSamples( )
    {
        return dataSamples_;
    }

    //! Set bandWidth_
    void setBandWidth( const Eigen::VectorXd& bandWidth )
    {
        bandWidth_ = bandWidth;
        generateKernelPointerMatrix( ); // Reset kernels
    }

    //! Get sample mean.
    Eigen::VectorXd getSampleMean( )
    {
        return sampleMean_;
    }

    //! Get sample variance.
    Eigen::VectorXd getSampleVariance( )
    {
        return sampleVariance_;
    }

    //! Get sample standard deviation.
    Eigen::VectorXd getSampleStandardDeviation( ){

        return sampleStandardDeviation_;
    }

    //! Get number of samples
    int getNumberOfSamples( )
    {
        return numberOfSamples_;
    }

protected:

private:

    //! Construct kernel
    boost::shared_ptr< ContinuousProbabilityDistribution< double > > constructKernel(
            const double mean, const double bandWidth)
    {
        boost::shared_ptr< ContinuousProbabilityDistribution< double > > kernelPointer;

        if( kernelType_ == KernelType::Epanechnikov )
        {
            kernelPointer = boost::make_shared< EpanechnikovKernelDistribution >( mean, bandWidth );
        }
        else
        {
            kernelPointer = createBoostRandomVariable(
                        normal_boost_distribution, { mean, bandWidth } );
        }

        return kernelPointer;
    }

    //! Set kernel pointer matrix
    void generateKernelPointerMatrix( );

    //! Compute sample mean
    void computeSampleMean( );

    //! Compute sample mean
    void computeSampleVariance( );

    //! Scale the samples such to the required variance
    void scaleSamplesWithVariance( const Eigen::VectorXd& standardDeviation );

    //! Get median from datasamples
    Eigen::VectorXd getMedian( std::vector< Eigen::VectorXd > samples );

    //! Datasamples
    std::vector< Eigen::VectorXd > dataSamples_;

    //! sample mean
    Eigen::VectorXd sampleMean_;

    //! sample variance
    Eigen::VectorXd sampleVariance_;

    //! sample standard deviation
    Eigen::VectorXd sampleStandardDeviation_;

    //! Kernel type
    KernelType kernelType_;

    //! BandWidth of Kernels
    Eigen::VectorXd bandWidth_;

    //! Optimal bandWidth of Kernels.
    Eigen::VectorXd optimalBandwidth_;

    //! Number of dimensions of the samples and distribution.
    int dimensions_;

    //! Number of datasamples
    int numberOfSamples_;

    //! Matrix of 1D kernel pointers
    std::vector< std::vector< boost::shared_ptr< ContinuousProbabilityDistribution< double > > > > kernelPointersMatrix_;

};

//! Pointer to Kernel Density distribution class
typedef boost::shared_ptr< KernelDensityDistribution > KernelDensityDistributionPointer;


} // Close Namespace statistics

} // Close Namespace tudat

#endif // TUDAT_KERNEL_DENSITY_DISTRIBUTION_H
