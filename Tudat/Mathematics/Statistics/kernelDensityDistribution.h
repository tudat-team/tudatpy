/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_KERNELDENSITYDISTRIBUTION_H
#define TUDAT_KERNELDENSITYDISTRIBUTION_H

#include <map>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/Statistics/continuousProbabilityDistributions.h"
#include "Tudat/Mathematics/Statistics/boostProbabilityDistributions.h"
namespace tudat
{

namespace statistics
{

//! Class for Probability distribution using a single Epanechnikov kernel
class EpanechnikovKernelDistribution: public tudat::statistics::ContinuousProbabilityDistribution< double >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param mean Kernel mean
     * \param bandWidth Kernel bandwidth
     */
    EpanechnikovKernelDistribution( const double mean, const double bandWidth ):
        mean_( mean ), bandWidth_( bandWidth ){ }

    //! Function to evaluate pdf of distribution
    /*!
     *  Function to evaluate probability distribution function at given independentVariable value.
     *  \param independentVariable Value of independent variable
     *  \return Evaluated pdf
     */
    double evaluatePdf( const double& independentVariable );

    //! Function to evaluate cdf of distribution
    /*!
     *  Function to evaluate cumulative distribution function at given independentVariable value.
     *  \param independentVariable Value of independent variable
     *  \return Evaluated cdf
     */
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
    gaussian_kernel = 0,
    epanechnikov_kernel = 1
};

//! Class that uses random samples to generate a multivariate probability distribution using Kernel Density distribution.
/*!
 *  Class that uses random samples to generate a multivariate probability distribution using Kernel Density distribution.
 *  The bandwidth of the kernels may be supplied by the user, or an optimal distrubution may be computed by this class
 *  At present, the user has the choice of a Gaussian or Epanechnikov distribution for the kernels.
 */
class KernelDensityDistribution: public tudat::statistics::ContinuousProbabilityDistribution< Eigen::VectorXd >
{


public:

    //! Constructor
    /*!
     * Constructor
     * \param samples Data samples from which the empirical distribution is to be created.
     * \param bandWidthFactor Factor by which to multiply the computed bandwidth values of the kernels.
     * \param kernel_type Type of distribution to be used for the kernels in the distribution.
     * \param manualStandardDeviation Standard deviations with which to scale the entries of each sample. By default this
     * vector is empty and not used. If it is provided, it is used to scale the entry i of each sample by
     * manualStandardDeviation( i ) / sampleStandardDeviation( i ), where sampleStandardDeviation is the standard
     * deviation computed from the samples/
     * \param manualBandwidth Vector of bandwidths for each dimension that is to be used in the kernels. By default this
     * vector is empty and not used. Optimal bandwidths (scaled by bandWidthFactor) are used in this default case.
     */
    KernelDensityDistribution(
            const std::vector< Eigen::VectorXd >& samples,
            const double bandWidthFactor = 1.0,
            const KernelType kernel_type = KernelType::gaussian_kernel,
            const Eigen::VectorXd& manualStandardDeviation = Eigen::VectorXd::Zero( 0 ),
            const Eigen::VectorXd& manualBandwidth = Eigen::VectorXd::Zero( 0 ) );

    //! Function to evaluate pdf of distribution
    /*!
     *  Function to evaluate probability distribution function at given independentVariable value.
     *  \param independentVariables Values of independent variable
     *  \return Evaluated pdf
     */
    double evaluatePdf( const Eigen::VectorXd& independentVariables );

    //! Function to evaluate cdf of distribution
    /*!
     *  Function to evaluate cumulative distribution function at given independentVariable value.
     *  \param independentVariables Values of independent variable
     *  \return Evaluated cdf
     */
    double evaluateCdf( const Eigen::VectorXd& independentVariables );

    //! Function to evaluate probability density of marginal (in one or more dimensions) distribution.
    /*!
     * Function to evaluate probability density of marginal (in one or more dimensions) distribution,
     * \param marginalDimensions Dimensions over which the marginal probaility is to be computed
     * \param independentVariables Values of independent variables (in marginalDimensions) at which marginal
     * probability is to be computed
     * \return Probability density of marginal distribution.
     */
    double evaluateMarginalProbabilityDensity(
            const std::vector< int >& marginalDimensions, const Eigen::VectorXd& independentVariables );

    //! Function to evaluate marginal distribution density at single dimension.
    /*!
     * Function to evaluate marginal distribution density at single dimension.
     *  \param marginalDimension Dimension in which the marginal distribution is to be computed
     *  \param independentVariable Value of independent variable in dimension marginalDimension
     *  \return Probability density of marginal distribution at given dimension and independent variable.
     */
    double evaluateMarginalProbabilityDensity( const int marginalDimension, const double independentVariable );

    //! Function to evaluate cumulative probability of marginal distribution at single dimension.
    /*!
     *  Function to evaluate cumulative probability of marginal distribution at single dimension.
     *  \param marginalDimension Dimension in whihc the marginal distribution is to be computed
     *  \param independentVariable Value of independent variable in dimension marginalDimension
     *  \return cumulative probability of marginal distribution at given dimension and independent variable.
     */
    double evaluateCumulativeMarginalProbability( const int marginalDimension, const double independentVariable );

    //! Function to evaluate cumulative probability of marginal distribution of conditional distribution at single dimension
    /*!
     * Function to evaluate cumulative probability of marginal distribution of conditional distribution at single dimension,
     * given conditions of cdf on other dimensions
     * \param conditionDimensions List of dimensions in which cdf is given
     * \param conditions Independent variables at which values of cdf, in dimensions listed in conditionDimensions,
     * are to be computed.
     * \param marginalDimension Dimension in which the marginal cdf is to be computed
     * \param independentVariable Independent variable at which conditional marginal cdf is to be computed
     * \return Cumulative conditional probability of marginal distribution at given settings.
     */
    double evaluateCumulativeConditionalMarginalProbability(
            const std::vector< int >& conditionDimensions,
            const std::vector< double >& conditions,
            const int marginalDimension, const double independentVariable );

    //! Function to evaluate probability density of marginal distribution of conditional distribution at single dimension
    /*!
     * Function to evaluate probability density of marginal distribution of conditional distribution at single dimension,
     * given conditions of cdf on other dimensions
     * \param conditionDimensions List of dimensions in which pdf is given
     * \param conditions Independent variables at which values of pdf, in dimensions listed in conditionDimensions,
     * are to be computed.
     * \param marginalDimension Dimension in which the marginal pdf is to be computed
     * \param independentVariable Independent variable at which conditional marginal pdf is to be computed
     * \return Conditional pdf of marginal distribution at given settings.
     */
    double evaluateConditionalMarginalProbabilityDensity(
            const std::vector< int >& conditionDimensions,
            const std::vector< double >& conditions,
            const int marginalDimension, const double independentVariable );

    //! Compute optimal bandwidth using the provided samples.
    void computeOptimalBandWidth( );

    //! Function to retrieve the optimal bandwidth
    /*!
     * Function to retrieve the optimal bandwidth
     * \return Optimal bandwidth
     */
    Eigen::VectorXd getOptimalBandWidth( )
    {
        return optimalBandwidth_;
    }

    //! Function to retrieve the bandwidth
    /*!
     * Function to retrieve the bandwidth
     * \return Bandwidth
     */
    Eigen::VectorXd getBandWidth( )
    {
        return bandWidth_;
    }

    //! Function to retrieve the data samples used to generate kernel density distribution.
    /*!
     * Function to retrieve the data samples used to generate kernel density distribution.
     * \return Data samples used to generate kernel density distribution.
     */
    std::vector< Eigen::VectorXd > getSamples( )
    {
        return dataSamples_;
    }

    //! Function to manually reset the bandwidth
    /*!
     * Function to manually reset the bandwidth. Calling this function automatically recomputes the kernels.
     * \param bandWidth New bandwidth to be used for each dimension.
     */
    void setBandWidth( const Eigen::VectorXd& bandWidth )
    {
        bandWidth_ = bandWidth;
        generateKernelPointerMatrix( ); // Reset kernels
    }

    //! Function to retrieve the sample mean.
    /*!
     * Function to retrieve the sample mean.
     * \return Sample mean.
     */
    Eigen::VectorXd getSampleMean( )
    {
        return sampleMean_;
    }

    //! Function to retrieve the sample variance.
    /*!
     * Function to retrieve the sample variance.
     * \return Sample variance.
     */
    Eigen::VectorXd getSampleVariance( )
    {
        return sampleVariance_;
    }

    //! Function to retrieve the sample standard deviation.
    /*!
     * Function to retrieve the sample standard deviation.
     * \return Sample standard deviation.
     */
    Eigen::VectorXd getSampleStandardDeviation( )
    {

        return sampleStandardDeviation_;
    }

    //! Function to retrieve the number of samples.
    /*!
     * Function to retrieve the number of samples.
     * \return Number of samples.
     */
    int getNumberOfSamples( )
    {
        return numberOfSamples_;
    }

protected:

private:

    //! Function that constructs a single kernel.
    /*!
     * Function that constructs a single kernel, called iteratively to construct the full multivariate kernel density.
     * \param mean Mean value of kernel that is to be created.
     * \param bandWidth Bandwidth of kernel thatr is to be created.
     * \return Probability distribution for single kernel.
     */
    boost::shared_ptr< ContinuousProbabilityDistribution< double > > constructKernel(
            const double mean, const double bandWidth )
    {
        boost::shared_ptr< ContinuousProbabilityDistribution< double > > kernelPointer;

        // Check kernel type and create selected kernel.
        if( kernelType_ == KernelType::epanechnikov_kernel )
        {
            kernelPointer = boost::make_shared< EpanechnikovKernelDistribution >( mean, bandWidth );
        }
        else if( kernelType_ == KernelType::gaussian_kernel )
        {
            kernelPointer = createBoostRandomVariable(
                        normal_boost_distribution, { mean, bandWidth } );
        }
        else
        {
            throw std::runtime_error( "Error when constructing probability kernels, kernel type not recognized" );
        }

        return kernelPointer;
    }

    //! Function that generates the kernel density distribution based on the samples and kernel type that is provided
    /*!
     *  Function that generates the kernel density distribution based on the samples and kernel type that is provided.
     *  This function iteratively calls the constructKernel function to create a kernel for each entry in each sample.
     */
    void generateKernelPointerMatrix( );

    //! Function that computes and sets the sample mean.
    /*!
     *  Function that computes and sets the sample mean (sampleMean_ variable).
     */
    void computeSampleMean( );

    //! Function that computes and sets the sample variance and standard deviation.
    /*!
     *  Function that computes and sets the sample variance and standard deviation
     *  (sampleVariance_ and sampleStandardDeviation_ variable).
     */
    void computeSampleVariance( );

    //! Function to scale the samples to the required standard deviation
    /*!
     *  Function to scale the samples to the required standard deviation
     *  \param standardDeviation New standard deviation that the data is to have.
     */
    void scaleSamplesWithVariance( const Eigen::VectorXd& standardDeviation );

    //! Function to retrieve the sample set median
    /*!
     * Function to retrieve the sample set median
     * \param samples Set of samples from which the (per entry) median is to be computed.
     * \return Sample set median
     */
    Eigen::VectorXd getMedian( const std::vector< Eigen::VectorXd >& samples );

    //! Datasamples
    std::vector< Eigen::VectorXd > dataSamples_;

    //! Sample mean (per entry)
    Eigen::VectorXd sampleMean_;

    //! Sample variance (per entry)
    Eigen::VectorXd sampleVariance_;

    //! Sample standard deviation (per entry)
    Eigen::VectorXd sampleStandardDeviation_;

    //! Kernel type
    KernelType kernelType_;

    //! Bandwidth of Kernels
    Eigen::VectorXd bandWidth_;

    //! Optimal bandwidth of Kernels.
    Eigen::VectorXd optimalBandwidth_;

    //! Number of dimensions of the samples and distribution.
    int dimensions_;

    //! Number of datasamples
    int numberOfSamples_;

    //! Matrix (vector of vectors) of 1D kernel pointers that define full kernel density distribution.
    std::vector< std::vector< boost::shared_ptr< ContinuousProbabilityDistribution< double > > > > kernelPointersMatrix_;

};

//! Pointer to Kernel Density distribution class
typedef boost::shared_ptr< KernelDensityDistribution > KernelDensityDistributionPointer;


} // namespace statistics

} // namespace tudat

#endif // TUDAT_KERNELDENSITYDISTRIBUTION_H
