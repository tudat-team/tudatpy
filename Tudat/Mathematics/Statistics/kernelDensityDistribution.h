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
    //! Default constructor
    EpanechnikovKernelDistribution(double Mean, double BandWidth){
        mean = Mean;
        bandWidth = BandWidth;
    }

    //! Get probability density
    double evaluatePdf( const double& x );

    //! Get cumulative probability
    double evaluateCdf(const double &x );

protected:

private:

    //! Mean of the kernel
    double mean;

    //! Bandwidth of the kernel
    double bandWidth;
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
private:

    //! Construct kernel
    boost::shared_ptr< ContinuousProbabilityDistribution< double > > constructKernel(double Mean, double BandWidth)
    {
        boost::shared_ptr< ContinuousProbabilityDistribution< double > > kernelPointer;

        if( kernelType == KernelType::Epanechnikov )
        {
            kernelPointer = boost::make_shared< EpanechnikovKernelDistribution >(Mean, BandWidth);
        }
        else
        { // Gaussian
            kernelPointer = createBoostRandomVariable(
                        normal_boost_distribution, { Mean, BandWidth } );
        }

        return kernelPointer;
    }

    //! Set kernel pointer matrix
    void generateKernelPointerMatrix();

    //! Compute sample mean
    void computeSampleMean();

    //! Compute sample mean
    void computeSampleVariance();

    //! Scale the samples such to the required variance
    void scaleSamplesWithVariance();

    //! Get median from datasamples
    Eigen::VectorXd getMedian( std::vector< Eigen::VectorXd > samples );

public:

    //! Constructor
    KernelDensityDistribution( std::vector< Eigen::VectorXd > samples,
                               double BandWidthFactor = 1.0,
                               KernelType kernel_type = KernelType::Gaussian,
                               bool adjustStandardDeviation = false,
                               Eigen::VectorXd standardDeviation = Eigen::VectorXd::Zero(2,1),
                               Eigen::VectorXd manualBandwidth = Eigen::VectorXd::Zero( 0 ) );

    //! Get probability density
    double evaluatePdf(const Eigen::VectorXd& x);

    //! Get cumulative probability ( F(x) = P(x<X) )
    double evaluateCdf( const Eigen::VectorXd& x );

    //! Get probability density of marginal distribution
    double getMarginalProbabilityDensity( std::vector< int > marginalDimensions, Eigen::VectorXd x );

    //! Get probability density of marginal distribution
    double getMarginalProbabilityDensityd( int marginalDimension, double x );

    //! Get cumulative probability of marginal distribution
    double getCumulativeMarginalProbability( int marginalDimension, double x );

    //! get cumulative probability of marginal distribution of conditional distribution
    double getCumulativeConditionalMarginalProbability( std::vector<int> conditionDimensions ,
                                                        std::vector<double> conditions,
                                                        int marginalDimension , double x );

    //! get probability density of marginal distribution of conditional distribution
    double getConditionalMarginalProbabilityDensity( std::vector<int> conditionDimensions ,
                                                        std::vector<double> conditions,
                                                        int marginalDimension , double x );

    //! Compute optimal bandwidth
    void computeOptimalBandWidth();

    //! Get optimal bandWidth
    Eigen::VectorXd getOptimalBandWidth(){
        return optimalBandwidth;
    }

    //! Get bandWidth
    Eigen::VectorXd getBandWidth(){
        return bandWidth;
    }

    //! Get samples
    std::vector< Eigen::VectorXd > getSamples(){
        return dataSamples;
    }

    //! Set bandWidth
    void setBandWidth( Eigen::VectorXd BandWidth ){
        bandWidth = BandWidth;
        generateKernelPointerMatrix(); // Reset kernels
    }

    //! Get sample mean.
    Eigen::VectorXd getSampleMean(){
        return sampleMean;
    }

    //! Get sample variance.
    Eigen::VectorXd getSampleVariance(){
        return sampleVariance;
    }

    //! Get sample standard deviation.
    Eigen::VectorXd getSampleStandardDeviation(){
        return sampleStandardDeviation;
    }

    //! Get number of samples
    int getNumberOfSamples(){
        return numberOfSamples;
    }

protected:

private:

    //! Datasamples
    std::vector< Eigen::VectorXd > dataSamples;

    //! sample mean
    Eigen::VectorXd sampleMean;

    //! sample variance
    Eigen::VectorXd sampleVariance;

    //! sample standard deviation
    Eigen::VectorXd sampleStandardDeviation;

    //! Standard deviation after correction
    Eigen::VectorXd standardDeviation_;

    //! Kernel type
    KernelType kernelType;

    //! BandWidth of Kernels
    Eigen::VectorXd bandWidth;

    //! Optimal bandWidth of Kernels.
    Eigen::VectorXd optimalBandwidth;

    //! Factor to multiply with optimal bandwidth.
    double bandWidthFactor;

    //! Number of dimensions of the samples and distribution.
    double dimensions;

    //! Number of datasamples
    double numberOfSamples;

    //! Matrix of 1D kernel pointers
    std::vector< std::vector< boost::shared_ptr< ContinuousProbabilityDistribution< double > > > > kernelPointersMatrix;

};

//! Pointer to Kernel Density distribution class
typedef boost::shared_ptr< KernelDensityDistribution > KernelDensityDistributionPointer;


} // Close Namespace statistics

} // Close Namespace tudat

#endif // TUDAT_KERNEL_DENSITY_DISTRIBUTION_H
