/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      160429    R. Hoogendoorn    File created.
 *
 *    References
 *
 *    Notes
 *
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
double EpanechnikovKernelDistribution::evaluatePdf(const double &x )
{
    if( (x-mean)>=(-bandWidth) && (x-mean)<=(bandWidth) )
    {
        return ( 3.0/(4.0*bandWidth) )*( 1.0 - std::pow( (x-mean) / bandWidth , 2.0 ) ) ;
    }
    else
    {
        return 0.0 ;
    }
}

//! Get probability mass
double EpanechnikovKernelDistribution::evaluateCdf( const double& x )
{
    if( (x-mean)>=(-bandWidth) && (x-mean)<=(bandWidth) )
    {
        return ( 3.0/(4.0*bandWidth) )*( (x-mean) - ( std::pow( (x-mean) , 3.0 )
                                                      / (3.0*std::pow( bandWidth , 2.0 )) ) ) + 0.5  ;
    }
    else if( (x-mean)<(-bandWidth) )
    {
        return 0.0 ;
    }
    else
    {
        return 1.0 ;
    }

}

//! Constructor
KernelDensityDistribution::KernelDensityDistribution(
        std::vector< Eigen::VectorXd > samples,
        double BandWidthFactor,
        KernelType kernel_type, bool adjustStandardDeviation,
        Eigen::VectorXd standardDeviation,
        Eigen::VectorXd manualBandwidth )
{
    // Load data
    dataSamples = samples;
    dimensions = static_cast< double >( samples[0].rows() );
    numberOfSamples = static_cast< double >( samples.size() );
    standardDeviation_ = standardDeviation;

    // Compute datasample properties
    computeSampleMean();
    computeSampleVariance();

    // Adjust samples with standard deviation and recompute datasample properties
    if( adjustStandardDeviation == true )
    {
        scaleSamplesWithVariance();
        computeSampleMean();
        computeSampleVariance();
    }


    // Compute bandwidth
    bandWidthFactor = BandWidthFactor ;

    computeOptimalBandWidth();

    bandWidth = optimalBandwidth * bandWidthFactor;

    // Construct kernel matrix rows: samples , cols: dimensions
    kernelType = kernel_type;

    if( manualBandwidth.rows( ) > 0 )
    {
        bandWidth = manualBandwidth;
    }
    generateKernelPointerMatrix();
}

//! Set kernel pointer matrix
void KernelDensityDistribution::generateKernelPointerMatrix()
{
    std::vector< std::vector< boost::shared_ptr< ContinuousProbabilityDistribution< double > > > > matrix(0);
    kernelPointersMatrix = matrix;

    for( unsigned int i = 0; i < bandWidth.rows( ); i++ )
    {
        if( bandWidth( i ) < 10.0 * std::numeric_limits< double >::epsilon( ) )
        {
            std::string errorMessage = "Error in kernel density distribution, index " +
                    boost::lexical_cast< std::string >( i ) + "has bandwidth " +
                    boost::lexical_cast< std::string >( bandWidth( i ) );
            throw std::runtime_error( errorMessage );
        }
    }
    // Fill kernel pointer matrix with distribution pointer objects
    std::vector< boost::shared_ptr< ContinuousProbabilityDistribution< double > > > vector( dataSamples[0].rows() );
    for( unsigned int i = 0 ; i < dataSamples.size() ; i++ )
    {
        for( unsigned int j = 0 ; j < vector.size() ; j++)
        {
            vector[j] = constructKernel( dataSamples[i](j) , bandWidth(j) );
        }
        kernelPointersMatrix.push_back( vector );
    }
}

void KernelDensityDistribution::computeSampleMean()
{
    sampleMean = Eigen::VectorXd::Zero( static_cast< int >(dimensions) , 1 ) ;
    for( unsigned int i = 0 ; i < dataSamples.size() ; i++ )
    {
        sampleMean +=  dataSamples[i] ;
    }
    sampleMean = sampleMean / numberOfSamples ;
}

void KernelDensityDistribution::computeSampleVariance()
{
    sampleVariance = Eigen::VectorXd::Zero( static_cast< int >(dimensions) ) ;
    for( unsigned int i = 0 ; i < dataSamples.size() ; i++ )
    {
        sampleVariance += ( dataSamples[i] - sampleMean ).cwiseAbs2() ;
    }
    sampleVariance = sampleVariance / ( numberOfSamples - 1.0 ) ;
    sampleStandardDeviation = sampleVariance.cwiseSqrt();
}

//! Scale the samples such to the required variance
void KernelDensityDistribution::scaleSamplesWithVariance()
{
    for( unsigned int i = 0 ; i < dataSamples.size() ; i++ )
    {
        for( int j = 0 ; j < dataSamples[i].rows() ; j++ )
        {
            dataSamples[i](j) = dataSamples[i](j) / sampleStandardDeviation(j) * standardDeviation_(j) ;
        }
    }
}

//! Compute the optimal bandwidth
void KernelDensityDistribution::computeOptimalBandWidth()
{
    int dims = static_cast< int >(dimensions) ;

    optimalBandwidth = Eigen::VectorXd::Zero( dims , 1);

    // Calculate sigma (median absolute deviation estimator)
    Eigen::VectorXd medianOfSamples = getMedian( dataSamples );

    std::vector< Eigen::VectorXd > dataSamples2( dataSamples.size() );
    for( unsigned int i = 0 ; i < dataSamples.size() ; i++)
    {
        dataSamples2[i] = dataSamples[i] - medianOfSamples ;
        dataSamples2[i] = dataSamples2[i].cwiseAbs();
    }
    Eigen::VectorXd medianOfSamples2 = getMedian( dataSamples2 );
    Eigen::VectorXd sigma = medianOfSamples2 / 0.6745;

    // Calculate optimal Bandwidth
    optimalBandwidth = std::pow( 4.0 / (( dimensions + 2.0 )*numberOfSamples) , 1.0/( dimensions + 4.0 ) ) *
            sigma;
}

//! Get median from datasamples
Eigen::VectorXd KernelDensityDistribution::getMedian( std::vector< Eigen::VectorXd > samples )
{
    int dims = samples[0].rows() ;
    std::vector< std::vector< double > > dataSamplesVectors( dims );

    std::vector< double > data(0);
    for( int i = 0 ; i < dims ; i++ )
    {
        for( unsigned int j = 0 ; j < samples.size() ; j++)
        {
            data.push_back( samples[j](i) );
        }
        dataSamplesVectors[i] = data ;
        data.clear();
    }

    Eigen::VectorXd medianOfSamples( dims );
    for( int i = 0 ; i < dims ; i++)
    {
        medianOfSamples[i] = tudat::statistics::computeSampleMedian( dataSamplesVectors[i] );
    }
    return medianOfSamples;
}

//! Get probability density of the kernel density distribution
double KernelDensityDistribution::evaluatePdf( const Eigen::VectorXd& x )
{
    double propbabilityDensity = 0.0 ;
    double value = 1.0 ;
    for( int i = 0 ; i < numberOfSamples ; i++ )
    {
        for( int dim = 0 ; dim < dimensions ; dim++ )
        {
            value = value * kernelPointersMatrix[ i ][ dim ]->evaluatePdf( x( dim ) ) ;
        }
        propbabilityDensity += value ;
        value = 1.0 ;
    }
    return propbabilityDensity / numberOfSamples ;
}

//! Get cumulative probability of the kernel density distribution
double KernelDensityDistribution::evaluateCdf( const Eigen::VectorXd& x )
{
    double cumulativeProbability = 0.0 ;
    double value = 1.0 ;
    for( int i = 0 ; i < numberOfSamples ; i++ )
    {
        for( int dim = 0 ; dim < dimensions ; dim++ )
        {
            value = value * kernelPointersMatrix[ i ][ dim ]->evaluateCdf( x( dim ) ) ;
        }
        cumulativeProbability += value ;
        value = 1.0 ;
    }
    return cumulativeProbability / numberOfSamples ;
}

//! Get cumulative probability of marginal distribution
double KernelDensityDistribution::getCumulativeMarginalProbability(int marginalDimension, double x )
{
    double cumulativeProbability = 0.0 ;
    for( int i = 0 ; i < numberOfSamples ; i++ )
    {
        cumulativeProbability += kernelPointersMatrix[ i ][ marginalDimension ]->evaluateCdf( x ) ;
    }
    return cumulativeProbability / numberOfSamples ;
}

//! Get probability density of marginal distribution
double KernelDensityDistribution::getMarginalProbabilityDensity( std::vector< int > marginalDimensions, Eigen::VectorXd x )
{
    double probabilityDensity = 0.0 ;
    double value = 1.0 ;
    for( int i = 0 ; i < numberOfSamples ; i++ )
    {
        for( unsigned int j = 0 ; j < marginalDimensions.size() ; j++ )
        {
            value = value * kernelPointersMatrix[ i ][ marginalDimensions[j] ]->evaluatePdf( x(j) ) ;
        }
        probabilityDensity += value ;
        value = 1.0 ;
    }
    return probabilityDensity / numberOfSamples ;
}

//! Get probability density of marginal distribution
double KernelDensityDistribution::getMarginalProbabilityDensityd( int marginalDimension, double x )
{
    double probabilityDensity = 0.0 ;
    for( int i = 0 ; i < numberOfSamples ; i++ )
    {
        probabilityDensity += kernelPointersMatrix[ i ][ marginalDimension ]->evaluatePdf( x ) ;
    }
    return probabilityDensity / numberOfSamples ;
}

//! get cumulative probability of marginal distribution of conditional distribution
double KernelDensityDistribution::getCumulativeConditionalMarginalProbability( std::vector< int > conditionDimensions ,
                                                    std::vector< double > conditions,
                                                    int marginalDimension , double x )
{
    double value = 1.0 ;
    double normalizationFactor = 0.0 ;
    double marginalValue = 0.0 ;

    // multiply all conditions * marginal , divide marginals of conditions
    for( int i = 0 ; i < numberOfSamples ; i++ )
    {
        value = 1.0 ;

        for( unsigned int j = 0 ; j < conditionDimensions.size() ; j++)
        {
            value = value * kernelPointersMatrix[ i ][ conditionDimensions[j] ]->evaluateCdf( conditions[j] ) ;
        }
        normalizationFactor += value ;

        value = value * kernelPointersMatrix[ i ][ marginalDimension ]->evaluateCdf( x ) ;
        marginalValue += value ;
    }

    return marginalValue / ( normalizationFactor ) ;
}

//! get probability density of marginal distribution of conditional distribution
double KernelDensityDistribution::getConditionalMarginalProbabilityDensity( std::vector< int > conditionDimensions ,
                                                    std::vector< double > conditions,
                                                    int marginalDimension , double x )
{
    double value = 1.0 ;
    double normalizationFactor = 0.0 ;
    double marginalValue = 0.0 ;

    // multiply all conditions * marginal , divide marginals of conditions
    for( int i = 0 ; i < numberOfSamples ; i++ )
    {
        value = 1.0 ;

        for( unsigned int j = 0 ; j < conditionDimensions.size() ; j++)
        {
            value = value * kernelPointersMatrix[ i ][ conditionDimensions[j] ]->evaluatePdf( conditions[j] ) ;
        }
        normalizationFactor += value ;

        value = value * kernelPointersMatrix[ i ][ marginalDimension ]->evaluatePdf( x ) ;
        marginalValue += value ;
    }

    return marginalValue / ( normalizationFactor ) ;
}

} // Close Namespace statistics

} // Close Namespace tudat

