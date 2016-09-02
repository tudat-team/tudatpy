/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      160902    R. Hoogendoorn        File created.

 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <iostream> // cout sometimes needs this
#include <vector>
#include <limits>

#include <Eigen/Core>

#include <boost/test/unit_test.hpp>
#include "tudat/Basics/testMacros.h"

#include "tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "tudat/InputOutput/basicInputOutput.h"
#include "tudat/Mathematics/Statistics/kernelDensityDistribution.h"

// Pseudo random sampler
#include <boost/random.hpp> // Random generator

namespace tudat
{
namespace unit_tests
{

//! Creates a std::vector of linear spaced values
std::vector<double> linspace(double start, double end, int N){

    std::vector <double> x(0) ;

    if ( N > 1 )
    {
        for( int i = 0 ; i < N ; i++ )
        {
            x.push_back( start + ( end - start )*i /( N - 1.0 ) ) ;
        }
    }
    return x;
}

//! Generator random vector using pseudo random generator
std::vector<Eigen::VectorXd> generateRandomVectorUniform(int seed, int numberOfSamples,
                                             Eigen::VectorXd lowerBound, Eigen::VectorXd upperBound)
{
    // Compute properties
    Eigen::VectorXd width = upperBound - lowerBound;
    Eigen::VectorXd average = (upperBound + lowerBound)/2.0 ;

    // Setup Random generator
    typedef boost::mt19937 RandomGeneratorType; // Mersenne Twister
    RandomGeneratorType randomGenerator(seed);              // Create random generator

    boost::uniform_real<> uniformDistribution(0.0 , 1.0); //
    boost::variate_generator< RandomGeneratorType, boost::uniform_real<> >
            Dice(randomGenerator, uniformDistribution); // define random generator

    std::vector< Eigen::VectorXd > randomSamples(numberOfSamples) ;
    Eigen::VectorXd randomSample( lowerBound.rows() ) ;

    // Sample
    for(int i = 0 ; i < numberOfSamples ; i++ ){ // Generate N samples
        for(int j = 0 ; j < randomSample.rows() ; j++){ // Generate vector of samples
            randomSample(j) = Dice() - 0.5 ;
        }
        randomSamples[i] = randomSample.cwiseProduct(width) + average ;
    }
    return randomSamples;
}

BOOST_AUTO_TEST_SUITE( test_Kernel_Density_Distribution )

//using tudat::mathematical_constants::PI;

BOOST_AUTO_TEST_CASE( test_Optimal_Bandwidth )
{

    using namespace tudat::statistics;

    std::vector< Eigen::VectorXd > data(0);
    Eigen::VectorXd sample(2);
    sample << 2.2538 , 1.177 ;
    data.push_back( sample ) ;
    sample << 0.76529 , 0.356 ;
    data.push_back( sample ) ;
    sample << 1.5179 , 1.14 ;
    data.push_back( sample ) ;
    sample << 2.0972 , 0.34093 ;
    data.push_back( sample ) ;
    sample << 2.6727 , 1.301 ;
    data.push_back( sample ) ;
    sample << 2.8779 , 0.48998 ;
    data.push_back( sample ) ;
    sample << 1.6416 , 0.27523 ;
    data.push_back( sample ) ;
    sample << 0.41587 , 0.35152 ;
    data.push_back( sample ) ;
    sample << 0.44788 , 0.86246 ;
    data.push_back( sample ) ;
    sample << 0.77252 , 0.6626 ;
    data.push_back( sample ) ;

    KernelDensityDistribution distribution( data );

    // Generated expected Bandwidth using MATLAB
    Eigen::VectorXd expectedBandwidth(2);
    expectedBandwidth << 0.819018940571854 , 0.263389630191125 ;

    Eigen::VectorXd computedOptimalBandwidth(2);
    computedOptimalBandwidth = distribution.getOptimalBandWidth() ;

    BOOST_CHECK_CLOSE_FRACTION( expectedBandwidth(0) , computedOptimalBandwidth(0) , 2E-5 );
    BOOST_CHECK_CLOSE_FRACTION( expectedBandwidth(1) , computedOptimalBandwidth(1) , 5E-6 );

}

BOOST_AUTO_TEST_CASE( test_Epanechnikov_Kernel )
{
    tudat::statistics::EpanechnikovKernelDistribution distribution(1.0, 2.0);

    BOOST_CHECK_CLOSE_FRACTION( distribution.getCumulativeProbability( 5.0 ) , 1.0 , 1E-15 );
    BOOST_CHECK_CLOSE_FRACTION( distribution.getCumulativeProbability( 3.0 ) , 1.0 , 1E-15 );
    BOOST_CHECK_CLOSE_FRACTION( distribution.getCumulativeProbability( 1.0 ) , 0.5 , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( distribution.getCumulativeProbability( -1.0 ) - 0.0 ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( distribution.getCumulativeProbability( -2.0 ) - 0.0 ) , 1E-15 );
    // Computed using MATLAB
    BOOST_CHECK_SMALL( std::fabs( distribution.getCumulativeProbability( 0.0 ) - 0.15625 ) , 1E-15 );

    BOOST_CHECK_SMALL( std::fabs( distribution.getProbabilityDensity( -1.0 ) - 0.0 ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( distribution.getProbabilityDensity( 3.0 ) - 0.0 ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( distribution.getProbabilityDensity( 4.0 ) - 0.0 ) , 1E-15 );
    BOOST_CHECK_SMALL( std::fabs( distribution.getProbabilityDensity( -3.0 ) - 0.0 ) , 1E-15 );

    // Computed using MATLAB
    BOOST_CHECK_SMALL( std::fabs( distribution.getProbabilityDensity( 0.0 ) - 0.28125 ) , 1E-15 );
}

BOOST_AUTO_TEST_CASE( test_Probability_Density_1D )
{
    using namespace tudat::statistics;

    Eigen::VectorXd sample(1);
    std::vector< Eigen::VectorXd > samples(0);

    sample(0) = 3.0 ;
    samples.push_back( sample );
    sample(0) = 1.0 ;
    samples.push_back( sample );
    sample(0) = 4.0 ;
    samples.push_back( sample );
    sample(0) = 7.0 ;
    samples.push_back( sample );

    KernelDensityDistribution distribution(samples , 1.0 , KernelType::Gaussian ) ;

    Eigen::VectorXd location(1);
    // Case 1
    {
        location << 5.396093156208099e-01 ;

        double computedDensity = distribution.getProbabilityDensity( location );

        // Data obtained from MATLAB ksdensity function
        double expectedDensity = 8.426926493705324e-02 ;

        // Check that the same bandwidths are used compared with MATLAB
        Eigen::VectorXd expectedBandwidth(1);
        expectedBandwidth << 1.785192502061299 ;
        Eigen::VectorXd computedBandwidth = distribution.getBandWidth() ;
        Eigen::VectorXd difference = expectedBandwidth - computedBandwidth ;
        BOOST_CHECK_SMALL( std::fabs( difference(0) ) , 1E-15 );

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );
    }

    // Case 2
    {
        Eigen::VectorXd location(1);
        location << 3.915600227210264 ;

        double computedDensity = distribution.getProbabilityDensity( location );
        double expectedDensity = 1.320686749820565e-01 ;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );
    }

    // Case 3
    {
        Eigen::VectorXd location(1);
        location << 8.135588866697081 ;

        double computedDensity = distribution.getProbabilityDensity( location );
        double expectedDensity = 5.036314749096298e-02 ;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );
    }
}

BOOST_AUTO_TEST_CASE( test_Probability_Density_2D )
{
    using namespace tudat::statistics;

    Eigen::VectorXd sample(2);
    std::vector< Eigen::VectorXd > samples(0);

    sample << 3.0 , 2E-2 ;
    samples.push_back( sample );
    sample << -1.0 , 5E-2 ;
    samples.push_back( sample );
    sample << 5.0 , 2E-1 ;
    samples.push_back( sample );
    sample << 2.0 , 1E-3 ;
    samples.push_back( sample );

    KernelDensityDistribution distribution(samples , 1.0 , KernelType::Gaussian ) ;

    // Case 1
    {
        Eigen::VectorXd location(2);
        location << 1.141869732328655 , 2.995235933344523e-02 ;

        double computedDensity = distribution.getProbabilityDensity( location );
        double expectedDensity = 1.136916246361795 ;

        // Check that the same bandwidths are used compared with MATLAB
        Eigen::VectorXd expectedBandwidth(2);
        expectedBandwidth << 1.765086418052112 , 2.882974482818450e-02 ;
        Eigen::VectorXd computedBandwidth = distribution.getBandWidth() ;
        Eigen::VectorXd difference = expectedBandwidth - computedBandwidth ;
        BOOST_CHECK_SMALL( std::fabs( difference(0) ) , 1E-15 );
        BOOST_CHECK_SMALL( std::fabs( difference(1) ) , 1E-15 );

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );
    }

    // Case 2
    {
        Eigen::VectorXd location(2);
        location << -2.862738183470956 , 1.582207969089993e-01 ;

        double computedDensity = distribution.getProbabilityDensity( location );
        double expectedDensity = 4.038687551312721e-04 ;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );
    }

    // Case 3
    {
        Eigen::VectorXd location(2);
        location << 6.862738183470957 , 2.995235933344523e-02 ;

        double computedDensity = distribution.getProbabilityDensity( location );
        double expectedDensity = 7.784149895640656e-02 ;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );
    }

    // Marginal
    {
        // Marginal CDF compare with MATLAB
        // Bandwidth slightly different in 1D
        Eigen::VectorXd bandwidth(2);
        bandwidth << 1.785192502061299 , 2.915814420033455e-02 ;
        distribution.setBandWidth( bandwidth );

        BOOST_CHECK_CLOSE_FRACTION(
                    distribution.getCumulativeMarginalProbability(0 , 2.276047714155371E-1 ),
                    2.446324562687310E-1 , 3E-3 ) ;
        BOOST_CHECK_CLOSE_FRACTION(
                    distribution.getCumulativeMarginalProbability(1 , 6.083875672099921e-02 ),
                    6.360523509878839e-01 , 2E-3 ) ;
        BOOST_CHECK_CLOSE_FRACTION(
                    distribution.getCumulativeMarginalProbability(1 , 9.861136936766661e-02 ),
                    7.371490745984073e-01 , 8E-4 ) ;

        BOOST_CHECK_SMALL( std::fabs(
            distribution.getCumulativeMarginalProbability(1 , 300.0 )
                               - 1.0 ) , 1E-15 );
        BOOST_CHECK_SMALL( std::fabs(
            distribution.getCumulativeMarginalProbability(1 , -200.0 )
                               - 0.0 ) , 1E-15 );
    }
}

BOOST_AUTO_TEST_CASE( test_standard_deviation_scaling )
{

    using namespace tudat::statistics;

    // Generate random datapoints
    Eigen::VectorXd mean(2);
    mean << 0.0 , 0.0 ;
    Eigen::VectorXd standardDeviation(2);
    standardDeviation << 2.0 , 0.5 ;
    Eigen::VectorXd lowerBound = mean - standardDeviation * std::sqrt(3.0) ;
    Eigen::VectorXd upperBound = mean + standardDeviation * std::sqrt(3.0) ;

    int numberOfSamples = 1E6;
    int seed = 100 ;
    std::vector<Eigen::VectorXd> samples = generateRandomVectorUniform(
                seed , numberOfSamples , lowerBound , upperBound );

    KernelDensityDistribution distribution2( samples );

    Eigen::VectorXd sampleMean = distribution2.getSampleMean() ;
    Eigen::VectorXd sampleStandardDeviation = distribution2.getSampleStandardDeviation() ;

    BOOST_CHECK_SMALL( std::fabs( sampleMean(0) - 0.0 ) , 2E-3 );
    BOOST_CHECK_SMALL( std::fabs( sampleMean(1) - 0.0 ) , 3E-4 );
    BOOST_CHECK_SMALL( std::fabs( sampleStandardDeviation(0) - 2.0 ) , 1E-3 );
    BOOST_CHECK_SMALL( std::fabs( sampleStandardDeviation(1) - 0.5 ) , 1E-4 );

    // Check get probability density approx uniform
    double probabilityDensity = 1.0 / ( ( upperBound(0) - lowerBound(0) )*( upperBound(1) - lowerBound(1) ) );

    Eigen::VectorXd x(2);
    x << 0.0 , 0.0 ;

    BOOST_CHECK_CLOSE_FRACTION( distribution2.getProbabilityDensity( x ) , probabilityDensity , 2E-2 );
    BOOST_CHECK_SMALL( std::fabs( distribution2.getProbabilityDensity( x ) - probabilityDensity ) , 2E-3);

    x << -1.0 , 0.3 ;
    BOOST_CHECK_CLOSE_FRACTION( distribution2.getProbabilityDensity( x ) , probabilityDensity , 2E-4 );

    Eigen::VectorXd standardDeviationAdjusted(2);
    standardDeviationAdjusted << 0.5 , 3.0 ;
    KernelDensityDistribution distribution( samples , 1.0 , KernelType::Gaussian , true , standardDeviationAdjusted );

    Eigen::VectorXd newSampleMean = distribution.getSampleMean() ;
    Eigen::VectorXd newSampleStandardDeviation = distribution.getSampleStandardDeviation() ;

    BOOST_CHECK_SMALL( std::fabs( newSampleMean(0) - 0.0 ) , 5E-3 );
    BOOST_CHECK_SMALL( std::fabs( newSampleMean(1) - 0.0 ) , 2E-3 );
    BOOST_CHECK_SMALL( std::fabs( newSampleStandardDeviation(0) - 0.5 ) , 1E-13 );
    BOOST_CHECK_SMALL( std::fabs( newSampleStandardDeviation(1) - 3.0 ) , 1E-13 );

    // Check get probability density approx uniform
    lowerBound = mean - standardDeviationAdjusted * std::sqrt(3.0) ;
    upperBound = mean + standardDeviationAdjusted * std::sqrt(3.0) ;
    probabilityDensity = 1.0 / ( ( upperBound(0) - lowerBound(0) )*( upperBound(1) - lowerBound(1) ) );
    x << 0.0 , 0.0 ;
    BOOST_CHECK_CLOSE_FRACTION( distribution.getProbabilityDensity( x ) , probabilityDensity , 2E-2 );
    BOOST_CHECK_SMALL( std::fabs( distribution.getProbabilityDensity( x ) - probabilityDensity ) , 8E-4);

    x << -0.2 , 2.3 ;
    BOOST_CHECK_CLOSE_FRACTION( distribution.getProbabilityDensity( x ) , probabilityDensity , 5E-4 );
}

BOOST_AUTO_TEST_CASE( test_Probability_Density_etc_3D )
{
    using namespace tudat::statistics;

    Eigen::VectorXd sample(3);
    std::vector< Eigen::VectorXd > samples(0);

    sample << 3.0 , 2E-2 , 0.1 ;
    samples.push_back( sample );
    sample << -1.0 , 5E-2 , 1.0;
    samples.push_back( sample );
    sample << 5.0 , 2E-1 , 3.0;
    samples.push_back( sample );
    sample << 2.0 , 1E-3 , 1.5;
    samples.push_back( sample );

    KernelDensityDistribution distribution(samples , 1.0 , KernelType::Gaussian ) ;

    // Test probability density
    {
        Eigen::VectorXd location(3);
        location <<  1.5 , 1E-2 , 0.8;

        double computedDensity = distribution.getProbabilityDensity( location );
        double expectedDensity = 4.653809309094347e-01 ;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );
    }

    // Test marginal cumulative probability
    {
        Eigen::VectorXd location(3);
        location <<  1.5 , 1E-2 , 0.8;

        int marginal = 0 ;
        double computedDensity = distribution.getCumulativeMarginalProbability( marginal , location(marginal) );
        double expectedDensity = 3.829580307170373e-01 ;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );

        marginal = 1 ;
        computedDensity = distribution.getCumulativeMarginalProbability( marginal , location(marginal) );
        expectedDensity = 2.674493565829487e-01 ;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );
    }

    // Test conditional marginal cumulative probability
    {
        Eigen::VectorXd location(3);
        location <<  1.5 , 1E-2 , 0.8;

        int marginal = 0 ;
        std::vector<int> conditionDimensions(0);
        conditionDimensions.push_back( 1 );
        std::vector<double> conditions(0);
        conditions.push_back( 0.1 );

        double computedDensity = distribution.getCumulativeConditionalMarginalProbability(
                    conditionDimensions, conditions , marginal , location(marginal) );
        double expectedDensity = 4.970339167925200e-01 ;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );

        conditionDimensions[0] = 0;
        conditions[0] = 0.4;
        marginal = 2 ;
        computedDensity = distribution.getCumulativeConditionalMarginalProbability(
                    conditionDimensions, conditions , marginal , location(marginal) );
        expectedDensity = 3.932443313127392e-01 ;

        // Check correct density
        BOOST_CHECK_CLOSE_FRACTION( computedDensity , expectedDensity , 1E-15 );
    }
}

BOOST_AUTO_TEST_CASE( test_Probability_Functions_uncorrelated_2D )
{
    using namespace tudat::statistics;

    Eigen::VectorXd sample(2);
    std::vector< Eigen::VectorXd > samples(0);

    sample << 3.0 , 0.1  ;
    samples.push_back( sample );
    sample << 5.0 , 0.1 ;
    samples.push_back( sample );
    sample << 7.0 , 0.1 ;
    samples.push_back( sample );

    KernelDensityDistribution distribution(samples , 1.0 , KernelType::Gaussian ) ;
    Eigen::VectorXd bandWidth(2);
    bandWidth << 1.0 , 0.3 ;
    distribution.setBandWidth( bandWidth );

    GaussianDistributiond gaussianDistribution1( samples[0](0) , 1.0 );
    GaussianDistributiond gaussianDistribution2( samples[1](0) , 1.0 );
    GaussianDistributiond gaussianDistribution3( samples[2](0) , 1.0 );
    GaussianDistributiond gaussianDistribution4( samples[0](1) , 0.3 );

    // Test probability density
    {
        Eigen::VectorXd location(2);
        location <<  1.5 , 1E-2 ;

        double computedDensity = distribution.getProbabilityDensity( location );
        double expectedDensity = (gaussianDistribution1.getProbabilityDensity( location(0) ) +
                                  gaussianDistribution2.getProbabilityDensity( location(0) ) +
                                  gaussianDistribution3.getProbabilityDensity( location(0) )) *
                gaussianDistribution4.getProbabilityDensity( location(1) ) /3.0;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );
    }

    // Test cumulative probability
    {
        Eigen::VectorXd location(2);
        location <<  2.5 , 1E-1 ;

        double computedDensity = distribution.getCumulativeProbability( location );
        double expectedDensity = (gaussianDistribution1.getCumulativeProbability( location(0) ) +
                                  gaussianDistribution2.getCumulativeProbability( location(0) ) +
                                  gaussianDistribution3.getCumulativeProbability( location(0) )) *
                gaussianDistribution4.getCumulativeProbability( location(1) ) /3.0;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );
    }

    // Test marginal cumulative probability
    {
        Eigen::VectorXd location(2);
        location <<  2.5 , 1E-1 ;
        int marginal = 0 ;

        double computedDensity = distribution.getCumulativeMarginalProbability( marginal , location(marginal) );
        double expectedDensity = (gaussianDistribution1.getCumulativeProbability( location(0) ) +
                                  gaussianDistribution2.getCumulativeProbability( location(0) ) +
                                  gaussianDistribution3.getCumulativeProbability( location(0) )) /3.0;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );

        // Test 2 other marginal
        marginal = 1 ;
        computedDensity = distribution.getCumulativeMarginalProbability( marginal , location(marginal) );
        expectedDensity = gaussianDistribution4.getCumulativeProbability( location(1) );

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );
    }

    // Test conditional marginal cumulative probability
    {
        double location;
        std::vector<int> conditionDimensions(0);
        std::vector<double> conditions(0);

        location = 2.5 ;
        int marginal = 0 ;
        conditionDimensions.push_back( 1 );
        conditions.push_back( 0.2 );

        double computedDensity = distribution.getCumulativeConditionalMarginalProbability(
                    conditionDimensions, conditions , marginal , location );
        double expectedDensity = (gaussianDistribution1.getCumulativeProbability( location ) +
                                  gaussianDistribution2.getCumulativeProbability( location ) +
                                  gaussianDistribution3.getCumulativeProbability( location )) /3.0;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );

        // Test 2 other marginal
        marginal = 1 ;
        conditionDimensions[0] = 0;
        conditions[0] = 2.5;
        location = 0.2 ;
        computedDensity = distribution.getCumulativeConditionalMarginalProbability(
                    conditionDimensions, conditions , marginal , location );
        expectedDensity = gaussianDistribution4.getCumulativeProbability( location );

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ) , 1E-15 );
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
