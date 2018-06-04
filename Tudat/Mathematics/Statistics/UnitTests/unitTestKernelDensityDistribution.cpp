/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <vector>
#include <limits>

#include <Eigen/Core>

#include <boost/test/unit_test.hpp>
#include <boost/random.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Mathematics/Statistics/kernelDensityDistribution.h"


namespace tudat
{
namespace unit_tests
{

//! Creates a std::vector of linear spaced values
std::vector< double > linspace( double start, double end, int N )
{
    std::vector < double > x( 0 );

    if ( N > 1 )
    {
        for( int i = 0; i < N; i++ )
        {
            x.push_back( start + ( end - start )*i /( N - 1.0 ) );
        }
    }
    return x;
}

//! Generator random vector using pseudo random generator
std::vector<Eigen::VectorXd> generateRandomVectorUniform(
        int seed, int numberOfSamples,
        Eigen::VectorXd lowerBound, Eigen::VectorXd upperBound)
{
    // Compute properties
    Eigen::VectorXd width = upperBound - lowerBound;
    Eigen::VectorXd average = (upperBound + lowerBound ) / 2.0;

    // Setup Random generator
    typedef boost::mt19937 RandomGeneratorType; // Mersenne Twister
    RandomGeneratorType randomGenerator(seed);              // Create random generator

    boost::uniform_real< > uniformDistribution( 0.0, 1.0 ); //
    boost::variate_generator< RandomGeneratorType, boost::uniform_real< > >
            Dice(randomGenerator, uniformDistribution); // define random generator

    std::vector< Eigen::VectorXd > randomSamples(numberOfSamples);
    Eigen::VectorXd randomSample( lowerBound.rows( ) );

    // Sample
    for(int i = 0; i < numberOfSamples; i++ ){ // Generate N samples
        for(int j = 0; j < randomSample.rows( ); j++){ // Generate vector of samples
            randomSample(j) = Dice( ) - 0.5;
        }
        randomSamples[i] = randomSample.cwiseProduct(width) + average;
    }
    return randomSamples;
}

BOOST_AUTO_TEST_SUITE( test_Kernel_Density_Distribution )

using tudat::mathematical_constants::PI;

//! Test whether optimal bandwidths are correctly computed (compared to Matlab results)
BOOST_AUTO_TEST_CASE( testOptimalKernelBandwidth )
{

    using namespace tudat::statistics;

    std::vector< Eigen::VectorXd > data( 0 );
    Eigen::VectorXd sample( 2 );
    sample << 2.2538, 1.177;
    data.push_back( sample );
    sample << 0.76529, 0.356;
    data.push_back( sample );
    sample << 1.5179, 1.14;
    data.push_back( sample );
    sample << 2.0972, 0.34093;
    data.push_back( sample );
    sample << 2.6727, 1.301;
    data.push_back( sample );
    sample << 2.8779, 0.48998;
    data.push_back( sample );
    sample << 1.6416, 0.27523;
    data.push_back( sample );
    sample << 0.41587, 0.35152;
    data.push_back( sample );
    sample << 0.44788, 0.86246;
    data.push_back( sample );
    sample << 0.77252, 0.6626;
    data.push_back( sample );

    KernelDensityDistribution distribution( data );

    // Generated expected Bandwidth using MATLAB
    Eigen::VectorXd expectedBandwidth( 2 );
    expectedBandwidth << 0.819018940571854, 0.263389630191125;

    Eigen::VectorXd computedOptimalBandwidth( 2 );
    computedOptimalBandwidth = distribution.getOptimalBandWidth( );

    // Compare Matlab against Tudal results; results are similar but not identical due to slightly different algorithms/
    BOOST_CHECK_CLOSE_FRACTION( expectedBandwidth( 0 ), computedOptimalBandwidth( 0 ), 2E-5 );
    BOOST_CHECK_CLOSE_FRACTION( expectedBandwidth( 1 ), computedOptimalBandwidth( 1 ), 2E-5 );

}

//! Test if the Epanechnikov distribution is correctly implemented.
BOOST_AUTO_TEST_CASE( testEpanechnikovKernel )
{
    tudat::statistics::EpanechnikovKernelDistribution distribution( 1.0, 2.0 );

    // Test theoretical values of Cdf
    BOOST_CHECK_CLOSE_FRACTION( distribution.evaluateCdf( 5.0 ), 1.0,
                                4.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( distribution.evaluateCdf( 3.0 ), 1.0,
                                4.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( distribution.evaluateCdf( 1.0 ), 0.5,
                                4.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( std::fabs( distribution.evaluateCdf( -1.0 ) - 0.0 ),
                       4.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( std::fabs( distribution.evaluateCdf( -2.0 ) - 0.0 ),
                       4.0 * std::numeric_limits< double >::epsilon( ) );

    // Test theoretical values of Pdf (zero outside of given range)
    BOOST_CHECK_SMALL( std::fabs( distribution.evaluatePdf( -1.0 ) - 0.0 ),
                       4.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( std::fabs( distribution.evaluatePdf( 3.0 ) - 0.0 ),
                       4.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( std::fabs( distribution.evaluatePdf( 4.0 ) - 0.0 ),
                       4.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( std::fabs( distribution.evaluatePdf( -3.0 ) - 0.0 ),
                       4.0 * std::numeric_limits< double >::epsilon( ) );

    // Compare pdf and cdf at center points, compared to Matlab implementation
    BOOST_CHECK_SMALL( std::fabs( distribution.evaluatePdf( 0.0 ) - 0.28125 ),
                       4.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( std::fabs( distribution.evaluateCdf( 0.0 ) - 0.15625 ),
                       4.0 * std::numeric_limits< double >::epsilon( ) );
}

//! Test 1-Dimensional kernel density distribution, compared against Matlab implementation.
BOOST_AUTO_TEST_CASE( testKernelProbabilityDensity1D )
{
    using namespace tudat::statistics;

    Eigen::VectorXd sample( 1 );
    std::vector< Eigen::VectorXd > samples( 0 );

    sample( 0 ) = 3.0;
    samples.push_back( sample );
    sample( 0 ) = 1.0;
    samples.push_back( sample );
    sample( 0 ) = 4.0;
    samples.push_back( sample );
    sample( 0 ) = 7.0;
    samples.push_back( sample );

    // Create kernel
    KernelDensityDistribution distribution(samples, 1.0, KernelType::gaussian_kernel );

    Eigen::VectorXd location( 1 );
    // Case 1
    {
        // Compure pdf
        location << 5.396093156208099e-01;
        double computedDensity = distribution.evaluatePdf( location );

        // Data obtained from MATLAB ksdensity function
        double expectedDensity = 8.426926493705324e-02;

        // Check that the same bandwidths are used compared with MATLAB
        Eigen::VectorXd expectedBandwidth( 1 );
        expectedBandwidth << 1.785192502061299;

        Eigen::VectorXd computedBandwidth = distribution.getBandWidth( );
        Eigen::VectorXd difference = expectedBandwidth - computedBandwidth;

        BOOST_CHECK_SMALL( std::fabs( difference( 0 ) ), 4.0 * std::numeric_limits< double >::epsilon( ) );

        // Check pdf against Matlab computation
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Case 2
    {
        Eigen::VectorXd location( 1 );
        location << 3.915600227210264;

        double computedDensity = distribution.evaluatePdf( location );
        double expectedDensity = 1.320686749820565e-01;

        // Check pdf against Matlab computation
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Case 3
    {
        Eigen::VectorXd location( 1 );
        location << 8.135588866697081;

        double computedDensity = distribution.evaluatePdf( location );
        double expectedDensity = 5.036314749096298e-02;

        // Check pdf against Matlab computation
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );
    }
}

//! Test 2-Dimensional kernel density distribution, compared against Matlab implementation.
BOOST_AUTO_TEST_CASE( testKernelProbabilityDensity2D )
{
    using namespace tudat::statistics;

    Eigen::VectorXd sample( 2 );
    std::vector< Eigen::VectorXd > samples( 0 );

    sample << 3.0, 2E-2;
    samples.push_back( sample );
    sample << -1.0, 5E-2;
    samples.push_back( sample );
    sample << 5.0, 2E-1;
    samples.push_back( sample );
    sample << 2.0, 1E-3;
    samples.push_back( sample );

    // Create kernel
    KernelDensityDistribution distribution( samples, 1.0, KernelType::gaussian_kernel );

    // Case 1
    {
        Eigen::VectorXd location( 2 );
        location << 1.141869732328655, 2.995235933344523e-02;

        double computedDensity = distribution.evaluatePdf( location );
        double expectedDensity = 1.136916246361795;

        // Check that the same bandwidths are used compared with MATLAB
        Eigen::VectorXd expectedBandwidth( 2 );
        expectedBandwidth << 1.765086418052112, 2.882974482818450e-02;
        Eigen::VectorXd computedBandwidth = distribution.getBandWidth( );
        Eigen::VectorXd difference = expectedBandwidth - computedBandwidth;
        BOOST_CHECK_SMALL( std::fabs( difference( 0 ) ), 4.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs( difference( 1 ) ), 4.0 * std::numeric_limits< double >::epsilon( ) );

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Case 2
    {
        Eigen::VectorXd location( 2 );
        location << -2.862738183470956, 1.582207969089993e-01;

        double computedDensity = distribution.evaluatePdf( location );
        double expectedDensity = 4.038687551312721e-04;

        // Check pdf against Matlab computation
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Case 3
    {
        Eigen::VectorXd location( 2 );
        location << 6.862738183470957, 2.995235933344523e-02;

        double computedDensity = distribution.evaluatePdf( location );
        double expectedDensity = 7.784149895640656e-02;

        // Check pdf against Matlab computation
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Compare marginal CDF against Matlab
    {
        // Marginal CDF compare with MATLAB
        Eigen::VectorXd bandwidth( 2 );
        bandwidth << 1.785192502061299, 2.915814420033455e-02;
        distribution.setBandWidth( bandwidth );

         // Compare against Matlab results,
        BOOST_CHECK_CLOSE_FRACTION(
                    distribution.evaluateCumulativeMarginalProbability( 0, 2.276047714155371E-1 ),
                    2.446324562687310E-1, 5E-3 );
        BOOST_CHECK_CLOSE_FRACTION(
                    distribution.evaluateCumulativeMarginalProbability( 1, 6.083875672099921e-02 ),
                    6.360523509878839e-01, 5E-3 );
        BOOST_CHECK_CLOSE_FRACTION(
                    distribution.evaluateCumulativeMarginalProbability( 1, 9.861136936766661e-02 ),
                    7.371490745984073e-01, 5E-3 );

        // Check theoretical CDF at edge of domain.
        BOOST_CHECK_SMALL( std::fabs(
            distribution.evaluateCumulativeMarginalProbability( 1, 300.0 )
                               - 1.0 ), 4.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs(
            distribution.evaluateCumulativeMarginalProbability( 1, -200.0 )
                               - 0.0 ), 4.0 * std::numeric_limits< double >::epsilon( ) );
    }
}

//! Test whether scaling of sample standard deviation is done correctly
BOOST_AUTO_TEST_CASE( testStandardDeviationScaling )
{

    using namespace tudat::statistics;

    // Generate random datapoints
    Eigen::VectorXd mean( 2 );
    mean << 0.0, 0.0;
    Eigen::VectorXd standardDeviation( 2 );
    standardDeviation << 2.0, 0.5;
    Eigen::VectorXd lowerBound = mean - standardDeviation * std::sqrt( 3.0 );
    Eigen::VectorXd upperBound = mean + standardDeviation * std::sqrt( 3.0 );

    int numberOfSamples = 1E6;
    int seed = 100;
    std::vector< Eigen::VectorXd > samples = generateRandomVectorUniform(
                seed, numberOfSamples, lowerBound, upperBound );

    // Create distribution
    KernelDensityDistribution distribution2( samples );

    Eigen::VectorXd sampleMean = distribution2.getSampleMean( );
    Eigen::VectorXd sampleStandardDeviation = distribution2.getSampleStandardDeviation( );

    // Check sample mean and standard deviation
    BOOST_CHECK_SMALL( std::fabs( sampleMean( 0 ) - 0.0 ), 5E-3 );
    BOOST_CHECK_SMALL( std::fabs( sampleMean( 1 ) - 0.0 ), 5E-3 );
    BOOST_CHECK_SMALL( std::fabs( sampleStandardDeviation( 0 ) - 2.0 ), 1E-3 );
    BOOST_CHECK_SMALL( std::fabs( sampleStandardDeviation( 1 ) - 0.5 ), 1E-3 );

    // Check whether probability density function is approximately uniform.
    double probabilityDensity = 1.0 / ( ( upperBound( 0 ) - lowerBound( 0 ) ) * ( upperBound( 1 ) - lowerBound( 1 ) ) );

    Eigen::VectorXd x( 2 );
    x << 0.0, 0.0;
    BOOST_CHECK_SMALL( std::fabs( distribution2.evaluatePdf( x ) - probabilityDensity ), 2E-3 );

    x << -1.0, 0.3;
    BOOST_CHECK_CLOSE_FRACTION( distribution2.evaluatePdf( x ), probabilityDensity, 2E-4 );

    // Set standard deviation scaling
    Eigen::VectorXd standardDeviationAdjusted( 2 );
    standardDeviationAdjusted << 0.5, 3.0;

    // Create scaled kernel distribution
    KernelDensityDistribution distribution( samples, 1.0, KernelType::gaussian_kernel, standardDeviationAdjusted );

    Eigen::VectorXd newSampleMean = distribution.getSampleMean( );
    Eigen::VectorXd newSampleStandardDeviation = distribution.getSampleStandardDeviation( );

    // Recheck mean.
    BOOST_CHECK_SMALL( std::fabs( newSampleMean( 0 ) - 0.0 ), 5E-3 );
    BOOST_CHECK_SMALL( std::fabs( newSampleMean( 1 ) - 0.0 ), 5E-3 );

    // Check whether sample standard deviation conforms to required values.
    BOOST_CHECK_SMALL( std::fabs( newSampleStandardDeviation( 0 ) - 0.5 ), 1E-13 );
    BOOST_CHECK_SMALL( std::fabs( newSampleStandardDeviation( 1 ) - 3.0 ), 1E-13 );

    // Check whether probability density function is approximately uniform.
    lowerBound = mean - standardDeviationAdjusted * std::sqrt( 3.0 );
    upperBound = mean + standardDeviationAdjusted * std::sqrt( 3.0 );
    probabilityDensity = 1.0 / ( ( upperBound( 0 ) - lowerBound( 0 ) ) * ( upperBound( 1 ) - lowerBound( 1 ) ) );

    x << 0.0, 0.0;
    BOOST_CHECK_SMALL( std::fabs( distribution.evaluatePdf( x ) - probabilityDensity ), 1.0E-3 );

    x << -0.2, 2.3;
    BOOST_CHECK_CLOSE_FRACTION( distribution.evaluatePdf( x ), probabilityDensity, 1.0E-3 );
}

//! Test 3-Dimensional kernel density distribution, compared against Matlab implementation.
BOOST_AUTO_TEST_CASE( testKernelProbabilityDensity3D )
{
    using namespace tudat::statistics;

    Eigen::VectorXd sample( 3 );
    std::vector< Eigen::VectorXd > samples( 0 );

    sample << 3.0, 2E-2, 0.1;
    samples.push_back( sample );
    sample << -1.0, 5E-2, 1.0;
    samples.push_back( sample );
    sample << 5.0, 2E-1, 3.0;
    samples.push_back( sample );
    sample << 2.0, 1E-3, 1.5;
    samples.push_back( sample );

    KernelDensityDistribution distribution(samples, 1.0, KernelType::gaussian_kernel );

    // Test probability density
    {
        Eigen::VectorXd location( 3 );
        location <<  1.5, 1E-2, 0.8;

        double computedDensity = distribution.evaluatePdf( location );
        double expectedDensity = 4.653809309094347e-01;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Test marginal cumulative probability
    {
        Eigen::VectorXd location( 3 );
        location <<  1.5, 1E-2, 0.8;

        int marginal = 0;
        double computedDensity = distribution.evaluateCumulativeMarginalProbability( marginal, location( marginal ) );
        double expectedDensity = 3.829580307170373e-01;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );

        marginal = 1;
        computedDensity = distribution.evaluateCumulativeMarginalProbability( marginal, location( marginal ) );
        expectedDensity = 2.674493565829487e-01;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Test conditional marginal cumulative probability
    {
        Eigen::VectorXd location( 3 );
        location <<  1.5, 1E-2, 0.8;

        int marginal = 0;
        std::vector< int > conditionDimensions( 0 );
        conditionDimensions.push_back( 1 );
        std::vector< double > conditions( 0 );
        conditions.push_back( 0.1 );

        double computedDensity = distribution.evaluateCumulativeConditionalMarginalProbability(
                    conditionDimensions, conditions, marginal, location( marginal ) );
        double expectedDensity = 4.970339167925200e-01;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );

        conditionDimensions[0] = 0;
        conditions[0] = 0.4;
        marginal = 2;
        computedDensity = distribution.evaluateCumulativeConditionalMarginalProbability(
                    conditionDimensions, conditions, marginal, location( marginal ) );
        expectedDensity = 3.932443313127392e-01;

        // Check correct density
        BOOST_CHECK_CLOSE_FRACTION( computedDensity, expectedDensity, 4.0 * std::numeric_limits< double >::epsilon( ) );
    }
}

//! Test kernel density function for uncorrelated input data.
BOOST_AUTO_TEST_CASE( testProbabilityFunctionsUncorrelated2D )
{
    using namespace tudat::statistics;

    Eigen::VectorXd sample( 2 );
    std::vector< Eigen::VectorXd > samples( 0 );

    // Define uncorrelated input data.
    sample << 3.0, 0.1 ;
    samples.push_back( sample );
    sample << 5.0, 0.1;
    samples.push_back( sample );
    sample << 7.0, 0.1;
    samples.push_back( sample );

    // Manually define bandwidth
    Eigen::VectorXd bandWidth( 2 );
    bandWidth << 1.0, 0.3;

    KernelDensityDistribution distribution( samples, 1.0, KernelType::gaussian_kernel, Eigen::VectorXd::Zero( 0 ), bandWidth );

    // Create manual Gaussian distributions
    std::shared_ptr< ContinuousProbabilityDistribution< double > > gaussianDistribution1 =
            createBoostRandomVariable( normal_boost_distribution, { samples[ 0 ]( 0 ), 1.0  } );

    std::shared_ptr< ContinuousProbabilityDistribution< double > > gaussianDistribution2 =
            createBoostRandomVariable( normal_boost_distribution, { samples[ 1 ]( 0 ), 1.0  } );

    std::shared_ptr< ContinuousProbabilityDistribution< double > > gaussianDistribution3 =
            createBoostRandomVariable( normal_boost_distribution, { samples[ 2 ]( 0 ), 1.0  } );

    std::shared_ptr< ContinuousProbabilityDistribution< double > > gaussianDistribution4 =
            createBoostRandomVariable( normal_boost_distribution, { samples[ 0 ]( 1 ), 0.3  } );

    // Test probability density
    {
        Eigen::VectorXd location( 2 );
        location <<  1.5, 1E-2;

        // Compute pdf from kernel and theoretical value.
        double computedDensity = distribution.evaluatePdf( location );
        double expectedDensity = ( gaussianDistribution1->evaluatePdf( location( 0 ) ) +
                                   gaussianDistribution2->evaluatePdf( location( 0 ) ) +
                                   gaussianDistribution3->evaluatePdf( location( 0 ) ) ) *
                gaussianDistribution4->evaluatePdf( location( 1 ) ) / 3.0;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Test cumulative probability
    {
        Eigen::VectorXd location( 2 );
        location <<  2.5, 1E-1;

        // Compute cdf from kernel and theoretical value.
        double computedDensity = distribution.evaluateCdf( location );
        double expectedDensity = (gaussianDistribution1->evaluateCdf( location( 0 ) ) +
                                  gaussianDistribution2->evaluateCdf( location( 0 ) ) +
                                  gaussianDistribution3->evaluateCdf( location( 0 ) ) ) *
                gaussianDistribution4->evaluateCdf( location( 1 ) ) / 3.0;

        // Check correct density
        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Test marginal cumulative probability
    {
        Eigen::VectorXd location( 2 );
        location <<  2.5, 1E-1;
        int marginal = 0;

        // Compute marginal cdf from kernel and theoretical value (marginalDimension = 0).
        double computedDensity = distribution.evaluateCumulativeMarginalProbability( marginal, location( marginal ) );
        double expectedDensity = ( gaussianDistribution1->evaluateCdf( location( 0 ) ) +
                                   gaussianDistribution2->evaluateCdf( location( 0 ) ) +
                                   gaussianDistribution3->evaluateCdf( location( 0 ) ) ) / 3.0;

        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );

        // Compute marginal cdf from kernel and theoretical value (marginalDimension = 1).
        marginal = 1;
        computedDensity = distribution.evaluateCumulativeMarginalProbability( marginal, location( marginal ) );
        expectedDensity = gaussianDistribution4->evaluateCdf( location( 1 ) );

        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Test conditional marginal cumulative probability
    {
        double location;
        std::vector< int > conditionDimensions( 0 );
        std::vector< double > conditions( 0 );

        location = 2.5;
        int marginal = 0;
        conditionDimensions.push_back( 1 );
        conditions.push_back( 0.2 );

        // Compute marginal cdf from kernel and theoretical value (marginalDimension = 0).
        double computedDensity = distribution.evaluateCumulativeConditionalMarginalProbability(
                    conditionDimensions, conditions, marginal, location );
        double expectedDensity = ( gaussianDistribution1->evaluateCdf( location ) +
                                   gaussianDistribution2->evaluateCdf( location ) +
                                   gaussianDistribution3->evaluateCdf( location ) ) / 3.0;

        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );

        // Compute conditional marginal cdf from kernel and theoretical value (marginalDimension = 1).
        marginal = 1;
        conditionDimensions[0] = 0;
        conditions[ 0 ]  = 2.5;
        location = 0.2;
        computedDensity = distribution.evaluateCumulativeConditionalMarginalProbability(
                    conditionDimensions, conditions, marginal, location );
        expectedDensity = gaussianDistribution4->evaluateCdf( location );

        BOOST_CHECK_SMALL( std::fabs( computedDensity - expectedDensity ), 4.0 * std::numeric_limits< double >::epsilon( ) );
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

