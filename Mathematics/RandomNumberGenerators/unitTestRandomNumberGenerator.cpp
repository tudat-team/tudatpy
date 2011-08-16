/*! \file unitTestRandomNumberGenerator.cpp
 *    Source file that defines a unit test that tests the random number
 *    generator included in Tudat.
 *
 *    Path              : /Mathematics/RandomNumberGenerators/
 *    Version           : 6
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : E.A.G. Heeren
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : E.A.G.Heeren@student.tudelft.nl
 *
 *    Date created      : 7 January, 2011
 *    Last modified     : 29 July, 2011
 *
 *    References
 *      Spiegel, M.R., Stephens, L.J. Statistics, Fourth Edition, Schaum's
 *          Outlines, McGraw-Hill, 2008.
 *      Walpole, R.E., et al. Probabilty and Statistics for Engineers &
 *          Scientists, Seventh Edition, Prentice Hall, NJ, 2002.
 *      Texas A&M University. Chi Square Calculator,
 *          http://www.stat.tamu.edu/~west/applets/chisqdemo.html, last
 *          accessed: 7 July, 2011.
 *      Faculty of Mathematics and Informatics, University of Sofia.
 *          http://www.fmi.uni-sofia.bg/vesta/virtual_labs/interval
 *          /interval6.html,last accessed: 7 July, 2011.
 *      Wikipedia. Variance: Population variance and sample variance,
 *          http://en.wikipedia.org/wiki/Variance#Population_variance_and
 *          _sample_variance, last accessed: 8 July, 2011.
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
 *
 *      The code tests for the quality of the random distributions by testing
 *      the sample mean and variance against the true mean and variance using
 *      standard error estimates for a normal distribution. This is permissible
 *      under Central Limit Theorem, so long as the number of samples is large.
 *      These estimates are stated with a confidence of 99.96%.
 *      If the unit test fails due to Test 5 for the uniform distribution test
 *      or Test 2 for the exponential distribution, it is recommended that the
 *      unit test is run again to rule out statistical anomalies.
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110107    K. Kumar          First creation of code.
 *      110207    K. Kumar          Updated code to conform to protocol; added
 *                                  cerr statements.
 *      110209    K. Kumar          Updated tests and added test for mean and
 *                                  variance of uniform distribution; added
 *                                  note and reference.
 *      110707    K. Kumar          Corrected variance calculation errors;
 *                                  updated references.
 *      110707    K. Kumar          Updated code to be compatible with file/
 *                                  class name change for uniform random
 *                                  number generator. Added unit test for
 *                                  exponential random number generator.
 *      110728    K. Kumar          Added unit test for normal random number
 *                                  generator.
 *      110729    E.A.G. Heeren     Minor changes to comments in normal random
 *                                  number generator.
 */

// Include statements.
#include "unitTestRandomNumberGenerator.h"

#include <fstream>

// Using declarations.
using mathematics::computeAbsoluteValue;
using mathematics::raiseToIntegerPower;
using mathematics::MACHINE_PRECISION_DOUBLES;
using mathematics::computeSampleMean;
using mathematics::computeSampleVariance;
using std::cerr;
using std::endl;
using std::map;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test implementation of uniform random number generator class.
bool testUniformRandomNumberGenerator( )
{
    // Summary of tests.
    // Test 1: Get a random integer ( 64-bit ).
    // Test 2: Get a random integer ( 32-bit ).
    // Test 3: Get a normalized random double in the interval [ 0, 1 ].
    // Test 4: Get a random plus/minus sign.
    // Test 5: Compute sample mean and standard distrubution of distribution
    //         and compare to analytical results.

    // Test result initialised to false.
    bool isUniformRandomNumberGeneratorErroneous = false;

    // Create uniform random number generator object.
    UniformRandomNumberGenerator uniformRandomNumbers( time( NULL ) );

    // Results computed using implementation of uniform random number generator
    // class.
    // Test 1: Get a random integer ( 64-bit ).
    int computedResultForTest1 = uniformRandomNumbers
            .getUniformlyDistributedRandomInteger( );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( fmod( computedResultForTest1,
               static_cast< double >( computedResultForTest1 ) )
         > MACHINE_PRECISION_DOUBLES )
    {
        isUniformRandomNumberGeneratorErroneous = true;

        cerr << "The computed value ( " << computedResultForTest1
                << " ) using the uniform 64-bit random integer generator "
                << "is not an integer" << endl;
    }

    // Test 2: Get a random integer ( 32-bit ).
    int computedResultForTest2 = uniformRandomNumbers
                                 .getUniformlyDistributedRandom32BitInteger( );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( fmod( computedResultForTest2,
               static_cast< double >( computedResultForTest2 ) )
         > MACHINE_PRECISION_DOUBLES )
    {
        isUniformRandomNumberGeneratorErroneous = true;

        cerr << "The computed value ( " << computedResultForTest2
                << " ) using the uniform 32-bit random integer generator "
                << "is not an integer" << endl;
    }

    // Test 3: Get a normalized random double in the interval [ 0, 1 ].
    double computedResultForTest3
            = uniformRandomNumbers
            .getUniformlyDistributedNormalizedRandomDouble( );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( computedResultForTest3 < 0.0 && computedResultForTest3 > 1.0 )
    {
        isUniformRandomNumberGeneratorErroneous = true;

        cerr << "The computed value ( " << computedResultForTest3
                << " ) using the uniform normalized random number generator "
                << "is not within the prescribed bounds [ 0, 1 ]" << endl;
    }

    // Test 4: Get a random plus/minus sign.
    double computedResultForTest4 = uniformRandomNumbers
            .getUniformlyDistributedRandomPlusMinusSign( );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( computeAbsoluteValue( computedResultForTest4 ) - 1.0
         > MACHINE_PRECISION_DOUBLES )
    {
        isUniformRandomNumberGeneratorErroneous = true;

        cerr << "The computed value ( " << computedResultForTest4
                << " ) using the uniform plus/minus sign random generator "
                << "is not equal to either +1 or -1" << endl;
    }

    // Test 5: Compute sample mean and variance of distribution and compare to
    //         analytical results.
    // Declare and initialize number of samples.
    unsigned int numberOfSamples = 100000;

    // Declare vector of sample uniform random values.
    vector< double > sampleOfUniformRandomValues;

    // Fill map of sample uniform random values.
    for ( unsigned int i = 0; i < numberOfSamples; i++ )
    {
        sampleOfUniformRandomValues.push_back(
                    uniformRandomNumbers
                    .getUniformlyDistributedNormalizedRandomDouble( ) );
    }

    // Estimate sample mean of the distribution.
    double sampleMean = computeSampleMean( sampleOfUniformRandomValues );

    // Estimate sample variance.
    double sampleVariance = computeSampleVariance(
                sampleOfUniformRandomValues );

    // Compute differences between computed and expected results.
    // Test sample mean versus true mean with a confidence of 99.96%.
    // This means that the test below guarentees that this unit test will
    // fail on average every 2500 times it is run.
    // Also test sample variance versus true variance with a confidence of
    // 99.96%. The true mean and true variance of a uniform distribution
    // from 0 to 1 are 1/2 and 1/12 respectively.
    // The tolerances were computed using confidence intervals for the sample
    // mean and variance of a normal distribution and table lookup.
    if ( computeAbsoluteValue( sampleMean - 0.5 )
         > 3.49 * 1.0 / sqrt( 12.0 )
           * 1.0 / sqrt( static_cast< double >( numberOfSamples ) )
        )
    {
        isUniformRandomNumberGeneratorErroneous = true;

        cerr << "The computed sample mean ( " << sampleMean
             << " ) using a sample size of " << numberOfSamples
             << " is not within the standard error of the true sample mean"
             << " ( 1/2 ) " << endl;
        cerr << "The distribution is not uniform to a confidence of 99.96%. "
             << "Run the unit test again to ensure error is not due to "
             << "statistical anomaly. " << endl;
    }

    if ( ( 1.0 / 12.0 ) / sampleVariance
         < static_cast< double >( numberOfSamples - 1 ) / 101590
         || ( 1.0 / 12.0 ) / sampleVariance
         > static_cast< double>( numberOfSamples - 1 ) / 98424 )
    {
        isUniformRandomNumberGeneratorErroneous = true;

        cerr << "The computed sample variance ( " << sampleVariance
             << " ) using a sample size of " << numberOfSamples
             << " is not within the standard error of the true sample variance"
             << "( 1/12 ) " << endl;
        cerr << "The distribution is not uniform to a confidence of 99.96%. "
             << "Run the unit test again to ensure error is not due to "
             << "statistical anomaly." << endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isUniformRandomNumberGeneratorErroneous;
}

//! Test implementation of exponential random number generator class.
bool testExponentialRandomNumberGenerator( )
{
    // Summary of tests.
    // Test 1: Get a normalized random double in the interval [ 0, infinity ).
    // Test 2: Compute sample mean and standard distrubution of distribution
    //         and compare to analytical results.

    // Test result initialised to false.
    bool isExponentialRandomNumberGeneratorErroneous = false;

    // Create exponential random number generator object.
    // Exponential parameter is set to 2.0.
    ExponentialRandomNumberGenerator exponentialRandomNumbers( 2.0,
                                                               time( NULL ) );

    // Test 1: Get a normalized random double in the interval [ 0, infinity ).
    double computedResultForTest1 = exponentialRandomNumbers
            .getExponentiallyDistributedNormalizedRandomDouble( );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( computedResultForTest1 < 0.0  )
    {
        isExponentialRandomNumberGeneratorErroneous = true;

        cerr << "The computed value ( " << computedResultForTest1
             << " ) using the normalized exponential random number generator "
             << "is not within the prescribed bounds [ 0, infinity )" << endl;
    }

    // Test 2: Compute sample mean and standard distrubution of distribution
    //         and compare to analytical results.
    // Declare and initialize number of samples.
    unsigned int numberOfSamples = 100000;

    // Declare and set size of vector of sample exponential random values.
    vector< double > sampleOfExponentialRandomValues;

    // Fill vector of sample exponential random values.
    for ( unsigned int i = 0; i < numberOfSamples; i++ )
    {
        sampleOfExponentialRandomValues.push_back(
                    exponentialRandomNumbers
                    .getExponentiallyDistributedNormalizedRandomDouble( ) );
    }

    // Estimate sample mean of the distribution.
    double sampleMean = computeSampleMean( sampleOfExponentialRandomValues );

    // Estimate sample variance.
    double sampleVariance = computeSampleVariance(
                sampleOfExponentialRandomValues );

    // Compute differences between computed and expected results.
    // Test sample mean versus true mean with a confidence of 99.96%.
    // This means that the test below guarentees that this unit test will
    // fail on average every 2500 times it is run.
    // Also test sample variance versus true variance with a confidence of
    // 99.96%. The true mean and true variance of an exponential distribution
    // with parameter equal to 2.0 are 1/2 and 1/4 respectively.
    // The tolerances were computed using confidence intervals for the sample
    // mean and variance of a normal distribution and table lookup.
    if ( computeAbsoluteValue( sampleMean - 0.5 )
         > 3.49 * 1.0 / sqrt( 4.0 )
           * 1.0 / sqrt( static_cast< double >( numberOfSamples ) )
        )
    {
        isExponentialRandomNumberGeneratorErroneous = true;

        cerr << "The computed sample mean ( " << sampleMean
             << " ) using a sample size of " << numberOfSamples
             << " is not within the standard error of the true sample mean"
             << " ( 1/2 ) " << endl;
        cerr << "The distribution is not exponential to a confidence of "
             << "99.96%. Run the unit test again to ensure error is not due "
             << "to statistical anomaly. " << endl;
    }

    if ( ( 1.0 / 4.0 ) / sampleVariance
         < static_cast< double >( numberOfSamples - 1 ) / 101590
         || ( 1.0 / 4.0 ) / sampleVariance
         > static_cast< double>( numberOfSamples - 1 ) / 98424 )
    {
        isExponentialRandomNumberGeneratorErroneous = true;

        cerr << "The computed sample variance ( " << sampleVariance
             << " ) using a sample size of " << numberOfSamples
             << " is not within the standard error of the true sample variance"
             << "( 1/4 ) " << endl;
        cerr << "The distribution is not exponential to a confidence of "
             << "99.96%. Run the unit test again to ensure error is not due "
             << "to statistical anomaly." << endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isExponentialRandomNumberGeneratorErroneous;
}

//! Test implementation of normal random number generator class.
bool testNormalRandomNumberGenerator( )
{
    // Summary of tests.
    // Test 1: Compute sample mean and standard deviation of distribution
    //         and compare to analytical results.

    // Test result initialized to false.
    bool isNormalRandomNumberGeneratorErroneous = false;

    // Create normal random number generator object.
    // Normal distribution mean set to 1.0, standard deviation set to 2.0.
    NormalRandomNumberGenerator normalRandomNumbers( 1.0, 2.0, time( NULL ) );

    // Test 1: Compute sample mean and standard deviation of distribution
    //         and compare to analytical results.
    // Declare and initialize number of samples.
    unsigned int numberOfSamples = 100000;

    // Declare and set size of vector of sample normal random values.
    vector< double > sampleOfNormalRandomValues;

    // Fill vector of sample normal random values.
    for ( unsigned int i = 0; i < numberOfSamples; i++ )
    {
        sampleOfNormalRandomValues.push_back(
                    normalRandomNumbers
                    .getNormallyDistributedNormalizedRandomDouble( ) );
    }

    // Estimate sample mean of the distribution.
    double sampleMean = computeSampleMean( sampleOfNormalRandomValues );

    // Estimate sample variance.
    double sampleVariance = computeSampleVariance(
                sampleOfNormalRandomValues );

    // Compute differences between computed and expected results.
    // Test sample mean versus true mean with a confidence of 99.96%.
    // This means that the test below guarentees that this unit test will
    // fail on average every 2500 times it is run.
    // Also test sample variance versus true variance with a confidence of
    // 99.96%. The true mean and true variance of the tested normal distribution
    // distribution is equal to 1 and 4 respectively.
    // The tolerances were computed using confidence intervals for the sample
    // mean and variance of a normal distribution and table lookup.
    if ( computeAbsoluteValue( sampleMean - 1.0 )
         > 3.49 * 4.0
           * 1.0 / sqrt( static_cast< double >( numberOfSamples ) )
        )
    {
        isNormalRandomNumberGeneratorErroneous = true;

        cerr << "The computed sample mean ( " << sampleMean
             << " ) using a sample size of " << numberOfSamples
             << " is not within the standard error of the true sample mean"
             << " ( 1/2 ) " << endl;
        cerr << "The distribution is not normal to a confidence of "
             << "99.96%. Run the unit test again to ensure error is not due "
             << "to statistical anomaly. " << endl;
    }

    if ( 4.0 / sampleVariance
         < static_cast< double >( numberOfSamples - 1 ) / 101590
         || 4.0 / sampleVariance
         > static_cast< double>( numberOfSamples - 1 ) / 98424 )
    {
        isNormalRandomNumberGeneratorErroneous = true;

        cerr << "The computed sample variance ( " << sampleVariance
             << " ) using a sample size of " << numberOfSamples
             << " is not within the standard error of the true sample variance"
             << "( 4 ) " << endl;
        cerr << "The distribution is not normal to a confidence of "
             << "99.96%. Run the unit test again to ensure error is not due "
             << "to statistical anomaly." << endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isNormalRandomNumberGeneratorErroneous;
}

}

// End of file.
