/*! \file unitTestRandomNumberGenerator.cpp
 *    Source file that defines a unit test that tests the random number
 *    generator included in Tudat.
 *
 *    Path              : /Mathematics/
 *    Version           : 3
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
 *    Date created      : 7 January, 2011
 *    Last modified     : 7 July, 2011
 *
 *    References
 *      Walpole, R.E., et al. Probabilty and Statistics for Engineers &
 *          Scientists, Seventh Edition, Prentice Hall, NJ, 2002.
 *      Texas A&M University. Chi Square Calculator,
 *          http://www.stat.tamu.edu/~west/applets/chisqdemo.html, last
 *          accessed: 7 July, 2011.
 *      Faculty of Mathematics and Informatics, University of Sofia.
 *          http://www.fmi.uni-sofia.bg/vesta/virtual_labs/interval/interval6.html,
 *          last accessed: 7 July, 2011.
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
 *
 *      The code tests for the quality of the uniform distribution by testing
 *      the sample mean and variance against the true mean and variance using
 *      standard error estimates. These estimates are stated with a
 *      confidence of 99.96%. If the unit test fails due to Test 5, it is
 *      recommended that the unit test is run again to rule out statistical
 *      anomalies.
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
 */

// Include statements.
#include "randomNumberGenerator.h"

// Using declarations.
using mathematics::computeAbsoluteValue;
using mathematics::raiseToIntegerPower;
using mathematics::MACHINE_PRECISION_DOUBLES;
using std::cerr;
using std::endl;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test of implementation of random number generator class.
bool testRandomNumberGenerator( )
{
    // Four tests.
    // Test 1: Get a random integer ( 64-bit ).
    // Test 2: Get a random integer ( 32-bit ).
    // Test 3: Get a normalized random double in the interval [ 0, 1 ].
    // Test 4: Get a random plus/minus sign.
    // Test 5: Compute sample mean and standard distrubution of distribution
    //         and compare to analytical results.

    // Test result initialised to false.
    bool isRandomNumberGeneratorErroneous = false;

    // Create random number generator object.
    RandomNumberGenerator randomNumbers( time( NULL ) );

    // Results computed using implementation of random number generator class.
    // Test 1: Get a random integer ( 64-bit ).
    int computedResultForTest1 = randomNumbers
                                 .getUniformlyDistributedRandom64BitInteger( );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( fmod( computedResultForTest1, floor( computedResultForTest1 ) )
         > MACHINE_PRECISION_DOUBLES )
    {
        isRandomNumberGeneratorErroneous = true;

        cerr << "The computed value ( " << computedResultForTest1
                << " ) using the 64-bit random integer generator "
                << "is not an integer" << endl;
    }

    // Test 2: Get a random integer ( 32-bit ).
    int computedResultForTest2 = randomNumbers
                                 .getUniformlyDistributedRandom32BitInteger( );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( fmod( computedResultForTest2, floor( computedResultForTest2 ) )
         > MACHINE_PRECISION_DOUBLES )
    {
        isRandomNumberGeneratorErroneous = true;

        cerr << "The computed value ( " << computedResultForTest2
                << " ) using the 32-bit random integer generator "
                << "is not an integer" << endl;
    }

    // Test 3: Get a normalized random double in the interval [ 0, 1 ].
    double computedResultForTest3
            = randomNumbers.getUniformlyDistributedNormalizedRandomDouble( );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( computedResultForTest3 < 0.0 && computedResultForTest3 > 1.0 )
    {
        isRandomNumberGeneratorErroneous = true;

        cerr << "The computed value ( " << computedResultForTest3
                << " ) using the normalized random number generator "
                << "is not within the prescribed bounds [ 0, 1 ]" << endl;
    }

    // Test 4: Get a random plus/minus sign.
    double computedResultForTest4 = randomNumbers.getRandomPlusMinusSign( );

    // Compute differences between computed and expected results and generate
    // cerr statement if test fails.
    if ( computeAbsoluteValue( computedResultForTest4 ) - 1.0
         > MACHINE_PRECISION_DOUBLES )
    {
        isRandomNumberGeneratorErroneous = true;

        cerr << "The computed value ( " << computedResultForTest4
                << " ) using the plus/minus sign random generator "
                << "is not equal to either +1 or -1" << endl;
    }

    // Test 5: Compute sample mean and variance of distribution and compare to
    //         analytical results.
    // Declare and initialize number of samples.
    unsigned int numberOfSamples = 100000;

    // Declare and set size of vector of sample random values.
    VectorXd sampleOfRandomValues = VectorXd( numberOfSamples );

    // Fill vector of sample random values.
    for ( unsigned int i = 0; i < numberOfSamples; i++ )
    {
        sampleOfRandomValues( i )
                = randomNumbers
                  .getUniformlyDistributedNormalizedRandomDouble( );
    }

    // Estimate sample mean of the distribution.
    double sampleMean = sampleOfRandomValues.sum( )
                        / numberOfSamples;

    // Estimate sample variance.
    double sumOfResidualsSquared = 0.0;
    for ( unsigned int i = 0; i < numberOfSamples; i++ )
    {
        sumOfResidualsSquared +=
                raiseToIntegerPower( sampleOfRandomValues( i )
                                     - sampleMean, 2 );
    }
    double sampleVariance
            = 1.0 / ( numberOfSamples - 1.0 ) * sumOfResidualsSquared;

    // Compute differences between computed and expected results.
    // Test sample mean versus true mean with a confidence of 99.96%.
    // This means that the test below guarentees that this unit test will
    // fail on average every 2500 times it is run.
    // Also test sample variance versus true variance with a confidence of
    // 99.96%. The true mean and true variance of a uniform distribution
    // from 0 to 1 are 1/2 and 1/12 respectively.
    // The tolerances were computed using confidence intervals for the sample
    // mean and variance and table lookup.
    if ( computeAbsoluteValue( sampleMean - 0.5 )
         > 3.49 * 1.0 / sqrt( 12.0 )
           * 1.0 / static_cast< double >( sqrt( numberOfSamples ) )
        )
    {
        isRandomNumberGeneratorErroneous = true;

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
        isRandomNumberGeneratorErroneous = true;

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
    return isRandomNumberGeneratorErroneous;
}

}

// End of file.
