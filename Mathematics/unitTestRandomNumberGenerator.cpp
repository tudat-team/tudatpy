/*! \file unitTestRandomNumberGenerator.cpp
 *    Source file that defines a unit test that tests the random number
 *    generator included in Tudat.
 *
 *    Path              : /Mathematics/
 *    Version           : 2
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
 *    Last modified     : 7 Feburary, 2011
 *
 *    References
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
 *
 *      It does not test whether the uniform distribution between 0 and 1
 *      has mean 0.5 and variance 1/12. This should be added in the future.
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
 */

// Include statements.
#include "randomNumberGenerator.h"

// Using declarations.
using mathematics::computeAbsoluteValue;
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

    // Test result initialised to false.
    bool isRandomNumberGeneratorErroneous = false;

    // Create random number generator object.
    RandomNumberGenerator randomNumbers( time( NULL ) );

    // Results computed using implementation of random number generator class.
    int computedResultForTest1 = randomNumbers
                                 .getUniformlyDistributedRandom64BitInteger( );
    int computedResultForTest2 = randomNumbers
                                 .getUniformlyDistributedRandom32BitInteger( );
    double computedResultForTest3
            = randomNumbers.getUniformlyDistributedNormalizedRandomDouble( );
    double computedResultForTest4 = randomNumbers.getRandomPlusMinusSign( );

    // Compute differences between computed and expected results.
    Vector4d differenceBetweenResults;
    differenceBetweenResults.setConstant( 1.0 );

    if ( fmod( computedResultForTest1, floor( computedResultForTest1 ) )
         < MACHINE_PRECISION_DOUBLES )
    {
        differenceBetweenResults( 0 ) = 0.0;
    }

    if ( fmod( computedResultForTest2, floor( computedResultForTest2 ) )
         < MACHINE_PRECISION_DOUBLES )
    {
        differenceBetweenResults( 1 ) = 0.0;
    }

    if ( computedResultForTest3 > 0.0 && computedResultForTest3 < 1.0 )
    {
        differenceBetweenResults( 2 ) = 0.0;
    }

    if ( computeAbsoluteValue( computedResultForTest4 ) - 1.0
         < MACHINE_PRECISION_DOUBLES )
    {
        differenceBetweenResults( 3 ) = 0.0;
    }

    // Set test result to false if the test does not match the expected result.
    if ( computeAbsoluteValue( differenceBetweenResults.norm( ) )
        > MACHINE_PRECISION_DOUBLES )
    {
        isRandomNumberGeneratorErroneous = true;

        // Generate cerr statements.
        if ( computeAbsoluteValue( differenceBetweenResults( 0 ) ) > 0.0 )
        {
            cerr << "The computed value ( " << computedResultForTest1
                 << " ) using the 64-bit random integer generator "
                 << " is not an integer" << endl;
        }

        if ( computeAbsoluteValue( differenceBetweenResults( 1 ) ) > 0.0 )
        {
            cerr << "The computed value ( " << computedResultForTest2
                 << " ) using the 32-bit random integer generator "
                 << " is not an integer" << endl;
        }

        if ( computeAbsoluteValue( differenceBetweenResults( 2 ) ) > 0.0 )
        {
            cerr << "The computed value ( " << computedResultForTest3
                 << " ) using the normalized random number generator "
                 << " is not within the prescribed bounds [ 0, 1 ]" << endl;
        }

        if ( computeAbsoluteValue( differenceBetweenResults( 3 ) ) > 0.0 )
        {
            cerr << "The computed value ( " << computedResultForTest4
                 << " ) using the plus/minus sign random generator "
                 << " is not equal to either +1 or -1" << endl;
        }
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isRandomNumberGeneratorErroneous;
}

}

// End of file.
