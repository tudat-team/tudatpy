/*! \file unitTestKeplerianElements.cpp
 *    Source file for a unit test that tests the implementation of the
 *    Keplerian elements state class in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 2
 *    Check status      : Unchecked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Date created      : 10 January, 2010
 *    Last modified     : 21 January, 2010
 *
 *    References
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
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
 *      YYMMDD    author              comment
 *      110110    K. Kumar            File created.
 *      110121    K. Kumar            Updated to comply with new protocol.
 */

// Include statements.
#include "unitTestKeplerianElements.h"

// Using statements.
using std::cerr;
using std::endl;

//! Namespace for all unit tests.
namespace unit_tests
{

bool testKeplerianElements( )
{
    // Three tests.
    // Test 1: Set Keplerian elements using individual set functions; and get
    //         Keplerian elements using individual get functions.
    // Test 2: Set Keplerian elements using setState() function; and get
    //         Keplerian elements using getState() function.
    // Test 3: Get auxilliary Keplerian elements.

    // Test result initialised to true.
    bool isKeplerianElements = false;

    // Create Keplerian elements state objects.
    KeplerianElements keplerianElementsStateTest1;
    KeplerianElements keplerianElementsStateTest2;

    // Expected test result.
    // Create vector of Keplerian elements: semi-major axis, eccentricity,
    // inclination, argument of periapsis, right ascension of the ascending
    // node, true anomaly.
    VectorXd vectorOfKeplerianElements( 6 );
    vectorOfKeplerianElements( 0 ) = 2.5e6;
    vectorOfKeplerianElements( 1 ) = 0.1;
    vectorOfKeplerianElements( 2 ) = convertDegreesToRadians( 102.3 );
    vectorOfKeplerianElements( 3 ) = convertDegreesToRadians( 125.7 );
    vectorOfKeplerianElements( 4 ) = convertDegreesToRadians( 215.34 );
    vectorOfKeplerianElements( 5 ) = convertDegreesToRadians( 0.0 );

    // Create vector of auxilliary Keplerian elements: longitude of periapsis,
    // true longitude, argument of latitude, semi-latus rectum.
    VectorXd vectorOfAuxilliaryKeplerianElements( 4 );
    vectorOfAuxilliaryKeplerianElements( 0 ) = vectorOfKeplerianElements( 3 )
                                               + vectorOfKeplerianElements( 4 );
    vectorOfAuxilliaryKeplerianElements( 1 ) = vectorOfKeplerianElements( 3 )
                                               + vectorOfKeplerianElements( 4 )
                                               + vectorOfKeplerianElements( 5 );
    vectorOfAuxilliaryKeplerianElements( 2 ) = vectorOfKeplerianElements( 3 )
                                               + vectorOfKeplerianElements( 5 );
    vectorOfAuxilliaryKeplerianElements( 3 ) = 2.0e6;

    // Test 1: Set Keplerian elements in state object.
    keplerianElementsStateTest1.setSemiMajorAxis(
            vectorOfKeplerianElements( 0 ) );
    keplerianElementsStateTest1.setEccentricity(
            vectorOfKeplerianElements( 1 ) );
    keplerianElementsStateTest1.setInclination(
            vectorOfKeplerianElements( 2 ) );
    keplerianElementsStateTest1.setArgumentOfPeriapsis(
            vectorOfKeplerianElements( 3 ) );
    keplerianElementsStateTest1.setRightAscensionOfAscendingNode(
            vectorOfKeplerianElements( 4 ) );
    keplerianElementsStateTest1.setTrueAnomaly(
            vectorOfKeplerianElements( 5 ) );

    // Test 1: Get Keplerian elements and store in a state vector.
    VectorXd keplerianElementsStateVectorTest1( 6 );
    keplerianElementsStateVectorTest1( 0 ) = keplerianElementsStateTest1
                                             .getSemiMajorAxis( );
    keplerianElementsStateVectorTest1( 1 ) = keplerianElementsStateTest1
                                             .getEccentricity( );
    keplerianElementsStateVectorTest1( 2 ) = keplerianElementsStateTest1
                                             .getInclination( );
    keplerianElementsStateVectorTest1( 3 ) = keplerianElementsStateTest1
                                             .getArgumentOfPeriapsis( );
    keplerianElementsStateVectorTest1( 4 ) = keplerianElementsStateTest1
                                          .getRightAscensionOfAscendingNode( );
    keplerianElementsStateVectorTest1( 5 ) = keplerianElementsStateTest1
                                             .getTrueAnomaly( );

    // Test 1: Difference between setting each Keplerian element and the
    // expected values.
    VectorXd differenceBetweenResultsTest1( 6 );
    differenceBetweenResultsTest1 =  keplerianElementsStateVectorTest1
                                     - vectorOfKeplerianElements;

    // Test 2: Set Keplerian elements using setState() function.
    keplerianElementsStateTest2.setState( vectorOfKeplerianElements );

    // Test 2: Difference between setting the Keplerian state as a whole and
    // the expected values.
    VectorXd differenceBetweenResultsTest2;
    differenceBetweenResultsTest2 = keplerianElementsStateTest2.getState( )
                                    - vectorOfKeplerianElements;

    // Test 3: Set eccentricity and semi-latus rectum ( only for parabolic
    // orbits ).
    keplerianElementsStateTest2.setEccentricity( 1.0 );
    keplerianElementsStateTest2.setSemiLatusRectum(
            vectorOfAuxilliaryKeplerianElements( 3 ) );

    // Test 3: Get auxilliary Keplerian elements and store in a vector.
    VectorXd auxilliaryKeplerianElementsVectorTest3( 4 );
    auxilliaryKeplerianElementsVectorTest3( 0 ) = keplerianElementsStateTest2
                                                  .getLongitudeOfPeriapsis( );
    auxilliaryKeplerianElementsVectorTest3( 1 ) = keplerianElementsStateTest2
                                                  .getTrueLongitude( );
    auxilliaryKeplerianElementsVectorTest3( 2 ) = keplerianElementsStateTest2
                                                  .getArgumentOfLatitude( );
    auxilliaryKeplerianElementsVectorTest3( 3 ) = keplerianElementsStateTest2
                                                  .getSemiLatusRectum( );

    // Test 3: Difference between getting each auxilliary Keplerian element and
    // the expected values.
    VectorXd differenceBetweenResultsTest3( 4 );
    differenceBetweenResultsTest3 = vectorOfAuxilliaryKeplerianElements
                                    - auxilliaryKeplerianElementsVectorTest3;

    // Set test result to true if the test does not match the expected result.
    if ( differenceBetweenResultsTest1.norm( )
         > mathematics::MACHINE_PRECISION_DOUBLES
         || differenceBetweenResultsTest2.norm( )
         > mathematics::MACHINE_PRECISION_DOUBLES
         || differenceBetweenResultsTest3.norm( )
         > mathematics::MACHINE_PRECISION_DOUBLES )
    {
        isKeplerianElements = true;

        if ( differenceBetweenResultsTest1.norm( )
             > mathematics::MACHINE_PRECISION_DOUBLES )
        {
            // Generate error statements.
            cerr << "The computed value " << endl;
            cerr << "( " << keplerianElementsStateVectorTest1 << " )"
                 << endl;
            cerr << "using indivudual set/get functions for the "
                 << "KeplerianElements class does not match the expected "
                 << "solution " << endl;
            cerr << "( " << vectorOfKeplerianElements << " )." << endl;
            cerr << "The difference is: "
                 << vectorOfKeplerianElements
                    - keplerianElementsStateVectorTest1 << endl;
        }

        if ( differenceBetweenResultsTest2.norm( )
             > mathematics::MACHINE_PRECISION_DOUBLES )
        {
            // Generate error statements.
            cerr << "The computed value " << endl;
            cerr << "( " << keplerianElementsStateTest2.getState( ) << " )"
                 << endl;
            cerr << "using set/get functions for the whole state for the "
                 << "KeplerianElements class does not match the expected "
                 << "solution " << endl;
            cerr << "( " << vectorOfKeplerianElements << " )." << endl;
            cerr << "The difference is: "
                 << vectorOfKeplerianElements
                    - keplerianElementsStateTest2.getState( ) << endl;
        }

        if ( differenceBetweenResultsTest3.norm( )
             > mathematics::MACHINE_PRECISION_DOUBLES )
        {
            // Generate error statements.
            cerr << "The computed value " << endl;
            cerr << "( " << auxilliaryKeplerianElementsVectorTest3 << " )"
                 << endl;
            cerr << "using indivudual set/get functions for the "
                 << "auxilliary elements of the KeplerianElements class does "
                 << "not match the expected solution " << endl;
            cerr << "( " << vectorOfAuxilliaryKeplerianElements << " )."
                    << endl;
            cerr << "The difference is: "
                 << vectorOfAuxilliaryKeplerianElements
                    - auxilliaryKeplerianElementsVectorTest3 << endl;
        }
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isKeplerianElements;
}

}

// End of file.
