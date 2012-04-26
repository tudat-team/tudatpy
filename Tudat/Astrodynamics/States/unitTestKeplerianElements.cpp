/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      110110    K. Kumar          File created.
 *      110121    K. Kumar          Updated to comply with new protocol.
 *      110310    K. Kumar          Changed right ascension of ascending node to longitude of
 *                                  ascending node.
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// Test runs code and verifies result against expected value.
// If the tested code is erroneous, the test function returns a boolean
// true; if the code is correct, the function returns a boolean false.
// 

#include <iostream>
#include <limits>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include "Tudat/Astrodynamics/States/keplerianElements.h"

//! Test Keplerian elements state class.
int main( )
{
    // Using declarations.
    using std::cerr;
    using std::endl;
    using tudat::unit_conversions::convertDegreesToRadians;
    using namespace tudat;

    // Three tests.
    // Test 1: Set Keplerian elements using individual set functions; and get
    //         Keplerian elements using individual get functions.
    // Test 2: Set Keplerian elements using setState( ) function; and get
    //         Keplerian elements using getState( ) function.
    // Test 3: Get auxilliary Keplerian elements.

    // Initialize unit test result to false.
    bool isKeplerianElementsErroneous = false;

    // Create Keplerian elements state objects.
    using tudat::astrodynamics::states::KeplerianElements;
    KeplerianElements keplerianElementsStateTest1_;
    KeplerianElements keplerianElementsStateTest2_;

    // Expected test result.
    // Create vector of Keplerian elements: semi-major axis, eccentricity,
    // inclination, argument of periapsis, longitude of the ascending node,
    // true anomaly.
    Eigen::VectorXd vectorOfKeplerianElements_( 6 );
    vectorOfKeplerianElements_( 0 ) = 2.5e6;
    vectorOfKeplerianElements_( 1 ) = 0.1;
    vectorOfKeplerianElements_( 2 ) = convertDegreesToRadians( 102.3 );
    vectorOfKeplerianElements_( 3 ) = convertDegreesToRadians( 125.7 );
    vectorOfKeplerianElements_( 4 ) = convertDegreesToRadians( 215.34 );
    vectorOfKeplerianElements_( 5 ) = convertDegreesToRadians( 0.0 );

    // Create vector of auxilliary Keplerian elements: longitude of periapsis,
    // true longitude, argument of latitude, semi-latus rectum.
    Eigen::VectorXd vectorOfAuxilliaryKeplerianElements_( 4 );
    vectorOfAuxilliaryKeplerianElements_( 0 ) = vectorOfKeplerianElements_( 3 )
            + vectorOfKeplerianElements_( 4 );
    vectorOfAuxilliaryKeplerianElements_( 1 ) = vectorOfKeplerianElements_( 3 )
            + vectorOfKeplerianElements_( 4 ) + vectorOfKeplerianElements_( 5 );
    vectorOfAuxilliaryKeplerianElements_( 2 )
            = vectorOfKeplerianElements_( 3 ) + vectorOfKeplerianElements_( 5 );
    vectorOfAuxilliaryKeplerianElements_( 3 ) = 2.0e6;

    // Test 1: Set Keplerian elements in state object.
    keplerianElementsStateTest1_.setSemiMajorAxis( vectorOfKeplerianElements_( 0 ) );
    keplerianElementsStateTest1_.setEccentricity( vectorOfKeplerianElements_( 1 ) );
    keplerianElementsStateTest1_.setInclination( vectorOfKeplerianElements_( 2 ) );
    keplerianElementsStateTest1_.setArgumentOfPeriapsis( vectorOfKeplerianElements_( 3 ) );
    keplerianElementsStateTest1_.setLongitudeOfAscendingNode( vectorOfKeplerianElements_( 4 ) );
    keplerianElementsStateTest1_.setTrueAnomaly( vectorOfKeplerianElements_( 5 ) );

    // Test 1: Get Keplerian elements and store in a state vector.
    Eigen::VectorXd keplerianElementsStateVectorTest1( 6 );
    keplerianElementsStateVectorTest1( 0 ) = keplerianElementsStateTest1_.getSemiMajorAxis( );
    keplerianElementsStateVectorTest1( 1 ) = keplerianElementsStateTest1_.getEccentricity( );
    keplerianElementsStateVectorTest1( 2 ) = keplerianElementsStateTest1_.getInclination( );
    keplerianElementsStateVectorTest1( 3 )
            = keplerianElementsStateTest1_.getArgumentOfPeriapsis( );
    keplerianElementsStateVectorTest1( 4 )
            = keplerianElementsStateTest1_.getLongitudeOfAscendingNode( );
    keplerianElementsStateVectorTest1( 5 ) = keplerianElementsStateTest1_.getTrueAnomaly( );

    // Test 1: Difference between setting each Keplerian element and the
    // expected values.
    Eigen::VectorXd differenceBetweenResultsTest1_( 6 );
    differenceBetweenResultsTest1_
            = keplerianElementsStateVectorTest1 - vectorOfKeplerianElements_;

    // Test 2: Set Keplerian elements using setState( ) function.
    keplerianElementsStateTest2_.state = vectorOfKeplerianElements_;

    // Test 2: Difference between setting the Keplerian state as a whole and
    // the expected values.
    Eigen::VectorXd differenceBetweenResultsTest2_;
    differenceBetweenResultsTest2_
            = keplerianElementsStateTest2_.state - vectorOfKeplerianElements_;

    // Test 3: Set eccentricity and semi-latus rectum ( only for parabolic
    // orbits ).
    keplerianElementsStateTest2_.setEccentricity( 1.0 );
    keplerianElementsStateTest2_.setSemiLatusRectum(
                vectorOfAuxilliaryKeplerianElements_( 3 ) );

    // Test 3: Get auxilliary Keplerian elements and store in a vector.
    Eigen::VectorXd auxilliaryKeplerianElementsVectorTest3_( 4 );
    auxilliaryKeplerianElementsVectorTest3_( 0 )
            = keplerianElementsStateTest2_.getLongitudeOfPeriapsis( );
    auxilliaryKeplerianElementsVectorTest3_( 1 )
            = keplerianElementsStateTest2_.getTrueLongitude( );
    auxilliaryKeplerianElementsVectorTest3_( 2 )
            = keplerianElementsStateTest2_.getArgumentOfLatitude( );
    auxilliaryKeplerianElementsVectorTest3_( 3 )
            = keplerianElementsStateTest2_.getSemiLatusRectum( );

    // Test 3: Difference between getting each auxilliary Keplerian element and
    // the expected values.
    Eigen::VectorXd differenceBetweenResultsTest3( 4 );
    differenceBetweenResultsTest3 = vectorOfAuxilliaryKeplerianElements_
            - auxilliaryKeplerianElementsVectorTest3_;

    // Set test result to true if the test does not match the expected result.
    if ( differenceBetweenResultsTest1_.norm( ) > std::numeric_limits< double >::epsilon( ) )
    {
        isKeplerianElementsErroneous = true;

        // Generate error statements.
        cerr << "The computed value " << endl;
        cerr << "( " << keplerianElementsStateVectorTest1 << " )" << endl;
        cerr << "using indivudual set/get functions for the "
             << "KeplerianElements class does not match the expected solution " << endl;
        cerr << "( " << vectorOfKeplerianElements_ << " )." << endl;
        cerr << "The difference is: "
             << vectorOfKeplerianElements_ - keplerianElementsStateVectorTest1 << endl;
    }

    if ( differenceBetweenResultsTest2_.norm( ) > std::numeric_limits< double >::epsilon( ) )
    {
        isKeplerianElementsErroneous = true;

        // Generate error statements.
        cerr << "The computed value " << endl;
        cerr << "( " << keplerianElementsStateTest2_.state << " )" << endl;
        cerr << "using set/get functions for the whole state for the "
             << "KeplerianElements class does not match the expected solution " << endl;
        cerr << "( " << vectorOfKeplerianElements_ << " )." << endl;
        cerr << "The difference is: "
             << vectorOfKeplerianElements_ - keplerianElementsStateTest2_.state << endl;
    }

    if ( differenceBetweenResultsTest3.norm( ) > std::numeric_limits< double >::epsilon( ) )
    {
        isKeplerianElementsErroneous = true;

        // Generate error statements.
        cerr << "The computed value " << endl;
        cerr << "( " << auxilliaryKeplerianElementsVectorTest3_ << " )" << endl;
        cerr << "using indivudual set/get functions for the auxilliary elements of the "
             << "KeplerianElements class does not match the expected solution " << endl;
        cerr << "( " << vectorOfAuxilliaryKeplerianElements_ << " )." << endl;
        cerr << "The difference is: "
             << vectorOfAuxilliaryKeplerianElements_ - auxilliaryKeplerianElementsVectorTest3_
             << endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isKeplerianElementsErroneous )
    {
        cerr << "testKeplerianElements failed!" << endl;
    }

    return isKeplerianElementsErroneous;
}
