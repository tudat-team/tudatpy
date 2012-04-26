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
 *      101215    K. Kumar          First creation of code.
 *      101216    K. Kumar          Updated to include test of getPotential( ).
 *      101230    K. Kumar          Updated to include test of getGradient( ) and getLaplacian( ).
 *      110104    J. Melman         Some minor comment and layout changes.
 *      110106    K. Kumar          Updated test using predefined gravity field and added machine
 *                                  precision variable.
 *      110107    K. Kumar          Updated call to predefined gravity field. Updated unit test to
 *                                  new protocol, added namespace, new filename, new file location.
 *      110113    K. Kumar          Added cerr statements.
 *      110115    J. Melman         Changed the error messages.
 *      110128    K. Kumar          Updated code to work with pointers.
 *      110202    K. Kumar          Updated code to work with State.
 *      110204    K. Kumar          Removed "vector" from naming.
 *      110310    K. Kumar          Changed naming from Laplacian to gradient tensor.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// Test runs code and verifies result against expected value.
// If the tested code is erroneous, the test function returns a boolean
// true; if the code is correct, the function returns a boolean false.
// 


#include <cmath>
#include <limits>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Bodies/planet.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"

//! Test implementation of spherical harmonics gravity field class.
int main( )
{
    using std::cerr;
    using std::endl;
    using std::pow;
    using std::fabs;
    using namespace tudat;

    // Summary of tests.
    // Test 1: Test setting and getting gravitational parameter.
    // Test 2: Test getting gravitational parameter for predefined Earth central body gravity
    //         field.
    // Test 3: Test getting potential given state.
    // Test 4: Test getting gradient of potential given state.
    // Test 5: Test getting gradient tensor of potential given state.

    // Test result initialised to false.
    bool isSphericalHarmonicsGravityFieldErroneous = false;

    // Set identity matrix.
    Eigen::Matrix3d identityMatrix;
    identityMatrix.setIdentity( 3, 3 );

    // Create gravity field for myPlanet.
    using astrodynamics::gravitation::SphericalHarmonicsGravityField;
    SphericalHarmonicsGravityField myPlanetGravityField;

    // Set gravitational parameter of myPlanet.
    double gravitationalParameterOfMyPlanet = 22032.00;
    myPlanetGravityField.setGravitationalParameter( gravitationalParameterOfMyPlanet );

    // Set origin of gravity field of myPlanet with respect to geometric
    // center.
    myPlanetGravityField.setOrigin( Eigen::Vector3d::Zero( ) );

    // Set position with respect to geometric center.
    Eigen::Vector3d cartesianPosition;
    cartesianPosition << 5.0e6, 3.0e6, 1.0e6;

    // Create predefined Earth central gravity field.
    using astrodynamics::gravitation::CentralGravityField;
    CentralGravityField predefinedEarthCentralGravityField( CentralGravityField::earth );

    // Expected test results.
    double expectedResultForTest1 = gravitationalParameterOfMyPlanet;
    double expectedResultForTest2 = 3.9859383624e+14;
    double expectedResultForTest3 = gravitationalParameterOfMyPlanet / cartesianPosition.norm( );
    Eigen::VectorXd expectedResultForTest4 =  -gravitationalParameterOfMyPlanet
            / pow( cartesianPosition.norm( ), 3.0 ) * cartesianPosition;
    Eigen::Matrix3d expectedResultForTest5 = gravitationalParameterOfMyPlanet
            / pow( cartesianPosition.norm( ), 5.0 )
            * ( ( 3.0 * cartesianPosition * cartesianPosition.transpose( ) )
                - ( cartesianPosition.squaredNorm( ) * identityMatrix ) );

    // Results computed using implementation of spherical harmonics gravity
    // field class.
    double computedResultForTest1 = myPlanetGravityField.getGravitationalParameter( );
    double computedResultForTest2 =
            predefinedEarthCentralGravityField.getGravitationalParameter( );
    double computedResultForTest3 = myPlanetGravityField.getPotential(
                cartesianPosition );
    Eigen::Vector3d computedResultForTest4 = myPlanetGravityField.getGradientOfPotential(
                cartesianPosition );
    Eigen::Matrix3d computedResultForTest5 = myPlanetGravityField.getGradientTensorOfPotential(
                cartesianPosition );

    // Compute differences between computed and expected results.
    Eigen::VectorXd differenceBetweenResults( 5 );
    Eigen::Vector3d differenceBetweenResultsForTest4;
    Eigen::Matrix3d differenceBetweenResultsForTest5;
    differenceBetweenResults( 0 ) = fabs( computedResultForTest1 - expectedResultForTest1 );
    differenceBetweenResults( 1 ) = fabs( computedResultForTest2 - expectedResultForTest2 );
    differenceBetweenResults( 2 ) = fabs( computedResultForTest3 - expectedResultForTest3 );
    differenceBetweenResultsForTest4 = computedResultForTest4 - expectedResultForTest4;
    differenceBetweenResults( 3 ) =  differenceBetweenResultsForTest4.norm( );
    differenceBetweenResultsForTest5 = computedResultForTest5 - expectedResultForTest5;
    differenceBetweenResults( 4 ) =  differenceBetweenResultsForTest5.norm( );

    // Set test result to false if the test does not match the expected result.
    if ( differenceBetweenResults.norm( ) > std::numeric_limits< double >::epsilon( ) )
    {
        // Check if the source of the error is Test 1.
        if ( differenceBetweenResults( 0 ) > std::numeric_limits< double >::epsilon( ) )
        {
            cerr << "The computed value ( " << computedResultForTest1
                 << " ) for a user-defined gravitational parameter does not "
                 << "match the set value ( " << expectedResultForTest1
                 << " )." << endl;
            cerr << "The difference is: "
                 << differenceBetweenResults( 0 ) << endl;
        }

        // Check if the source of the error is Test 2.
        if ( differenceBetweenResults( 1 ) > std::numeric_limits< double >::epsilon( ) )
        {
            cerr << "The retrieved value ( " << computedResultForTest2
                 << " ) for a predefined gravitational parameter does not "
                 << "match the set value ( " << expectedResultForTest2
                 << " )." << endl;
            cerr << "The difference is: "
                 << differenceBetweenResults( 1 ) << endl;
        }

        // Check if the source of the error is Test 3.
        if ( differenceBetweenResults( 2 ) > std::numeric_limits< double >::epsilon( ) )
        {
            cerr << "The computed value ( " << computedResultForTest3
                 << " ) for the gravitational potential does not match "
                 << "the expected solution ( " << expectedResultForTest3
                 << " )." << endl;
            cerr << "The difference is: "
                 << differenceBetweenResults( 2 ) << endl;
        }

        // Check if the source of the error is Test 4.
        if ( differenceBetweenResults( 3 ) > std::numeric_limits< double >::epsilon( ) )
        {
            cerr << "The computed value ( " << computedResultForTest3
                 << " ) for the gravitational gradient does not match "
                 << "the expected solution  (" << expectedResultForTest3
                 << " )." << endl;
            cerr << "The difference is: "
                 << differenceBetweenResults( 4 ) << endl;
        }

        // Check if the source of the error is Test 5.
        if ( differenceBetweenResults( 4 ) > std::numeric_limits< double >::epsilon( ) )
        {
            cerr << "The computed value ( " << computedResultForTest5
                 <<  ") for the gravitational gradient tensor does not match "
                 << "the expected solution ( " << expectedResultForTest5 << " )." << endl;
            cerr << "The difference is: " << differenceBetweenResults( 4 ) << endl;
        }

        isSphericalHarmonicsGravityFieldErroneous = true;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isSphericalHarmonicsGravityFieldErroneous )
    {
        cerr << "testSphericalHarmonicsGravityField failed!" << std::endl;
    }

    return isSphericalHarmonicsGravityFieldErroneous;
}
