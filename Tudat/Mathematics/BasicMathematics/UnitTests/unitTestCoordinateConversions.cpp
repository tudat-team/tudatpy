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
 *      110207    B. Romgens        File created.
 *      110215    K. Kumar          Minor modifications to layout, comments
 *                                  and variable-naming.
 *      110411    K. Kumar          Added unit test for
 *                                  convertCartesianToSpherical() function.
 *      110701    K. Kumar          Updated failing tests with relative errors.
 *      110708    K. Kumar          Added unit tests for computeSampleMean()
 *                                  and computeSampleVariance() functions.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      111111    K. Kumar          Strange error with convertCylindricalToCartesian function;
 *                                  achieved precision of results is less than machine precision,
 *                                  fixed by using slightly larger precision tolerance.
 *      120127    D. Dirkx          Moved unit to separate file from basic mathematics test; moved
 *                                  tests to separate functions; moved to Tudat Core.
 *      120127    K. Kumar          Transferred unit tests over to Boost unit test framework.
 *      120128    K. Kumar          Changed BOOST_CHECK to BOOST_CHECK_CLOSE_FRACTION and
 *                                  BOOST_CHECK_SMALL for unit test comparisons.
 *      120118    D. Gondelach      Added unit tests for convertCylindricalToCartesian.
 *                                  Removed unit test for old convertCylindricalToCartesian
 *                                  function.
 *      120214    K. Kumar          Branched from old Tudat trunk for new coordinate conversions.
 *
 *    References
 *
 */

#include <cmath>
#include <iostream>
#include <Eigen/Core>
#include <limits>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

int main( )
{
    using std::cerr;
    using std::endl;
    using tudat::mathematics::PI;

    double epsilon = std::numeric_limits< double >::epsilon( );

    // Declare and initialize test result to false.
    bool isBasicMathematicsFunctionsErroneous = false;

    // Test conversion from cylindrical (r, theta, z) to Cartesian (x, y, z) coordinates.
    // Test 33: Test conversion of ( 0.0, pi, 1.2 ).
    // Test 34: Test conversion of ( 2.3, -pi/2, -3.5 ).

    // Cylindrical coordinates (r,theta,z).
    Eigen::Vector3d cylindricalCoordinatesTest33, cylindricalCoordinatesTest34;
    cylindricalCoordinatesTest33 << 0.0, PI, 1.2;
    cylindricalCoordinatesTest34 << 2.3, -PI / 2.0, -3.5;

    // Expected Cartesian coordinates (x, y, z).
    Eigen::Vector3d expectedCartesianCoordinatesTest33, expectedCartesianCoordinatesTest34;
    expectedCartesianCoordinatesTest33 << 0.0, 0.0, 1.2;
    expectedCartesianCoordinatesTest34 << 2.3 * cos( -PI / 2.0 ), 2.3*sin( -PI / 2.0 ), -3.5;

    // Test Cartesian coordinates (x, y, z).
    Eigen::Vector3d cartesianCoordinatesTest33 = tudat::mathematics::coordinate_conversions::
            convertCylindricalToCartesian( cylindricalCoordinatesTest33 );
    Eigen::Vector3d cartesianCoordinatesTest34 = tudat::mathematics::coordinate_conversions::
            convertCylindricalToCartesian( cylindricalCoordinatesTest34 );

    if ( fabs( cartesianCoordinatesTest33(0) - expectedCartesianCoordinatesTest33(0) )
         > epsilon ||
         fabs( cartesianCoordinatesTest33(1) - expectedCartesianCoordinatesTest33(1) )
         > epsilon ||
         fabs( cartesianCoordinatesTest33(2) - expectedCartesianCoordinatesTest33(2) )
         / cartesianCoordinatesTest33(2) > epsilon )
    {
        cerr << "The convertCylindricalToCartesian function does not "
             << "function correctly, as the computed coordinates: ( "
             << cartesianCoordinatesTest33(0) << ", " << cartesianCoordinatesTest33(1)
             << ", " << cartesianCoordinatesTest33(2)
             << " ) do not match the expected coordinates: ( "
             << expectedCartesianCoordinatesTest33(0) << ", "
             << expectedCartesianCoordinatesTest33(1) << ", "
             << expectedCartesianCoordinatesTest33(2) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    if ( fabs( cartesianCoordinatesTest34(0) - expectedCartesianCoordinatesTest34(0) )
          > epsilon ||
         fabs( cartesianCoordinatesTest34(1) - expectedCartesianCoordinatesTest34(1) )
          / expectedCartesianCoordinatesTest34(1) > epsilon ||
         fabs( cartesianCoordinatesTest34(2) - expectedCartesianCoordinatesTest34(2) )
          / expectedCartesianCoordinatesTest34(2) > epsilon )
    {
        cerr << "The convertCylindricalToCartesian function does not "
             << "function correctly, as the computed coordinates: ( "
             << cartesianCoordinatesTest34(0) << ", " << cartesianCoordinatesTest34(1)
             << ", " << cartesianCoordinatesTest34(2)
             << " ) do not match the expected coordinates: ( "
             << expectedCartesianCoordinatesTest34(0) << ", "
             << expectedCartesianCoordinatesTest34(1) << ", "
             << expectedCartesianCoordinatesTest34(2) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    // Test conversion from cylindrical (r, theta, z, Vr, Vtheta, Vz) to
    // Cartesian (x, y, z, xdot, ydot, zdot) state.
    // Test 35: Test conversion of (2.1, pi/2.0, 1.2, 5.4, 4.5, -3.9)
    // Test 36: Test conversion of (0.0, 8.2*pi/3.0, -2.5, -5.8, 0.0, 1.7)

    // Cylindrical states (r, theta, z, Vr, Vtheta, Vz).
    Eigen::VectorXd cylindricalStateTest35( 6 ), cylindricalStateTest36( 6 );
    cylindricalStateTest35 << 2.1, PI / 2.0, 1.2, 5.4, 4.5, -3.9;
    cylindricalStateTest36 << 0.0, 8.2 * PI / 3.0, -2.5, -5.8, 0.0, 1.7;

    // Expected Cartesian states (x, y, z, xdot, ydot, zdot).
    Eigen::VectorXd expectedCartesianStateTest35( 6 ), expectedCartesianStateTest36( 6 );
    expectedCartesianStateTest35 << 2.1 * cos( PI / 2.0 ),
                                    2.1 * sin( PI / 2.0 ),
                                    1.2,
                                    5.4 * cos( PI / 2.0 ) - 4.5 * sin( PI / 2.0 ),
                                    5.4 * sin( PI / 2.0 ) + 4.5 * cos( PI / 2.0 ),
                                    -3.9;
    expectedCartesianStateTest36 << 0.0,
                                    0.0,
                                    -2.5,
                                    -5.8 * cos( 8.2 * PI / 3.0 ),
                                    -5.8 * sin( 8.2 * PI / 3.0 ),
                                    1.7;

    // Test Cartesian states (x, y, z, xdot, ydot, zdot).
    Eigen::VectorXd cartesianStateTest35 = tudat::mathematics::coordinate_conversions::
            convertCylindricalToCartesian( cylindricalStateTest35 );
    Eigen::VectorXd cartesianStateTest36 = tudat::mathematics::coordinate_conversions::
            convertCylindricalToCartesian( cylindricalStateTest36 );

    if ( fabs( cartesianStateTest35(0) - expectedCartesianStateTest35(0) )
          > epsilon ||
         fabs( cartesianStateTest35(1) - expectedCartesianStateTest35(1) )
          / expectedCartesianStateTest35(1) > epsilon ||
         fabs( cartesianStateTest35(2) - expectedCartesianStateTest35(2) )
          / expectedCartesianStateTest35(2) > epsilon ||
         fabs( cartesianStateTest35(3) - expectedCartesianStateTest35(3) )
          / expectedCartesianStateTest35(3) > epsilon ||
         fabs( cartesianStateTest35(4) - expectedCartesianStateTest35(4) )
          / expectedCartesianStateTest35(4) > epsilon ||
         fabs( cartesianStateTest35(5) - expectedCartesianStateTest35(5) )
          / expectedCartesianStateTest35(5) > epsilon )
    {
        cerr << "The convertCylindricalToCartesianState function does not "
             << "function correctly, as the computed sate: ( "
             << cartesianStateTest35(0) << ", " << cartesianStateTest35(1) << ", "
             << cartesianStateTest35(2) << ", " << cartesianStateTest35(3) << ", "
             << cartesianStateTest35(4) << ", " << cartesianStateTest35(5)
             << " ) do not match the expected coordinates: ( "
             << expectedCartesianStateTest35(0) << ", " << expectedCartesianStateTest35(1) << ", "
             << expectedCartesianStateTest35(2) << ", " << expectedCartesianStateTest35(3) << ", "
             << expectedCartesianStateTest35(4) << ", " << expectedCartesianStateTest35(5)
             << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    if ( fabs( cartesianStateTest36(0) - expectedCartesianStateTest36(0) )
         > epsilon ||
         fabs( cartesianStateTest36(1) - expectedCartesianStateTest36(1) )
         > epsilon ||
         fabs( cartesianStateTest36(2) - expectedCartesianStateTest36(2) )
          / expectedCartesianStateTest36(2) > epsilon ||
         fabs( cartesianStateTest36(3) - expectedCartesianStateTest36(3) )
          / expectedCartesianStateTest36(3) > epsilon ||
         fabs( cartesianStateTest36(4) - expectedCartesianStateTest36(4) )
          / expectedCartesianStateTest36(4) > epsilon ||
         fabs( cartesianStateTest36(5) - expectedCartesianStateTest36(5) )
          / expectedCartesianStateTest36(5) > epsilon )
    {
        cerr << "The convertCylindricalToCartesianState function does not "
             << "function correctly, as the computed sate: ( "
             << cartesianStateTest36(0) << ", " << cartesianStateTest36(1) << ", "
             << cartesianStateTest36(2) << ", " << cartesianStateTest36(3) << ", "
             << cartesianStateTest36(4) << ", " << cartesianStateTest36(5)
             << " ) do not match the expected coordinates: ( "
             << expectedCartesianStateTest36(0) << ", " << expectedCartesianStateTest36(1) << ", "
             << expectedCartesianStateTest36(2) << ", " << expectedCartesianStateTest36(3) << ", "
             << expectedCartesianStateTest36(4) << ", " << expectedCartesianStateTest36(5)
             << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    // Test conversion from Cartesian (x, y, z) to cylindrical (r, theta, z) coordinates.
    // Test 37: Test conversion of ( 0.0, 0.0, 1.0 ).
    // Test 38: Test conversion of ( 0.0, 2.0, 1.0 ).
    // Test 39: Test conversion of ( 0.0, -2.0, -1.0 ).
    // Test 40: Test conversion of ( -5.0, -8.0, 5.0 ).

    // Cartesian coordinates (x, y, z).
    Eigen::Vector3d cartesianCoordinatesTest37, cartesianCoordinatesTest38,
                    cartesianCoordinatesTest39, cartesianCoordinatesTest40;
    cartesianCoordinatesTest37 << 0.0, 0.0, 0.0;
    cartesianCoordinatesTest38 << 0.0, 2.0, 1.0;
    cartesianCoordinatesTest39 << 0.0, -2.0, -1.0;
    cartesianCoordinatesTest40 << -5.0, -8.0, 5.0;

    // Expected cylindrical coordinates (r, theta, z).
    Eigen::Vector3d expectedCylindricalCoordinatesTest37, expectedCylindricalCoordinatesTest38,
                    expectedCylindricalCoordinatesTest39, expectedCylindricalCoordinatesTest40;
    expectedCylindricalCoordinatesTest37 << 0.0, 0.0, 0.0;
    expectedCylindricalCoordinatesTest38 << 2.0, PI / 2.0, 1.0;
    expectedCylindricalCoordinatesTest39 << 2.0, 3.0 * PI / 2.0, -1.0;
    expectedCylindricalCoordinatesTest40 << sqrt( 25.0 + 64.0 ), atan2(-8.0,-5.0) + 2.0 * PI, 5.0;

    // Test cylindrical coordinates (r, theta, z).
    Eigen::Vector3d cylindricalCoordinatesTest37 = tudat::mathematics::coordinate_conversions::
            convertCartesianToCylindrical( cartesianCoordinatesTest37 );
    Eigen::Vector3d cylindricalCoordinatesTest38 = tudat::mathematics::coordinate_conversions::
            convertCartesianToCylindrical( cartesianCoordinatesTest38 );
    Eigen::Vector3d cylindricalCoordinatesTest39 = tudat::mathematics::coordinate_conversions::
            convertCartesianToCylindrical( cartesianCoordinatesTest39 );
    Eigen::Vector3d cylindricalCoordinatesTest40 = tudat::mathematics::coordinate_conversions::
            convertCartesianToCylindrical( cartesianCoordinatesTest40 );

    if ( fabs( cylindricalCoordinatesTest37(0) - expectedCylindricalCoordinatesTest37(0) )
         > epsilon ||
         fabs( cylindricalCoordinatesTest37(1) - expectedCylindricalCoordinatesTest37(1) )
         > epsilon ||
         fabs( cylindricalCoordinatesTest37(2) - expectedCylindricalCoordinatesTest37(2) )
         > epsilon )
    {
        cerr << "The tudat::mathematics::coordinate_conversions::convertCartesianToCylindrical function does not "
             << "function correctly, as the computed coordinates: ( "
             << cylindricalCoordinatesTest37(0) << ", " << cylindricalCoordinatesTest37(1)
             << ", " << cylindricalCoordinatesTest37(2)
             << " ) do not match the expected coordinates: ( "
             << expectedCylindricalCoordinatesTest37(0) << ", "
             << expectedCylindricalCoordinatesTest37(1) << ", "
             << expectedCylindricalCoordinatesTest37(2) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    if ( fabs( cylindricalCoordinatesTest38(0) - expectedCylindricalCoordinatesTest38(0) )
          / expectedCylindricalCoordinatesTest38(0) > epsilon ||
         fabs( cylindricalCoordinatesTest38(1) - expectedCylindricalCoordinatesTest38(1) )
          / expectedCylindricalCoordinatesTest38(1) > epsilon ||
         fabs( cylindricalCoordinatesTest38(2) - expectedCylindricalCoordinatesTest38(2) )
          / expectedCylindricalCoordinatesTest38(2) > epsilon )
    {
        cerr << "The tudat::mathematics::coordinate_conversions::convertCartesianToCylindrical function does not "
             << "function correctly, as the computed coordinates: ( "
             << cylindricalCoordinatesTest38(0) << ", " << cylindricalCoordinatesTest38(1)
             << ", " << cylindricalCoordinatesTest38(2)
             << " ) do not match the expected coordinates: ( "
             << expectedCylindricalCoordinatesTest38(0) << ", "
             << expectedCylindricalCoordinatesTest38(1) << ", "
             << expectedCylindricalCoordinatesTest38(2) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    if ( fabs( cylindricalCoordinatesTest39(0) - expectedCylindricalCoordinatesTest39(0) )
          / expectedCylindricalCoordinatesTest39(0) > epsilon ||
         fabs( cylindricalCoordinatesTest39(1) - expectedCylindricalCoordinatesTest39(1) )
          / expectedCylindricalCoordinatesTest39(1) > epsilon ||
         fabs( cylindricalCoordinatesTest39(2) - expectedCylindricalCoordinatesTest39(2) )
          / expectedCylindricalCoordinatesTest39(2) > epsilon )
    {
        cerr << "The tudat::mathematics::coordinate_conversions::convertCartesianToCylindrical function does not "
             << "function correctly, as the computed coordinates: ( "
             << cylindricalCoordinatesTest39(0) << ", " << cylindricalCoordinatesTest39(1)
             << ", " << cylindricalCoordinatesTest39(2)
             << " ) do not match the expected coordinates: ( "
             << expectedCylindricalCoordinatesTest39(0) << ", "
             << expectedCylindricalCoordinatesTest39(1) << ", "
             << expectedCylindricalCoordinatesTest39(2) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    if ( fabs( cylindricalCoordinatesTest40(0) - expectedCylindricalCoordinatesTest40(0) )
          / expectedCylindricalCoordinatesTest40(0) > epsilon ||
         fabs( cylindricalCoordinatesTest40(1) - expectedCylindricalCoordinatesTest40(1) )
          / expectedCylindricalCoordinatesTest40(1) > epsilon ||
         fabs( cylindricalCoordinatesTest40(2) - expectedCylindricalCoordinatesTest40(2) )
          / expectedCylindricalCoordinatesTest40(2) > epsilon )
    {
        cerr << "The tudat::mathematics::coordinate_conversions::convertCartesianToCylindrical function does not "
             << "function correctly, as the computed coordinates: ( "
             << cylindricalCoordinatesTest40(0) << ", " << cylindricalCoordinatesTest40(1)
             << ", " << cylindricalCoordinatesTest40(2)
             << " ) do not match the expected coordinates: ( "
             << expectedCylindricalCoordinatesTest40(0) << ", "
             << expectedCylindricalCoordinatesTest40(1) << ", "
             << expectedCylindricalCoordinatesTest40(2) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    // Test conversion from Cartesian (x, y, z, xdot, ydot, zdot) to
    // cylindrical (r, theta, z, Vr, Vtheta, Vz) state.
    // Test 41: Test conversion of ( 0.0, 0.0, 1.0, 5.0, 6.0, -9.0 ).
    // Test 42: Test conversion of ( 2.0, 0.0, -5.0, -4.0, 6.0, -6.0 ).
    // Test 43: Test conversion of ( -7.0, -4.0, 3.0, 5.0, -3.0, 7.0 ).

    // Cartesian states (x,y,z,xdot,ydot,zdot).
    Eigen::VectorXd cartesianStateTest41( 6 ), cartesianStateTest42( 6 ),
                    cartesianStateTest43( 6 );
    cartesianStateTest41 << 0.0, 0.0, 1.0, 5.0, 6.0, -9.0;
    cartesianStateTest42 << 2.0, 0.0, -5.0, -4.0, 6.0, -6.0;
    cartesianStateTest43 << -7.0, -4.0, 3.0, 5.0, -3.0, 7.0;

    // Expected cylindrical states (r, theta, z, Vr, Vtheta, Vz).
    Eigen::VectorXd expectedCylindricalStateTest41( 6 ), expectedCylindricalStateTest42( 6 ),
                    expectedCylindricalStateTest43( 6 );
    expectedCylindricalStateTest41 << 0.0, 0.0, 1.0, sqrt( 25.0 + 36.0 ), 0.0, -9.0;
    expectedCylindricalStateTest42 << 2.0, 0.0, -5.0, -4.0, 6.0, -6.0;
    expectedCylindricalStateTest43 << sqrt(49.0+16.0), atan2( -4.0, -7.0 ) + 2.0 * PI, 3.0,
            ( -7.0 * 5.0 + ( -4.0 ) * -3.0 ) / sqrt( 49.0 + 16.0),
            ( -7.0 * -3.0 - ( -4.0 ) * 5.0 ) / sqrt( 49.0 + 16.0), 7.0;

    // Test cylindrical states (r, theta, z, Vr, Vtheta, Vz).
    Eigen::VectorXd cylindricalStateTest41 = tudat::mathematics::coordinate_conversions::
            convertCartesianToCylindrical( cartesianStateTest41 );
    Eigen::VectorXd cylindricalStateTest42 = tudat::mathematics::coordinate_conversions::
            convertCartesianToCylindrical( cartesianStateTest42 );
    Eigen::VectorXd cylindricalStateTest43 = tudat::mathematics::coordinate_conversions::
            convertCartesianToCylindrical( cartesianStateTest43 );

    if ( fabs( cylindricalStateTest41(0) - expectedCylindricalStateTest41(0) )
         > epsilon ||
         fabs( cylindricalStateTest41(1) - expectedCylindricalStateTest41(1) )
         > epsilon ||
         fabs( cylindricalStateTest41(2) - expectedCylindricalStateTest41(2) )
          / expectedCylindricalStateTest41(2) > epsilon ||
         fabs( cylindricalStateTest41(3) - expectedCylindricalStateTest41(3) )
          / expectedCylindricalStateTest41(3) > epsilon ||
         fabs( cylindricalStateTest41(4) - expectedCylindricalStateTest41(4) )
         > epsilon ||
         fabs( cylindricalStateTest41(5) - expectedCylindricalStateTest41(5) )
          / expectedCylindricalStateTest41(5) > epsilon )
    {
        cerr << "The tudat::mathematics::coordinate_conversions::convertCartesianToCylindrical function does not "
             << "function correctly, as the computed sate: ( "
             << cylindricalStateTest41(0) << ", " << cylindricalStateTest41(1) << ", "
             << cylindricalStateTest41(2) << ", " << cylindricalStateTest41(3) << ", "
             << cylindricalStateTest41(4) << ", " << cylindricalStateTest41(5)
             << " ) do not match the expected coordinates: ( "
             << expectedCylindricalStateTest41(0) << ", "
             << expectedCylindricalStateTest41(1) << ", "
             << expectedCylindricalStateTest41(2) << ", "
             << expectedCylindricalStateTest41(3) << ", "
             << expectedCylindricalStateTest41(4) << ", "
             << expectedCylindricalStateTest41(5) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    if ( fabs( cylindricalStateTest42(0) - expectedCylindricalStateTest42(0) )
          / expectedCylindricalStateTest42(0) > epsilon ||
         fabs( cylindricalStateTest42(1) - expectedCylindricalStateTest42(1) )
          > epsilon ||
         fabs( cylindricalStateTest42(2) - expectedCylindricalStateTest42(2) )
          / expectedCylindricalStateTest42(2) > epsilon ||
         fabs( cylindricalStateTest42(3) - expectedCylindricalStateTest42(3) )
          / expectedCylindricalStateTest42(3) > epsilon ||
         fabs( cylindricalStateTest42(4) - expectedCylindricalStateTest42(4) )
          / expectedCylindricalStateTest42(4) > epsilon ||
         fabs( cylindricalStateTest42(5) - expectedCylindricalStateTest42(5) )
          / expectedCylindricalStateTest42(5) > epsilon )
    {
        cerr << "The tudat::mathematics::coordinate_conversions::convertCartesianToCylindrical function does not "
             << "function correctly, as the computed sate: ( "
             << cylindricalStateTest42(0) << ", " << cylindricalStateTest42(1) << ", "
             << cylindricalStateTest42(2) << ", " << cylindricalStateTest42(3) << ", "
             << cylindricalStateTest42(4) << ", " << cylindricalStateTest42(5)
             << " ) do not match the expected coordinates: ( "
             << expectedCylindricalStateTest42(0) << ", "
             << expectedCylindricalStateTest42(1) << ", "
             << expectedCylindricalStateTest42(2) << ", "
             << expectedCylindricalStateTest42(3) << ", "
             << expectedCylindricalStateTest42(4) << ", "
             << expectedCylindricalStateTest42(5) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    if ( fabs( cylindricalStateTest43(0) - expectedCylindricalStateTest43(0) )
          / expectedCylindricalStateTest43(0) > epsilon ||
         fabs( cylindricalStateTest43(1) - expectedCylindricalStateTest43(1) )
          / expectedCylindricalStateTest43(1) > epsilon ||
         fabs( cylindricalStateTest43(2) - expectedCylindricalStateTest43(2) )
          / expectedCylindricalStateTest43(2) > epsilon ||
         fabs( cylindricalStateTest43(3) - expectedCylindricalStateTest43(3) )
          / expectedCylindricalStateTest43(3) > epsilon ||
         fabs( cylindricalStateTest43(4) - expectedCylindricalStateTest43(4) )
          / expectedCylindricalStateTest43(4) > epsilon ||
         fabs( cylindricalStateTest43(5) - expectedCylindricalStateTest43(5) )
          / expectedCylindricalStateTest43(5) > epsilon )
    {
        cerr << "The tudat::mathematics::coordinate_conversions::convertCartesianToCylindrical function does not "
             << "function correctly, as the computed sate: ( "
             << cylindricalStateTest43(0) << ", " << cylindricalStateTest43(1) << ", "
             << cylindricalStateTest43(2) << ", " << cylindricalStateTest43(3) << ", "
             << cylindricalStateTest43(4) << ", " << cylindricalStateTest43(5)
             << " ) do not match the expected coordinates: ( "
             << expectedCylindricalStateTest43(0) << ", "
             << expectedCylindricalStateTest43(1) << ", "
             << expectedCylindricalStateTest43(2) << ", "
             << expectedCylindricalStateTest43(3) << ", "
             << expectedCylindricalStateTest43(4) << ", "
             << expectedCylindricalStateTest43(5) << " )" << endl;

        isBasicMathematicsFunctionsErroneous = true;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isBasicMathematicsFunctionsErroneous )
    {
        cerr << "testBasicMathematicsFunctions failed!" << endl;
    }

    return isBasicMathematicsFunctionsErroneous;
}
