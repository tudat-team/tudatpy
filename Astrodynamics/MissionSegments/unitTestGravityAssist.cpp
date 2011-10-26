/*! \file unitTestGravityAssist.cpp
 *    Source file of unit test file of gravity assist code. This unit test file will test the
 *    gravity assist code.
 *
 *    Path              : /Astrodynamics/MissionSegments/
 *    Version           : 7
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 17 January, 2011
 *    Last modified     : 27 June, 2011
 *
 *    References
 *
 *    Notes
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
 *      A test should be added together with results from the Lambert targeter
 *      and the Ephemeris class.
 *      The expected result for the current test should be calculated inside
 *      this code, not outside by a calculator.
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110117    E. Iorfida        First creation of the code.
 *      110202    J. Melman         Some small changes to layout and
 *                                  added some notes and comments.
 *      110203    E. Iorfida        Changed the code respect to the
 *                                  modification inside the main file.
 *      110208    E. Iorfida        Update file with CartesianVelocityElements
 *                                  pointers as input.
 *      110212    J. Melman         Added comments to clarify test case.
 *      110512    K. Kumar          Updated code to not use dynamic memory
 *                                  allocation and new
 *                                  createPredefinedCelestialBody() function.
 *      110627    K. Kumar          Updated to use new predefined planets code.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include "Astrodynamics/Bodies/planet.h"
#include "Astrodynamics/MissionSegments/gravityAssist.h"
#include "Mathematics/unitConversions.h"

//! Test of gravity assist code.
int main( )
{
    // Using declarations.
    using std::endl;
    using std::cerr;
    using std::fabs;

    // Test result initialised to false.
    bool isGravityAssistErroneous = false;

    // Tolerances.
    // Expected velocity output is defined with an accuracy of 1 m/s.
    double velocityTolerance = 1.0;

    // In the first test case, the incoming and outgoing inertial
    // velocities are defined such that the hyperbolic excess velocities
    // are equal. In that way, a delta-V is only needed to rotate the
    // velocity vector, which has been calculated by hand.
    // Expected delta-V for a powered swing-by around Mars.
    double expectedDeltaV = 3.652e3;

    // Define body that is swung by.
    Planet predefinedMars;
    predefinedMars.setPredefinedPlanetSettings( Planet::mars );

    // Define Sun gravitational parameter.
    double gravitationalParameterSun = 1.32712440018e20;

    // Define planet-Sun distance.
    double distanceMarsToSun = unit_conversions::convertAstronomicalUnitsToMeters( 1.5 );

    // Define smallest periapsis distance factor.
    double marsSmallestPeriapsisDistanceFactor = 1.076;

    // Declare GravityAssist object.
    GravityAssist myGravityAssist;

    // Define planet heliocentric velocity vector.
    // The orbit is considered to be circular.
    Vector3d marsVelocity;
    marsVelocity.x( ) = 0.0;
    marsVelocity.y( ) = sqrt( gravitationalParameterSun / distanceMarsToSun );
    marsVelocity.z( ) = 0.0;

    // Define pointer to satellite incoming vector.
    CartesianVelocityElements incomingVelocityTest;
    incomingVelocityTest.setCartesianElementXDot( -25.0e3 * sin( M_PI / 6.0 ) );
    incomingVelocityTest.setCartesianElementYDot( 25.0e3 * cos( M_PI / 6.0 ) );
    incomingVelocityTest.setCartesianElementZDot( 0.0 );

    // Define pointer to satellite outgoing vector.
    CartesianVelocityElements outgoingVelocityTest;
    outgoingVelocityTest.setCartesianElementXDot(
                incomingVelocityTest.getCartesianElementXDot( ) );
    outgoingVelocityTest.setCartesianElementYDot(
                2.0 * marsVelocity.y( ) - incomingVelocityTest.getCartesianElementYDot( ) );
    outgoingVelocityTest.setCartesianElementZDot( 0.0 );

    // Set values to compute gravity assist code.
    SphereSegment sphereSegment;
    sphereSegment.setRadius( 3398.0e3 );
    predefinedMars.setShapeModel( &sphereSegment );
    myGravityAssist.setCentralBody( &predefinedMars );
    myGravityAssist.setCentralBodyVelocity( marsVelocity );
    myGravityAssist.setSmallestPeriapsisDistanceFactor( marsSmallestPeriapsisDistanceFactor );
    myGravityAssist.setPointerToIncomingVelocity( &incomingVelocityTest );
    myGravityAssist.setPointerToOutgoingVelocity( &outgoingVelocityTest );

    // Create pointers to new Newton Raphson objects.
    NewtonRaphson myNewtonRaphsonGravityAssist;
    myGravityAssist.setNewtonRaphsonMethod( &myNewtonRaphsonGravityAssist );

    // Compute powered gravity-assist implementation.
    double deltaV = myGravityAssist.computeDeltaV( );

    // Set test result to true if the test does not match the expected results.
    if ( fabs( deltaV - expectedDeltaV ) >= velocityTolerance )
    {
        // Set error flag to true.
        isGravityAssistErroneous = true;

        // Generate error statements.
        cerr << "The computed value of delta-V for the "
                "case of equal hyperbolic excess velocities ( " << deltaV
             << " ) using the powered gravity-assist algorithm "
             << "does not match the expected solution ( " << expectedDeltaV << " )." << endl;
        cerr << "The relative error is: "  << fabs( deltaV - expectedDeltaV ) << endl;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    if ( isGravityAssistErroneous )
    {
        cerr << "testGravityAssist failed!" << endl;
    }

    return isGravityAssistErroneous;
}

// End of file.
