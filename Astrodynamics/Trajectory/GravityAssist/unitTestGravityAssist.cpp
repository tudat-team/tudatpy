/*! \file unitTestGravityAssist.cpp
 *    Source file of unit test file of gravity assist code.
 *    This unit test file will test the gravity assist code.
 *
 *    Path              : /Astrodynamics/Trajectory/GravityAssist/
 *    Version           : 5
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
 *    Last modified     : 12 February, 2011
 *
 *    References
 *
 *    Notes
 *      A test should be added together with results from the Lambert targeter
 *      and the Ephemeris class.
 *      The expected result for the current test should be calculated inside
 *      this code, not outside by a calculator.
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
 *      110117    E. Iorfida          First creation of the code.
 *      110202    J. Melman           Some small changes to layout and
 *                                    added some notes and comments.
 *      110203    E. Iorfida          Changed the code respect to the
 *                                    modification inside the main file.
 *      110208    E. Iorfida          Update file with CartesianVelocityElements
 *                                    pointers as input.
 *      110212    J. Melman           Added comments to clarify test case.
 */

// Include statements.
#include "unitTestGravityAssist.h"
#include "predefinedPlanets.h"

// Using directives.
using mathematics::computeAbsoluteValue;
using predefined_planets::createPredefinedPlanet;
using std::endl;
using std::cerr;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test of gravity assist code.
bool testGravityAssist( )
{
    // Test result initialised to false.
    bool isGravityAssistErroneous = false;

    // Tolerances.
    // Velocity outputs have a tolerance of 1 m/s.
    double eccentricityTolerance = 1.0e-13;
    double velocityTolerance = 1.0;

    // In the first test case, the incoming and outgoing inertial
    // velocities are defined such that the hyperbolic excess velocities
    // are equal. In that way, a delta-V is only needed to rotate the
    // velocity vector, which has been calculated by hand.
    // Expected delta-V for a powered swing-by around Mars.
    double expectedDeltaV = 3.652e3;

    // Define body that is swung by.
    CelestialBody* pointerToMars = new CelestialBody;
    pointerToMars = createPredefinedPlanet( predefined_planets::mars );

    // Define Sun gravitational parameter.
    double gravitationalParameterSun = 1.32712440018e20;

    // Define planet-Sun distance.
    double distanceMarsToSun =
            unit_conversions::convertAstronomicalUnitsToMeters( 1.5 );

    // Define planet radius.
    double marsRadius = 3398.0e3;

    // Define smallest periapsis distance factor.
    double marsSmallestPeriapsisDistanceFactor = 1.076;

    // Declare GravityAssist object.
    GravityAssist myGravityAssist_;

    // Define planet heliocentric velocity vector.
    // The orbit is considered to be circular.
    Vector3d marsVelocity_;
    marsVelocity_.x( ) = 0.0;
    marsVelocity_.y( ) =
            sqrt( gravitationalParameterSun / distanceMarsToSun );
    marsVelocity_.z( ) = 0.0;

    // Define pointer to satellite incoming vector.
    CartesianVelocityElements* pointerToIncomingVelocityTest =
            new CartesianVelocityElements;
    pointerToIncomingVelocityTest->
            setCartesianElementXDot( -25.0e3 * sin( M_PI / 6.0 ) );
    pointerToIncomingVelocityTest->
            setCartesianElementYDot( 25.0e3 * cos( M_PI / 6.0 ) );
    pointerToIncomingVelocityTest->
            setCartesianElementZDot( 0.0 );

    // Define pointer to satellite outgoing vector.
    CartesianVelocityElements* pointerToOutgoingVelocityTest =
            new CartesianVelocityElements;
    pointerToOutgoingVelocityTest->setCartesianElementXDot(
            pointerToIncomingVelocityTest->getCartesianElementXDot( ) );
    pointerToOutgoingVelocityTest->setCartesianElementYDot(
            2.0 * marsVelocity_.y( ) -
            pointerToIncomingVelocityTest->getCartesianElementYDot( ) );
    pointerToOutgoingVelocityTest->setCartesianElementZDot( 0.0 );

    // Set values to compute gravity assist code.
    myGravityAssist_.setCentralBody( pointerToMars );
    myGravityAssist_.setCentralBodyRadius( marsRadius );
    myGravityAssist_.setCentralBodyVelocity( marsVelocity_ );
    myGravityAssist_.setSmallestPeriapsisDistanceFactor(
            marsSmallestPeriapsisDistanceFactor );
    myGravityAssist_.setPointerToIncomingVelocity(
            pointerToIncomingVelocityTest );
    myGravityAssist_.setPointerToOutgoingVelocity(
            pointerToOutgoingVelocityTest );

    // Create pointers to new Newton Raphson objects.
    NewtonRaphson* pointerToMyNewtonRaphsonGravityAssist =
            new NewtonRaphson;
    myGravityAssist_.setNewtonRaphsonMethod(
            pointerToMyNewtonRaphsonGravityAssist );

    // Set boolean to true, to apply Newton-Raphson method anyway,
    // although the excess velocities are equal, and thus normally
    // no iterator would be required.
    myGravityAssist_.isRootFinderRequiredForChecking = true;

    // Compute powered gravity-assist implementation.
    double deltaV = myGravityAssist_.computeDeltaV( );

    // Set test result to true if the test does not match the expected results.
    if ( computeAbsoluteValue( deltaV - expectedDeltaV ) >= velocityTolerance ||
         computeAbsoluteValue( myGravityAssist_.incomingEccentricity -
            myGravityAssist_.outgoingEccentricity ) >= eccentricityTolerance )
    {
        // Set error flag to true.
        isGravityAssistErroneous = true;

        if ( computeAbsoluteValue( deltaV - expectedDeltaV ) >=
             velocityTolerance )
        {
            // Generate error statements.
            cerr << "The computed value of delta-V for the "
                    "case of equal hyperbolic excess velocities ( "
                 << deltaV
                 << " ) using the powered gravity-assist algorithm "
                 << "does not match the expected solution ( "
                 << expectedDeltaV << " )." << endl;
            cerr << "The relative error is: "
                 << computeAbsoluteValue( deltaV - expectedDeltaV ) << endl;
        }

        if ( computeAbsoluteValue(
                myGravityAssist_.incomingEccentricity -
                myGravityAssist_.outgoingEccentricity )
            >= eccentricityTolerance )
        {
            // Generate error statements.
            cerr << "The value of eccentricity of the incoming hyperbolic leg ( "
                 << myGravityAssist_.incomingEccentricity
                 << " ) using the powered gravity-assist algorithm "
                 << "does not match the value of eccentricity of the"
                    " outgoing hyperbolic leg ( "
                 << myGravityAssist_.outgoingEccentricity << " )."
                 << endl;
            cerr << "The relative error is: "
                 << computeAbsoluteValue(
                    myGravityAssist_.incomingEccentricity -
                    myGravityAssist_.outgoingEccentricity ) << endl;
        }
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isGravityAssistErroneous;
}

}

// End of file.

int main( ) {
return unit_tests::testGravityAssist( );
}

