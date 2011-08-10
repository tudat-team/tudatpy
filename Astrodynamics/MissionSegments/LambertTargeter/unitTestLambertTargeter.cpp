/*! \file unitTestLambertTargeter.cpp
 *    Source file of unit test file of Lambert targeting algorithm code.
 *    This unit test file will test the Lambert targeting algorithm code for
 *    both the hyperbolic case and the elliptical case.
 *
 *    Path              : /Astrodynamics/MissionSegments/LambertTargeter/
 *    Version           : 8
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
 *    Date created      : 13 January, 2011
 *    Last modified     : 27 June, 2011
 *
 *    References        :
 *      Mengali, G., and A.A. Quarta, Fondamenti di Meccanica del volo Spaziale.
 *      Noomen, R., Lambert targeter Excel file.
 *
 *    Notes
 *      DISCLAIMER: At the moment, the Lambert targeter only converges for
 *      about half of the cases. This is not evident from the tests below, but
 *      it was observed during simulations carried out by the author. The
 *      reason might very well be an erroneous definition of the starters.
 *
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
 *
 *      The elliptical case was taken from Example 6.1, page 159-162 of
 *      ( Mengali, Quarta ). The hyperbolic case was taken from ( Noomen, R. ).
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
 *      110113    E. Iorfida        First creation of the code.
 *      110126    E. Iorfida        Added test for velocities.
 *      110130    J. Melman         Simplified variable names, e.g.,
 *                                  'normOfVelocityVector' became 'speed'.
 *                                  Also corrected 'tangential' to
 *                                  Suggested less trivial test case.
 *      110201    E. Iorfida        Added pointerToCelestialBody and modified
 *                                  variable names (from heliocentric, to
 *                                  inertial). Added elliptical case and check
 *                                  for the anti-clockwise direction of motion.
 *      110207    E. Iorfida        Added a more stric value for tolerance of
 *                                  semi-major axis for hyperbolic case.
 *      110208    E. Iorfida        Added CartesianPositionElements objects as
 *                                  input and CartesianVelocityElements objects
 *                                  as output.
 *      110512    K. Kumar          Updated code not to use dynamic memory
 *                                  allocation.
 *      110627    K. Kumar          Updated to use new predefined planets code.
 */

// Include statements.
#include "unitTestLambertTargeter.h"

// Using directives.
using mathematics::computeAbsoluteValue;
using unit_conversions::convertAstronomicalUnitsToMeters;
using std::endl;
using std::cerr;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test of Lambert targeting algorithm code.
bool testLambertTargeter( )
{
    // Test result initialised to false.
    bool isLambertTargeterErroneous = false;

    // Expected test result in meters.
    // Hyperbolic test case (results from excel file [1]).
    double expectedValueOfSemiMajorAxisHyperbola = -1270129.3602e3;
    double expectedValueOfRadialSpeedAtDepartureHyperbola = -0.74546e3;
    double expectedValueOfRadialSpeedAtArrivalHyperbola = 0.69321e3;
    double expectedValueOfTransverseSpeedAtDepartureHyperbola = 0.15674e3;
    double expectedValueOfTransverseSpeedAtArrivalHyperbola = 0.10450e3;

    // Elliptical test case (results from example 6.1 page 159-162 [2]).
    // Set canonical units for Earth (see page 29 [2]).
    double distanceUnit = 6.378136e6;
    double timeUnit = 806.78;

    double expectedValueOfSemiMajorAxisEllipse = 5.4214 * distanceUnit;
    double expectedValueOfRadialSpeedAtDepartureEllipse = 2.73580e3;
    double expectedValueOfRadialSpeedAtArrivalEllipse = 2.97503e3;
    double expectedValueOfTransverseSpeedAtDepartureEllipse = 6.59430e3;
    double expectedValueOfTransverseSpeedAtArrivalEllipse = 3.29715e3;

    // Tolerance in absolute units.
    double toleranceSemiMajorAxisHyperbola = 1.0e2;
    double toleranceSemiMajorAxisEllipse = 1.0e4;
    double toleranceVelocity = 1.0e-02;

    // Time conversions.
    double timeOfFlightInDaysHyperbola = 100.0;
    double timeOfFlightHyperbola = unit_conversions::convertJulianDaysToSeconds(
            timeOfFlightInDaysHyperbola );
    double timeOfFlightEllipse = 5.0 * timeUnit;

    // Central bodies parameters.
    Planet predefinedEarth;
    predefinedEarth.setPredefinedPlanetSettings( Planet::earth );
    GravityFieldModel* pointerToEarthGravityField
            = predefinedEarth.getGravityFieldModel( );
    pointerToEarthGravityField->setGravitationalParameter( 398600.4418e9 );

    // Compute Lambert targeting algorithm.

    // Hyperbolic orbit case.
    LambertTargeter lambertTargeterHyperbola;

    // The starting point is twice as far as L1 and L2, which is not really
    // realistic, but it is not about the case, but about the verification.
    CartesianPositionElements positionAtDepartureHyperbola;
    positionAtDepartureHyperbola.setCartesianElementX(
            convertAstronomicalUnitsToMeters( 0.02 ) );
    positionAtDepartureHyperbola.setCartesianElementY( 0.0 );
    positionAtDepartureHyperbola.setCartesianElementZ( 0.0 );

    CartesianPositionElements positionAtArrivalHyperbola;
    positionAtArrivalHyperbola.setCartesianElementX( 0.0 );
    positionAtArrivalHyperbola.setCartesianElementY(
            convertAstronomicalUnitsToMeters( -0.03 ) );
    positionAtArrivalHyperbola.setCartesianElementZ( 0.0 );

    lambertTargeterHyperbola.setPositionAtDeparture(
            &positionAtDepartureHyperbola );
    lambertTargeterHyperbola.setPositionAtArrival(
            &positionAtArrivalHyperbola );
    lambertTargeterHyperbola.setNumberOfRevolutions( 0 );
    lambertTargeterHyperbola.setTimeOfFlight(
            timeOfFlightHyperbola );
    lambertTargeterHyperbola.setCentralBody( &predefinedEarth );

    // Create pointers to new NewtonRaphson object.
    NewtonRaphson newtonRaphsonLambertHyperbola;
    lambertTargeterHyperbola.setNewtonRaphsonMethod(
            &newtonRaphsonLambertHyperbola );

    // Elliptical orbit case.
    LambertTargeter lambertTargeterEllipse;

    CartesianPositionElements positionAtDepartureEllipse;
    positionAtDepartureEllipse.setCartesianElementX( 2.0 * distanceUnit );
    positionAtDepartureEllipse.setCartesianElementY( 0.0 );
    positionAtDepartureEllipse.setCartesianElementZ( 0.0 );

    CartesianPositionElements positionAtArrivalEllipse;
    positionAtArrivalEllipse.setCartesianElementX( 2.0 * distanceUnit );
    positionAtArrivalEllipse.setCartesianElementY( 2.0 * sqrt( 3.0 )
                                                   * distanceUnit );
    positionAtArrivalEllipse.setCartesianElementZ( 0.0 );

    lambertTargeterEllipse.setPositionAtDeparture(
            &positionAtDepartureEllipse );
    lambertTargeterEllipse.setPositionAtArrival(
            &positionAtArrivalEllipse );
    lambertTargeterEllipse.setNumberOfRevolutions( 0 );
    lambertTargeterEllipse.setTimeOfFlight( timeOfFlightEllipse );
    lambertTargeterEllipse.setCentralBody( &predefinedEarth );

    // Create pointers to new NewtonRaphson object.
    NewtonRaphson newtonRaphsonLambertEllipse;
    lambertTargeterEllipse.setNewtonRaphsonMethod(
            &newtonRaphsonLambertEllipse );

    // Compute Lambert targeting algorithms.
    lambertTargeterHyperbola.execute( );
    lambertTargeterEllipse.execute( );

    // Create local vectors for position and velocity.
    Vector3d positionDepartureHyperbola =
            positionAtDepartureHyperbola.state;

    Vector3d velocityDepartureHyperbola =
            lambertTargeterHyperbola.getInertialVelocityAtDeparture( )->state;

    Vector3d positionDepartureEllipse =
            positionAtDepartureEllipse.state;

    Vector3d velocityDepartureEllipse =
            lambertTargeterEllipse.getInertialVelocityAtDeparture( )->state;

    // Set test result to true if the test does not match the expected result.
    if ( ( computeAbsoluteValue(
            lambertTargeterHyperbola.getLambertSemiMajorAxis( ) -
            expectedValueOfSemiMajorAxisHyperbola ) >=
           toleranceSemiMajorAxisHyperbola ) ||
         ( computeAbsoluteValue(
            lambertTargeterHyperbola.getRadialSpeedAtDeparture( ) -
            expectedValueOfRadialSpeedAtDepartureHyperbola ) >=
           toleranceVelocity ) ||
         ( computeAbsoluteValue(
            lambertTargeterHyperbola.getRadialSpeedAtArrival( )-
            expectedValueOfRadialSpeedAtArrivalHyperbola ) >=
           toleranceVelocity ) ||
         ( computeAbsoluteValue(
            lambertTargeterHyperbola.getTransverseSpeedAtDeparture( ) -
            expectedValueOfTransverseSpeedAtDepartureHyperbola ) >=
           toleranceVelocity ) ||
         ( computeAbsoluteValue(
            lambertTargeterHyperbola.getTransverseSpeedAtArrival( ) -
            expectedValueOfTransverseSpeedAtArrivalHyperbola ) >=
           toleranceVelocity ) ||
         ( computeAbsoluteValue(
            lambertTargeterEllipse.getLambertSemiMajorAxis( ) -
            expectedValueOfSemiMajorAxisEllipse ) >=
           toleranceSemiMajorAxisEllipse ) ||
         ( computeAbsoluteValue(
            lambertTargeterEllipse.getRadialSpeedAtDeparture( ) -
            expectedValueOfRadialSpeedAtDepartureEllipse ) >=
           toleranceVelocity ) ||
         ( computeAbsoluteValue(
            lambertTargeterEllipse.getRadialSpeedAtArrival( ) -
            expectedValueOfRadialSpeedAtArrivalEllipse ) >=
           toleranceVelocity ) ||
         ( computeAbsoluteValue(
            lambertTargeterEllipse.getTransverseSpeedAtDeparture( ) -
            expectedValueOfTransverseSpeedAtDepartureEllipse ) >=
           toleranceVelocity ) ||
         ( computeAbsoluteValue(
            lambertTargeterEllipse.getTransverseSpeedAtArrival( ) -
            expectedValueOfTransverseSpeedAtArrivalEllipse ) >=
           toleranceVelocity ) ||

         // Check anti-clockwise direction of the computed orbit.
         ( positionDepartureHyperbola.cross(
                 velocityDepartureHyperbola ).z( ) < 0.0 ) ||
         ( positionDepartureEllipse.cross(
                 velocityDepartureEllipse).z( ) < 0.0 ) )
    {
        // Set error flag to true.
        isLambertTargeterErroneous = true;

        if ( computeAbsoluteValue(
                lambertTargeterHyperbola.getLambertSemiMajorAxis( ) -
                expectedValueOfSemiMajorAxisHyperbola ) >=
                toleranceSemiMajorAxisHyperbola )
        {
            // Generate error statements.
            cerr << "The computed value of the semi-major axis ( "
                 << lambertTargeterHyperbola.
                    getLambertSemiMajorAxis( )
                 << " ) using the Lambert targeting algorithm "
                 << "does not match the expected solution ( "
                 << expectedValueOfSemiMajorAxisHyperbola << " )." << endl;
            cerr << "The error is: "
                 << computeAbsoluteValue( lambertTargeterHyperbola.
                    getLambertSemiMajorAxis( ) -
                    expectedValueOfSemiMajorAxisHyperbola ) << endl;
        }

        if ( computeAbsoluteValue( lambertTargeterHyperbola.
             getRadialSpeedAtDeparture( ) -
             expectedValueOfRadialSpeedAtDepartureHyperbola ) >=
             toleranceVelocity )
        {
            // Generate error statements.
            cerr << " The computed value of the radial speed at departure ( "
                 << lambertTargeterHyperbola.
                    getRadialSpeedAtDeparture( ) << " ) using "
                    " the Lambert targeting algorithm does not match the"
                    " expected solution (" <<
                    expectedValueOfRadialSpeedAtDepartureHyperbola <<
                    " )." << endl;
            cerr << "The error is: "
                 << computeAbsoluteValue(
                    lambertTargeterHyperbola.
                    getRadialSpeedAtDeparture( ) -
                    expectedValueOfRadialSpeedAtDepartureHyperbola ) << endl;
        }

        if ( computeAbsoluteValue( lambertTargeterHyperbola.
             getRadialSpeedAtArrival( ) -
             expectedValueOfRadialSpeedAtArrivalHyperbola ) >=
             toleranceVelocity )
        {
            // Generate error statements.
            cerr << "The computed value of the radial speed at arrival ( "
                 << lambertTargeterHyperbola.
                    getRadialSpeedAtArrival( ) << " ) using "
                    " the Lambert targeting algorithm does not match the"
                    " expected solution (" <<
                    expectedValueOfRadialSpeedAtArrivalHyperbola <<
                    " )." << endl;
            cerr << "The error is: "
                 << computeAbsoluteValue(
                    lambertTargeterHyperbola.
                    getRadialSpeedAtArrival( ) -
                    expectedValueOfRadialSpeedAtArrivalHyperbola ) << endl;
        }

        if ( computeAbsoluteValue( lambertTargeterHyperbola.
             getTransverseSpeedAtDeparture( ) -
             expectedValueOfTransverseSpeedAtDepartureHyperbola ) >=
             toleranceVelocity )
        {
            // Generate error statements.
            cerr << "The computed value of the tangential speed at departure ( "
                 << lambertTargeterHyperbola.
                    getTransverseSpeedAtDeparture( ) << " ) using "
                    " the Lambert targeting algorithm does not match the"
                    " expected solution (" <<
                    expectedValueOfTransverseSpeedAtDepartureHyperbola <<
                    " )." << endl;
            cerr << "The error is: "
                 << computeAbsoluteValue(
                    lambertTargeterHyperbola.
                    getTransverseSpeedAtDeparture( ) -
                    expectedValueOfTransverseSpeedAtDepartureHyperbola ) << endl;
        }

        if ( computeAbsoluteValue( lambertTargeterHyperbola.
             getTransverseSpeedAtArrival( ) -
             expectedValueOfTransverseSpeedAtArrivalHyperbola ) >=
             toleranceVelocity )
        {
            // Generate error statements.
            cerr << "The computed value of the radial speed at arrival ( "
                 << lambertTargeterHyperbola.
                    getTransverseSpeedAtArrival( ) << " ) using "
                    " the Lambert targeting algorithm does not match the"
                    " expected solution (" <<
                    expectedValueOfTransverseSpeedAtArrivalHyperbola <<
                    " )." << endl;
            cerr << "The error is: "
                 << computeAbsoluteValue(
                    lambertTargeterHyperbola.
                    getTransverseSpeedAtArrival( ) -
                    expectedValueOfTransverseSpeedAtArrivalHyperbola ) << endl;
        }

        if ( computeAbsoluteValue(
                lambertTargeterEllipse.getLambertSemiMajorAxis( ) -
                expectedValueOfSemiMajorAxisEllipse ) >=
                toleranceSemiMajorAxisEllipse )
        {
            // Generate error statements.
            cerr << "The computed value of the semi-major axis ( "
                 << lambertTargeterEllipse.
                    getLambertSemiMajorAxis( )
                 << " ) using the Lambert targeting algorithm "
                 << "does not match the expected solution ( "
                 << expectedValueOfSemiMajorAxisEllipse << " )." << endl;
            cerr << "The error is: "
                 << computeAbsoluteValue( lambertTargeterEllipse.
                    getLambertSemiMajorAxis( ) -
                    expectedValueOfSemiMajorAxisEllipse ) << endl;
        }

        if ( computeAbsoluteValue( lambertTargeterEllipse.
             getRadialSpeedAtDeparture( ) -
             expectedValueOfRadialSpeedAtDepartureEllipse ) >=
             toleranceVelocity )
        {
            // Generate error statements.
            cerr << " The computed value of the radial speed at departure ( "
                 << lambertTargeterEllipse.
                    getRadialSpeedAtDeparture( ) << " ) using "
                    " the Lambert targeting algorithm does not match the"
                    " expected solution (" <<
                    expectedValueOfRadialSpeedAtDepartureEllipse <<
                    " )." << endl;
            cerr << "The error is: "
                 << computeAbsoluteValue(
                    lambertTargeterEllipse.
                    getRadialSpeedAtDeparture( ) -
                    expectedValueOfRadialSpeedAtDepartureEllipse ) << endl;
        }

        if ( computeAbsoluteValue( lambertTargeterEllipse.
             getRadialSpeedAtArrival( ) -
             expectedValueOfRadialSpeedAtArrivalEllipse ) >=
             toleranceVelocity )
        {
            // Generate error statements.
            cerr << "The computed value of the radial speed at arrival ( "
                 << lambertTargeterEllipse.
                    getRadialSpeedAtArrival( ) << " ) using "
                    " the Lambert targeting algorithm does not match the"
                    " expected solution (" <<
                    expectedValueOfRadialSpeedAtArrivalEllipse <<
                    " )." << endl;
            cerr << "The error is: "
                 << computeAbsoluteValue(
                    lambertTargeterEllipse.
                    getRadialSpeedAtArrival( ) -
                    expectedValueOfRadialSpeedAtArrivalEllipse ) << endl;
        }

        if ( computeAbsoluteValue( lambertTargeterEllipse.
             getTransverseSpeedAtDeparture( ) -
             expectedValueOfTransverseSpeedAtDepartureEllipse ) >=
             toleranceVelocity )
        {
            // Generate error statements.
            cerr << "The computed value of the tangential speed at departure ( "
                 << lambertTargeterEllipse.
                    getTransverseSpeedAtDeparture( ) << " ) using "
                    " the Lambert targeting algorithm does not match the"
                    " expected solution (" <<
                    expectedValueOfTransverseSpeedAtDepartureEllipse <<
                    " )." << endl;
            cerr << "The error is: "
                 << computeAbsoluteValue(
                    lambertTargeterEllipse.
                    getTransverseSpeedAtDeparture( ) -
                    expectedValueOfTransverseSpeedAtDepartureEllipse ) << endl;
        }

        if ( computeAbsoluteValue( lambertTargeterEllipse.
             getTransverseSpeedAtArrival( ) -
             expectedValueOfTransverseSpeedAtArrivalEllipse ) >=
             toleranceVelocity )
        {
            // Generate error statements.
            cerr << "The computed value of the radial speed at arrival ( "
                 << lambertTargeterEllipse.
                    getTransverseSpeedAtArrival( ) << " ) using "
                    " the Lambert targeting algorithm does not match the"
                    " expected solution (" <<
                    expectedValueOfTransverseSpeedAtArrivalEllipse <<
                    " )." << endl;
            cerr << "The error is: "
                 << computeAbsoluteValue(
                    lambertTargeterEllipse.
                    getTransverseSpeedAtArrival( ) -
                    expectedValueOfTransverseSpeedAtArrivalEllipse ) << endl;
        }

        if ( positionDepartureHyperbola.cross(
                velocityDepartureHyperbola ).z( ) < 0.0 )

        {
            // Generate error statements.
            cerr << "The computed hyperbolic orbit path does not follow "
                    "the standard anti-clockwise direction." << endl;
        }

        if ( positionDepartureEllipse.cross(
                velocityDepartureEllipse ).z( ) < 0.0 )
        {
            // Generate error statements.
            cerr << "The computed elliptical orbit path does not follow "
                    "the standard anti-clockwise direction." << endl;
        }

    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isLambertTargeterErroneous;

}

}

// End of file.
