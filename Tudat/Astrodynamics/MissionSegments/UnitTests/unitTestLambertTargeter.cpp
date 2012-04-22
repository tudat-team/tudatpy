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
 *      110113    E. Iorfida        First creation of the code.
 *      110126    E. Iorfida        Added test for velocities.
 *      110130    J. Melman         Simplified variable names, e.g., 'normOfVelocityVector' became
 *                                  'speed'. Also corrected 'tangential' to Suggested less trivial
 *                                  test case.
 *      110201    E. Iorfida        Added pointerToCelestialBody and modified variable names (from
 *                                  heliocentric, to inertial). Added elliptical case and check
 *                                  for the anti-clockwise direction of motion.
 *      110207    E. Iorfida        Added a more stric value for tolerance of semi-major axis for
 *                                  hyperbolic case.
 *      110208    E. Iorfida        Added CartesianPositionElements objects as input and
 *                                  CartesianVelocityElements objects as output.
 *      110512    K. Kumar          Updated code not to use dynamic memory allocation.
 *      110627    K. Kumar          Updated to use new predefined planets code.
 *      120416    T. Secretin       Boostified unit test.
 *
 *    References        :
 *      Mengali, G., and A.A. Quarta, Fondamenti di Meccanica del volo Spaziale.
 *      Noomen, R., Lambert targeter Excel file.
 *
 */

// Temporary notes (move to class/function doxygen):
// DISCLAIMER: At the moment, the Lambert targeter only converges for
// about half of the cases. This is not evident from the tests below, but
// it was observed during simulations carried out by the author. The
// reason might very well be an erroneous definition of the starters.
// 
// Test runs code and verifies result against expected value.
// If the tested code is erroneous, the test function returns a boolean
// true; if the code is correct, the function returns a boolean false.
// 
// The elliptical case was taken from Example 6.1, page 159-162 of
// ( Mengali, Quarta ). The hyperbolic case was taken from ( Noomen, R. ).
// 

#define BOOST_TEST_MAIN

#include <cmath>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include "Tudat/Astrodynamics/Bodies/planet.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertTargeter.h"
#include "Tudat/Astrodynamics/States/cartesianElements.h"

namespace tudat
{
namespace unit_tests
{

//! Test Lambert targeting algorithm code.
BOOST_AUTO_TEST_SUITE( test_lambert_targeter )

//! Test hyperbolic case.
BOOST_AUTO_TEST_CASE( testHyperbolicCase )
{
    // Expected test result in meters.
    // Hyperbolic test case (results from excel file [1]).
    const double expectedValueOfSemiMajorAxisHyperbola = -1270129.3602e3;
    const double expectedValueOfRadialSpeedAtDepartureHyperbola = -0.74546e3;
    const double expectedValueOfRadialSpeedAtArrivalHyperbola = 0.69321e3;
    const double expectedValueOfTransverseSpeedAtDepartureHyperbola = 0.15674e3;
    const double expectedValueOfTransverseSpeedAtArrivalHyperbola = 0.10450e3;

    // Tolerance in absolute units.
    const double toleranceSemiMajorAxisHyperbola = 1.0e2;
    const double toleranceVelocity = 1.0e-02;

    // Time conversions.
    const double timeOfFlightInDaysHyperbola = 100.0;
    const double timeOfFlightHyperbola = unit_conversions::convertJulianDaysToSeconds(
            timeOfFlightInDaysHyperbola );

    // Central bodies parameters.
    Planet predefinedEarth;
    predefinedEarth.setPredefinedPlanetSettings( Planet::earth );
    GravityFieldModel* pointerToEarthGravityField = predefinedEarth.getGravityFieldModel( );
    pointerToEarthGravityField->setGravitationalParameter( 398600.4418e9 );

    // Hyperbolic orbit case.
    LambertTargeter lambertTargeterHyperbola;

    // The starting point is twice as far as L1 and L2, which is not really
    // realistic, but it is not about the case, but about the verification.
    CartesianPositionElements positionAtDepartureHyperbola;
    positionAtDepartureHyperbola.setCartesianElementX(
                tudat::unit_conversions::convertAstronomicalUnitsToMeters( 0.02 ) );
    positionAtDepartureHyperbola.setCartesianElementY( 0.0 );
    positionAtDepartureHyperbola.setCartesianElementZ( 0.0 );

    CartesianPositionElements positionAtArrivalHyperbola;
    positionAtArrivalHyperbola.setCartesianElementX( 0.0 );
    positionAtArrivalHyperbola.setCartesianElementY(
                tudat::unit_conversions::convertAstronomicalUnitsToMeters( -0.03 ) );
    positionAtArrivalHyperbola.setCartesianElementZ( 0.0 );

    lambertTargeterHyperbola.setPositionAtDeparture( &positionAtDepartureHyperbola );
    lambertTargeterHyperbola.setPositionAtArrival( &positionAtArrivalHyperbola );
    lambertTargeterHyperbola.setNumberOfRevolutions( 0 );
    lambertTargeterHyperbola.setTimeOfFlight( timeOfFlightHyperbola );
    lambertTargeterHyperbola.setCentralBody( &predefinedEarth );

    // Create pointers to new NewtonRaphson object.
    NewtonRaphson newtonRaphsonLambertHyperbola;
    lambertTargeterHyperbola.setNewtonRaphsonMethod( &newtonRaphsonLambertHyperbola );

    // Compute Lambert targeting algorithms.
    lambertTargeterHyperbola.execute( );

    // Create local vectors for position and velocity.
    Eigen::Vector3d positionDepartureHyperbola = positionAtDepartureHyperbola.state;

    Eigen::Vector3d velocityDepartureHyperbola =
            lambertTargeterHyperbola.getInertialVelocityAtDeparture( )->state;

    // Test if the computed semi-major axis corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterHyperbola.getLambertSemiMajorAxis( ),
                                expectedValueOfSemiMajorAxisHyperbola,
                                toleranceSemiMajorAxisHyperbola );

    // Test if the computed velocity components corresponds to the expected value within the
    // specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterHyperbola.getRadialSpeedAtDeparture( ),
                                expectedValueOfRadialSpeedAtDepartureHyperbola,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterHyperbola.getRadialSpeedAtArrival( ),
                                expectedValueOfRadialSpeedAtArrivalHyperbola,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterHyperbola.getTransverseSpeedAtDeparture( ),
                                expectedValueOfTransverseSpeedAtDepartureHyperbola,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterHyperbola.getTransverseSpeedAtArrival( ),
                                expectedValueOfTransverseSpeedAtArrivalHyperbola,
                                toleranceVelocity );

    // Test if the computed solution is anti-clockwise, if the z-component of the angular momentum
    // (h = r \times v) is positive.
    BOOST_CHECK_GT( positionDepartureHyperbola.cross( velocityDepartureHyperbola ).z( ),
                    std::numeric_limits< double >::epsilon( ) );
}

//! Test elliptical case.
BOOST_AUTO_TEST_CASE( testEllipticalCase )
{
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
    double toleranceSemiMajorAxisEllipse = 1.0e4;
    double toleranceVelocity = 1.0e-02;

    // Time conversions.
    double timeOfFlightEllipse = 5.0 * timeUnit;

    // Central bodies parameters.
    Planet predefinedEarth;
    predefinedEarth.setPredefinedPlanetSettings( Planet::earth );
    GravityFieldModel* pointerToEarthGravityField = predefinedEarth.getGravityFieldModel( );
    pointerToEarthGravityField->setGravitationalParameter( 398600.4418e9 );

    // Elliptical orbit case.
    LambertTargeter lambertTargeterEllipse;

    CartesianPositionElements positionAtDepartureEllipse;
    positionAtDepartureEllipse.setCartesianElementX( 2.0 * distanceUnit );
    positionAtDepartureEllipse.setCartesianElementY( 0.0 );
    positionAtDepartureEllipse.setCartesianElementZ( 0.0 );

    CartesianPositionElements positionAtArrivalEllipse;
    positionAtArrivalEllipse.setCartesianElementX( 2.0 * distanceUnit );
    positionAtArrivalEllipse.setCartesianElementY( 2.0 * sqrt( 3.0 ) * distanceUnit );
    positionAtArrivalEllipse.setCartesianElementZ( 0.0 );

    lambertTargeterEllipse.setPositionAtDeparture( &positionAtDepartureEllipse );
    lambertTargeterEllipse.setPositionAtArrival( &positionAtArrivalEllipse );
    lambertTargeterEllipse.setNumberOfRevolutions( 0 );
    lambertTargeterEllipse.setTimeOfFlight( timeOfFlightEllipse );
    lambertTargeterEllipse.setCentralBody( &predefinedEarth );

    // Create pointers to new NewtonRaphson object.
    NewtonRaphson newtonRaphsonLambertEllipse;
    lambertTargeterEllipse.setNewtonRaphsonMethod( &newtonRaphsonLambertEllipse );

    // Compute Lambert targeting algorithms.
    lambertTargeterEllipse.execute( );

    // Create local vectors for position and velocity.
    Eigen::Vector3d positionDepartureEllipse = positionAtDepartureEllipse.state;
    Eigen::Vector3d velocityDepartureEllipse =
            lambertTargeterEllipse.getInertialVelocityAtDeparture( )->state;

    // Test if the computed semi-major axis corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterEllipse.getLambertSemiMajorAxis( ),
                                expectedValueOfSemiMajorAxisEllipse,
                                toleranceSemiMajorAxisEllipse );

    // Test if the computed velocity components corresponds to the expected value within the
    // specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterEllipse.getRadialSpeedAtDeparture( ),
                                expectedValueOfRadialSpeedAtDepartureEllipse,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterEllipse.getRadialSpeedAtArrival( ),
                                expectedValueOfRadialSpeedAtArrivalEllipse,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterEllipse.getTransverseSpeedAtDeparture( ),
                                expectedValueOfTransverseSpeedAtDepartureEllipse,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterEllipse.getTransverseSpeedAtArrival( ),
                                expectedValueOfTransverseSpeedAtArrivalEllipse,
                                toleranceVelocity );

    // Test if the computed solution is anti-clockwise, if the z-component of the angular momentum
    // (h = r \times v) is positive.
    BOOST_CHECK_GT( positionDepartureEllipse.cross( velocityDepartureEllipse ).z( ),
                    std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
