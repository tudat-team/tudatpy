/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Noomen, R., Lambert targeter Excel file.
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/MissionSegments/lambertTargeter.h"

namespace tudat
{
namespace mission_segments
{

// A dummy derived class to test the functionality of the Lambert targeter base class.
class LambertTargeterDummy : public LambertTargeter
{
public:

    // Default constructor.
    LambertTargeterDummy( const Eigen::Vector3d& cartesianPositionAtDeparture,
                          const Eigen::Vector3d& cartesianPositionAtArrival,
                          const double& timeOfFlight,
                          const double& gravitationalParameter )
        : LambertTargeter( cartesianPositionAtDeparture, cartesianPositionAtArrival, timeOfFlight,
                           gravitationalParameter )
    {
        execute( );
    }

protected:

    void execute( )
    {
        // Set inertial velocities.
        cartesianVelocityAtDeparture << 2735.8, 6594.3, 0.0;
        cartesianVelocityAtArrival << -1367.9, 4225.03, 0.0;
    }

};

} // namespace mission_segments

namespace unit_tests
{

//! Test the Izzo Lambert targeting algorithm code.
BOOST_AUTO_TEST_SUITE( test_lambert_targeter )

//! Test the return of the inertial velocities.
BOOST_AUTO_TEST_CASE( testGetInertialVelocity )
{
    // Set tolerance.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    // Set canonical units for Earth (see page 29 [1]).
    const double distanceUnit = 6.378136e6;
    const double timeUnit = 806.78;

    // Set expected inertial vectors.
    const Eigen::Vector3d expectedInertialVelocityAtDeparture( 2735.8, 6594.3, 0.0 ),
            expectedInertialVelocityAtArrival( -1367.9, 4225.03, 0.0 );

    // Time conversions.
    const double testTimeOfFlight = 5.0 * timeUnit;

    // Set central body graviational parameter.
    const double testGravitationalParameter = 398600.4418e9;

    // Set position at departure and arrival.
    const Eigen::Vector3d testCartesianPositionAtDeparture( 2.0 * distanceUnit, 0.0, 0.0 ),
            testCartesianPositionAtArrival( 2.0 * distanceUnit, 2.0 * sqrt( 3.0 ) * distanceUnit,
                                            0.0 );
    // Declare Lambert targeter object.
    mission_segments::LambertTargeterDummy
            testLambertTargeter( testCartesianPositionAtDeparture, testCartesianPositionAtArrival,
                                 testTimeOfFlight, testGravitationalParameter );

    // Check that returned vectors are equal to expected vectors.
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.x( ),
                                testLambertTargeter.getInertialVelocityAtDeparture( ).x( ),
                                tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.y( ),
                                testLambertTargeter.getInertialVelocityAtDeparture( ).y( ),
                                tolerance );
    BOOST_CHECK_SMALL( testLambertTargeter.getInertialVelocityAtDeparture( ).z( ), tolerance );

    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.x( ),
                                testLambertTargeter.getInertialVelocityAtArrival( ).x( ),
                                tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.y( ),
                                testLambertTargeter.getInertialVelocityAtArrival( ).y( ),
                                tolerance );
    BOOST_CHECK_SMALL( testLambertTargeter.getInertialVelocityAtArrival( ).z( ), tolerance );

    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.x( ),
                                testLambertTargeter.getInertialVelocityVectors( ).first.x( ),
                                tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.y( ),
                                testLambertTargeter.getInertialVelocityVectors( ).first.y( ),
                                tolerance );
    BOOST_CHECK_SMALL( testLambertTargeter.getInertialVelocityVectors( ).first.z( ), tolerance );

    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.x( ),
                                testLambertTargeter.getInertialVelocityVectors( ).second.x( ),
                                tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.y( ),
                                testLambertTargeter.getInertialVelocityVectors( ).second.y( ),
                                tolerance );
    BOOST_CHECK_SMALL( testLambertTargeter.getInertialVelocityVectors( ).second.z( ), tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
