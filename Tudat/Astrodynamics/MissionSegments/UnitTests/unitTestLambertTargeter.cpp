/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
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
 *      120418    T. Secretin       Adapted to new class implementation.
 *      120704    P. Musegaas       Changed tolerance.
 *      120705    T. Secretin       Updated Eigen::Vector3d constructors to const-correctness.
 *
 *    References
 *      Noomen, R., Lambert targeter Excel file.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Basics/testMacros.h>

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
        cartesianVelocityAtDeparture_ << 2735.8, 6594.3, 0.0;
        cartesianVelocityAtArrival_ << -1367.9, 4225.03, 0.0;
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
    tudat::mission_segments::LambertTargeterDummy
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
