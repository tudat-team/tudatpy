/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Musegaas, P. (2012). Optimization of Space Trajectories Including Multiple Gravity Assists
 *          and Deep Space Maneuvers. MSc Thesis, Delft University of Technology, Delft,
 *          The Netherlands.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <tudat/basics/testMacros.h>

#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/mission_segments/transferNode.h"

namespace tudat
{
namespace unit_tests
{

//! Test implementation of departure leg within MGA trajectory model.
BOOST_AUTO_TEST_SUITE( test_capture_leg )

//! Test delta-V computation for infinite parking orbit.
BOOST_AUTO_TEST_CASE( testVelocitiesInfiniteParkingOrbit )
{
    // Set tolerance.
    const double tolerance = 1.0e-13;

    // Expected test result based on the last part of the ideal Messenger trajectory as modelled by
    // GTOP software distributed and downloadable from the ESA website, or within the PaGMO
    // Astrotoolbox.
    // Note that this concerns a capture at infinite distance and hence is not the best test case.
    const double expectedDeltaV = 4632.24252029314;

    // Specify the required parameters.
    // Set the planetary positions and velocities.
    const Eigen::Vector6d planetState = ( Eigen::Vector6d( ) <<
                53627979831.9492, -5044669560.01491, -5339232305.54465,
                -4857.99954791498, 50668.339570669, 4579.44661303178 ).finished( );
    std::shared_ptr< ephemerides::Ephemeris > constantEphemeris =
            std::make_shared< ephemerides::ConstantEphemeris >( planetState );

    // Set velocity before capture body.
    const Eigen::Vector3d velocityBeforePlanet (
                -5080.6362408257, 55179.1205883308, 3549.4183219232 );

    // Set the gravitational parameters. The sun's is irrelevant for the deltaV consumption, but
    // required only for other functionality in the class.
    const double mercuryGravitationalParameter = 2.2321e13;

    // Set the capture orbit (at inifinity).
    const double semiMajorAxis = std::numeric_limits< double >::infinity( );
    const double eccentricity = 0.;

    // Set test case.
    using namespace tudat::mission_segments;
    CaptureWithFixedIncomingVelocityNode captureNode(
                constantEphemeris,
                mercuryGravitationalParameter, semiMajorAxis, eccentricity,
                [=]( ){return velocityBeforePlanet; } );
    captureNode.updateNodeParameters( ( Eigen::VectorXd( 1 )<<0.0 ).finished( ) );


    // Prepare the variables for the results.
    double resultingDeltaV = captureNode.getNodeDeltaV( );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
}

//! Test delta-V computation for circular parking orbit.
BOOST_AUTO_TEST_CASE( testVelocitiesCircularParkingOrbit )
{
    // Set tolerance.
    const double tolerance = 1.0e-13;

    // Expected test result based on Table 18.2 of [Wakker, 2007]. Capture velocity to arrive at
    // a parking orbit around Mars of 1.1 times the planetary radius starting from Earth. Result
    // calculated with more precision separately using Equation 18.25 of [Wakker, 2007].
    const double expectedDeltaV = 2087.1062716740798;

    // Specify the required parameters.
    // Set the planetary positions and velocities.
    const Eigen::Vector6d planetState = ( Eigen::Vector6d( ) <<
                 227936637698.942, 0.0, 0.0,
                 0.0, 24129.4836355380, 0.0 ).finished( );
    std::shared_ptr< ephemerides::Ephemeris > constantEphemeris =
            std::make_shared< ephemerides::ConstantEphemeris >( planetState );

    // Set velocity before capture body.
    const Eigen::Vector3d velocityBeforePlanet ( 0.0, 21480.6500358053, 0.0 );

    // Set the gravitational parameters. The Sun's is irrelevant for the deltaV consumption, but
    // required only for other functionality in the class.
    const double marsGravitationalParameter = 4.2830e13;

    // Set the capture orbit.
    const double semiMajorAxis = 1.1 * 3.3895e6;;
    const double eccentricity = 0.;

    // Set test case.
    using namespace tudat::mission_segments;
    CaptureWithFixedIncomingVelocityNode captureNode(
                constantEphemeris,
                marsGravitationalParameter, semiMajorAxis, eccentricity,
                [=]( ){return velocityBeforePlanet; } );
    captureNode.updateNodeParameters( ( Eigen::VectorXd( 1 )<<0.0 ) .finished( ) );


    // Prepare the variables for the results.
    double resultingDeltaV = captureNode.getNodeDeltaV( );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
