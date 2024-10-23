/*    Copyright (c) 2010-2020, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/ephemerides/tleEphemeris.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{
namespace unit_tests
{

//! Test the functionality of the two-line elements ephemeris
BOOST_AUTO_TEST_SUITE( test_two_line_elements_ephemeris )


//! Test the functionality of the two-line elements ephemeris
BOOST_AUTO_TEST_CASE( testTwoLineElementsEphemerisVallado )
{
	using namespace tudat;
	using namespace tudat::ephemerides;

	// Dummy two line element set from Vallado (2013), page 234
	std::string elements =  "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753\n"
							"2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

	// Create TLE object from element set
	std::shared_ptr< Tle > tlePtr = std::make_shared< Tle >( elements );

	// Create TLE Ephemeris object; output state in J2000 frame fixed at the centre of the Earth
	std::shared_ptr< TleEphemeris > tleEphemeris = std::make_shared< TleEphemeris >( "Earth", "J2000", tlePtr, false );

	// Propagate TLE for 3 Julian days to be able to compare state to the one given in Vallado
    double utcEpoch =  3.0 * tudat::physical_constants::JULIAN_DAY + tlePtr->getEpoch();
    double tdbEpoch = earth_orientation::createDefaultTimeConverter( )->getCurrentTime( basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, utcEpoch );

	Eigen::Vector6d propagatedState = tleEphemeris->getCartesianState( tdbEpoch );
	Eigen::Vector3d propagatedPosition = propagatedState.head( 3 );
	Eigen::Vector3d propagatedVelocity = propagatedState.tail( 3 );

	// Position and velocity vectors as found in Vallado
	Eigen::Vector3d verificationPosition;
	Eigen::Vector3d verificationVelocity;
	verificationPosition << -9059941.3786, 4659697.2000, 813958.8875;
	verificationVelocity << -2233.348094, -4110.136162, -3157.394074;

	// Check if difference within tolerances
    BOOST_CHECK_SMALL( ( verificationPosition - propagatedPosition ).norm( ), 50.0 );
    BOOST_CHECK_SMALL( ( propagatedVelocity - verificationVelocity ).norm( ), 0.05 );


}

BOOST_AUTO_TEST_CASE( testTwoLineElementsEphemerisGehly )
{
    using namespace tudat;
    using namespace tudat::ephemerides;

    //Retrieve the initial state of Delfi-C3 using Two-Line-Elements (TLEs)
    std::string elements = "1 32789U 07021G   08119.60740078 -.00000054  00000-0  00000+0 0  9999\n"
        "2 32789 098.0082 179.6267 0015321 307.2977 051.0656 14.81417433    68";

    std::shared_ptr< Tle > delfiTle = std::make_shared< Tle >( elements );

    // Create TLE Ephemeris object; output state in J2000 frame fixed at the centre of the Earth
    std::shared_ptr< TleEphemeris > tleEphemeris = std::make_shared< TleEphemeris >( "Earth", "J2000", delfiTle, false );


    basic_astrodynamics::DateTime currentDate = basic_astrodynamics::DateTime(
        2008, 4, 28, 14, 34, 39.427392 );

    // Create time scale converter object
    auto timeScaleConverter = earth_orientation::createDefaultTimeConverter( );

    // Set the epoch in UTC scale (for instance from the above example using DateTime)
    double tdbTime = timeScaleConverter->getCurrentTime< double >(
        basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale,
        currentDate.epoch< double >( ) );

    Eigen::VectorXd currentState  = tleEphemeris->getCartesianState( tdbTime );
    Eigen::VectorXd referenceState =
        ( Eigen::Vector6d( ) << -6.99682490e+06,  3.09944881e+04, -1.90322410e+05,
        -1.99894020e+02,  1.05296158e+03,  7.47516550e+03 ).finished( );

    // Check if difference within tolerances
    BOOST_CHECK_SMALL( ( currentState - referenceState ).segment( 0, 3 ).norm( ), 50.0 );
    BOOST_CHECK_SMALL( ( currentState - referenceState ).segment( 3, 3 ).norm( ), 0.05 );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
