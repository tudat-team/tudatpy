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

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/ephemerides/tleEphemeris.h"

namespace tudat
{
namespace unit_tests
{

//! Test the functionality of the two-line elements ephemeris
BOOST_AUTO_TEST_SUITE( test_two_line_elements_ephemeris )


//! Test the functionality of the two-line elements ephemeris
BOOST_AUTO_TEST_CASE( testTwoLineElementsEphemeris )
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
	Eigen::Vector6d propagatedState = tleEphemeris->getCartesianState( 3.0 * tudat::physical_constants::JULIAN_DAY + tlePtr->getEpoch() );
	Eigen::Vector3d propagatedPosition = propagatedState.head( 3 );
	Eigen::Vector3d propagatedVelocity = propagatedState.tail( 3 );

	// Position and velocity vectors as found in Vallado
	Eigen::Vector3d verificationPosition;
	Eigen::Vector3d verificationVelocity;
	verificationPosition << -9059941.3786, 4659697.2000, 813958.8875;
	verificationVelocity << -2233.348094, -4110.136162, -3157.394074;

	// Check if difference within tolerances
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( propagatedPosition, verificationPosition, 5.0e-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( propagatedVelocity, verificationVelocity, 5.0e-6);

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
