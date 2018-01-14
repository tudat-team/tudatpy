/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include "Tudat/JsonInterface/UnitTests/unitTestSupport.h"
#include "Tudat/JsonInterface/Propagation/referenceFrames.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_referenceFrames )

// Test 1: aerodynamic reference frames
BOOST_AUTO_TEST_CASE( test_json_referenceFrames_aerodynamic )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "aerodynamic" ),
                            reference_frames::aerodynamicsReferenceFrames,
                            reference_frames::unsupportedAerodynamicsReferenceFrames );
}

// Test 2: aerodynamic reference frame angles
BOOST_AUTO_TEST_CASE( test_json_referenceFrames_aerodynamicAngles )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "aerodynamicAngles" ),
                            reference_frames::aerodynamicsReferenceFrameAngles,
                            reference_frames::unsupportedAerodynamicsReferenceFrameAngles );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
