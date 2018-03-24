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
#include "Tudat/JsonInterface/Propagation/acceleration.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_acceleration )

// Test 1: acceleration types
BOOST_AUTO_TEST_CASE( test_json_acceleration_types )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "types" ),
                            basic_astrodynamics::accelerationTypes,
                            basic_astrodynamics::unsupportedAccelerationTypes );
}

// Test 2: sphericalHarmonicGravity
BOOST_AUTO_TEST_CASE( test_json_acceleration_sphericalHarmonicGravity )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::json_interface;

    // Create AccelerationSettings from JSON file
    const boost::shared_ptr< AccelerationSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< AccelerationSettings > >( INPUT( "sphericalHarmonicGravity" ) );

    // Create AccelerationSettings manually
    const boost::shared_ptr< AccelerationSettings > manualSettings =
            boost::make_shared< SphericalHarmonicAccelerationSettings >( 7, 2 );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 3: mutualSphericalHarmonicGravity
BOOST_AUTO_TEST_CASE( test_json_acceleration_mutualSphericalHarmonicGravity )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::json_interface;

    // Create AccelerationSettings from JSON file
    const boost::shared_ptr< AccelerationSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< AccelerationSettings > >( INPUT( "mutualSphericalHarmonicGravity" ) );

    // Create AccelerationSettings manually
    const unsigned int maximumDegreeOfBodyExertingAcceleration = 7;
    const unsigned int maximumOrderOfBodyExertingAcceleration = 0;
    const unsigned int maximumDegreeOfBodyUndergoingAcceleration = 3;
    const unsigned int maximumOrderOfBodyUndergoingAcceleration = 2;
    const unsigned int maximumDegreeOfCentralBody = 5;
    const unsigned int maximumOrderOfCentralBody = 4;
    const boost::shared_ptr< AccelerationSettings > manualSettings =
            boost::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                maximumDegreeOfBodyExertingAcceleration,
                maximumOrderOfBodyExertingAcceleration,
                maximumDegreeOfBodyUndergoingAcceleration,
                maximumOrderOfBodyUndergoingAcceleration,
                maximumDegreeOfCentralBody,
                maximumOrderOfCentralBody );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 4: relativisticCorrection
BOOST_AUTO_TEST_CASE( test_json_acceleration_relativisticCorrection )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::json_interface;

    // Create AccelerationSettings from JSON file
    const boost::shared_ptr< AccelerationSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< AccelerationSettings > >( INPUT( "relativisticCorrection" ) );

    // Create AccelerationSettings manually
    const bool calculateSchwarzschildCorrection = true;
    const bool calculateLenseThirringCorrection = false;
    const bool calculateDeSitterCorrection = true;
    const std::string primaryBody = "Mars";
    const Eigen::Vector3d centralBodyAngularMomentum = ( Eigen::Vector3d( ) << 7.0E-9, 8.0E-10, 5.0E-5 ).finished( );
    const boost::shared_ptr< AccelerationSettings > manualSettings =
            boost::make_shared< RelativisticAccelerationCorrectionSettings >( calculateSchwarzschildCorrection,
                                                                              calculateLenseThirringCorrection,
                                                                              calculateDeSitterCorrection,
                                                                              primaryBody,
                                                                              centralBodyAngularMomentum );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 5: empirical
BOOST_AUTO_TEST_CASE( test_json_acceleration_empirical )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::json_interface;

    // Create AccelerationSettings from JSON file
    const boost::shared_ptr< AccelerationSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< AccelerationSettings > >( INPUT( "empirical" ) );

    // Create AccelerationSettings manually
    const Eigen::Vector3d constantAcceleration = ( Eigen::Vector3d( ) << 0.4, -0.1, 0.05 ).finished( );
    const Eigen::Vector3d sineAcceleration = ( Eigen::Vector3d( ) << 0.0, 0.02, 0.0 ).finished( );
    const Eigen::Vector3d cosineAcceleration = ( Eigen::Vector3d( ) << -0.01, 0.0, 0.0 ).finished( );
    const boost::shared_ptr< AccelerationSettings > manualSettings =
            boost::make_shared< EmpiricalAccelerationSettings >( constantAcceleration,
                                                                 sineAcceleration,
                                                                 cosineAcceleration );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
