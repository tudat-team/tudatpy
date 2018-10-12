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
#include "Tudat/JsonInterface/Propagation/thrust.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_acceleration )

// Test 1: thrust direction types
BOOST_AUTO_TEST_CASE( test_json_acceleration_directionTypes )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "directionTypes" ),
                            simulation_setup::thrustDirectionTypes,
                            simulation_setup::unsupportedThrustDirectionTypes );
}

// Test 2: thrust magnitude types
BOOST_AUTO_TEST_CASE( test_json_acceleration_magnitudeTypes )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "magnitudeTypes" ),
                            simulation_setup::thrustMagnitudeTypes,
                            simulation_setup::unsupportedThrustMagnitudeTypes );
}

// Test 3: thrust frame types
BOOST_AUTO_TEST_CASE( test_json_acceleration_frameTypes )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "frameTypes" ),
                            simulation_setup::thrustFrameTypes,
                            simulation_setup::unsupportedThrustFrameTypes );
}

// Test 4: thrust from direction and magnitude
BOOST_AUTO_TEST_CASE( test_json_thrust_directionMagnitude )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::json_interface;

    // Create ThrustAccelerationSettings from JSON file
    const std::shared_ptr< ThrustAccelerationSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< ThrustAccelerationSettings > >( INPUT( "directionMagnitude" ) );

    // Create ThrustAccelerationSettings manually
    const std::shared_ptr< ThrustDirectionFromStateGuidanceSettings > directionSettings =
            std::make_shared< ThrustDirectionFromStateGuidanceSettings >( "Mercury", true, false );
    const std::shared_ptr< FromBodyThrustMagnitudeSettings > magnitudeSettings =
            std::make_shared< FromBodyThrustMagnitudeSettings >( false, "booster" );
    const std::shared_ptr< ThrustAccelerationSettings > manualSettings =
            std::make_shared< ThrustAccelerationSettings >( directionSettings,
                                                              magnitudeSettings);

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 5: interpolated thrust
BOOST_AUTO_TEST_CASE( test_json_thrust_interpolated )
{
    using namespace tudat::interpolators;
    using namespace tudat::simulation_setup;
    using namespace tudat::json_interface;

    // Create ThrustAccelerationSettings from JSON file
    const std::shared_ptr< ThrustAccelerationSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< ThrustAccelerationSettings > >( INPUT( "interpolated" ) );

    // Create ThrustAccelerationSettings manually
    const std::map< double, Eigen::Vector3d > mapDependentVariables =
    {
        { 0.0,    ( Eigen::Vector3d( ) << 0.0, 0.0, 5.0 ).finished( ) },
        { 6068.0, ( Eigen::Vector3d( ) << 0.0, 1.0, 5.0 ).finished( ) },
        { 6097.0, ( Eigen::Vector3d( ) << 1.0, 0.0, 5.0 ).finished( ) }
    };
    const std::vector< Eigen::Vector3d > vectorDependentVariablesDerivatives =
    {
        ( Eigen::Vector3d( ) <<  0.0,  0.1, 0.0 ).finished( ),
        ( Eigen::Vector3d( ) <<  0.1, -0.1, 0.0 ).finished( ),
        ( Eigen::Vector3d( ) << -0.02, 0.0, 0.0 ).finished( )
    };
    const std::shared_ptr< DataInterpolationSettings< double, Eigen::Vector3d > > dataInterpolation =
            std::make_shared< DataInterpolationSettings< double, Eigen::Vector3d > >(
                std::make_shared< HermiteDataSettings< double, Eigen::Vector3d > >(
                    mapDependentVariables, vectorDependentVariablesDerivatives ),
                std::make_shared< InterpolatorSettings >( hermite_spline_interpolator ) );
    const double specificImpulse = 3000.0;
    const ThrustFrames thrustFrame = inertial_thurst_frame;
    const std::string centralBody = "Moon";
    const std::shared_ptr< ThrustAccelerationSettings > manualSettings =
            std::make_shared< ThrustAccelerationSettings >( dataInterpolation,
                                                              specificImpulse,
                                                              thrustFrame,
                                                              centralBody );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
