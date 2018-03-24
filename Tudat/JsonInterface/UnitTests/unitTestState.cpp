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

#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.cpp"
#include "Tudat/JsonInterface/UnitTests/unitTestSupport.h"
#include "Tudat/JsonInterface/Propagation/state.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_state )

// Test 1: state types
BOOST_AUTO_TEST_CASE( test_json_state_types )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "types" ), json_interface::stateTypes, json_interface::unsupportedStateTypes );
}

// Test 2: direct Cartesian state
BOOST_AUTO_TEST_CASE( test_json_state_directCartesian )
{
    using namespace tudat::json_interface;

    const Eigen::Vector6d fromFileState = getCartesianState< double >( parseJSONFile( INPUT( "direct" ) ) );

    const Eigen::Vector6d manualState = ( Eigen::Vector6d( ) << 1.5, 0.0, 0.0, 0.0, -0.02, 0.0 ).finished( );

    BOOST_CHECK_EQUAL( fromFileState, manualState );
}

// Test 3: Cartesian state
BOOST_AUTO_TEST_CASE( test_json_state_cartesian )
{
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::json_interface;

    const Eigen::Vector6d fromFileState = getCartesianState< double >( parseJSONFile( INPUT( "cartesian" ) ) );

    Eigen::Vector6d manualState1;
    manualState1( xCartesianPositionIndex ) = 1.5;
    manualState1( yCartesianPositionIndex ) = 0.0;
    manualState1( zCartesianPositionIndex ) = 0.0;
    manualState1( xCartesianVelocityIndex ) = 0.0;
    manualState1( yCartesianVelocityIndex ) = -0.02;
    manualState1( zCartesianVelocityIndex ) = 0.0;

    BOOST_CHECK_EQUAL( fromFileState, manualState1 );
}

// Test 4: Keplerian state
BOOST_AUTO_TEST_CASE( test_json_state_keplerian )
{
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::json_interface;

    const double gravitationalParameter = 4.0E+8;


    // Semi-latus rectum

    const Eigen::Vector6d fromFileState0 = getCartesianState< double >( parseJSONFile( INPUT( "keplerian0" ) ) );

    Eigen::Vector6d manualState0 = Eigen::Vector6d::Zero( );
    manualState0( semiMajorAxisIndex ) = 3.0;
    manualState0( eccentricityIndex ) = 0.2;
    manualState0( inclinationIndex ) = 0.3;
    manualState0( argumentOfPeriapsisIndex ) = 0.4;
    manualState0( longitudeOfAscendingNodeIndex ) = 0.5;
    manualState0( trueAnomalyIndex ) = 0.6;
    manualState0 = convertKeplerianToCartesianElements( manualState0, gravitationalParameter );

    BOOST_CHECK_EQUAL( fromFileState0, manualState0 );


    // Only radius

    const Eigen::Vector6d fromFileState1 = getCartesianState< double >( parseJSONFile( INPUT( "keplerian1" ) ) );

    Eigen::Vector6d manualState1 = Eigen::Vector6d::Zero( );
    manualState1( semiMajorAxisIndex ) = 2.0;
    manualState1 = convertKeplerianToCartesianElements( manualState1, gravitationalParameter );

    BOOST_CHECK_EQUAL( fromFileState1, manualState1 );


    // Only altitude

    const Eigen::Vector6d fromFileState2 = getCartesianState< double >( parseJSONFile( INPUT( "keplerian2" ) ) );

    BOOST_CHECK_EQUAL( fromFileState2, manualState1 );


    // Only mean motion

    const Eigen::Vector6d fromFileState3 = getCartesianState< double >( parseJSONFile( INPUT( "keplerian3" ) ) );

    Eigen::Vector6d manualState3 = Eigen::Vector6d::Zero( );
    double meanMotion = 0.01;
    manualState3( semiMajorAxisIndex ) = std::pow( gravitationalParameter / std::pow( meanMotion, 2 ), 1.0/3.0 );
    manualState3 = convertKeplerianToCartesianElements( manualState3, gravitationalParameter );

    BOOST_CHECK_EQUAL( fromFileState3, manualState3 );


    // Only period

    const Eigen::Vector6d fromFileState4 = getCartesianState< double >( parseJSONFile( INPUT( "keplerian4" ) ) );

    Eigen::Vector6d manualState4 = Eigen::Vector6d::Zero( );
    meanMotion = 2.0 * mathematical_constants::PI / 0.05;
    manualState4( semiMajorAxisIndex ) = std::pow( gravitationalParameter / std::pow( meanMotion, 2 ), 1.0/3.0 );
    manualState4 = convertKeplerianToCartesianElements( manualState4, gravitationalParameter );

    BOOST_CHECK_EQUAL( fromFileState4, manualState4 );


    // Peri/apo distances

    const Eigen::Vector6d fromFileState5 = getCartesianState< double >( parseJSONFile( INPUT( "keplerian5" ) ) );

    Eigen::Vector6d manualState5 = Eigen::Vector6d::Zero( );
    const double ra = 3.0;
    const double rp = 2.0;
    manualState5( semiMajorAxisIndex ) = ( ra + rp ) / 2.0;
    manualState5( eccentricityIndex ) = ( ra - rp ) / ( ra + rp );
    manualState5( inclinationIndex ) = 0.3;
    manualState5 = convertKeplerianToCartesianElements( manualState5, gravitationalParameter );

    BOOST_CHECK_EQUAL( fromFileState5, manualState5 );


    // Peri/apo altitudes

    const Eigen::Vector6d fromFileState6 = getCartesianState< double >( parseJSONFile( INPUT( "keplerian6" ) ) );

    Eigen::Vector6d manualState6 = Eigen::Vector6d::Zero( );
    manualState6( semiMajorAxisIndex ) = ( ra + rp ) / 2.0;
    manualState6( eccentricityIndex ) = ( ra - rp ) / ( ra + rp );
    manualState6( argumentOfPeriapsisIndex ) = 0.4;
    manualState6( longitudeOfAscendingNodeIndex ) = 0.5;
    manualState6( trueAnomalyIndex ) = 0.6;
    manualState6 = convertKeplerianToCartesianElements( manualState6, gravitationalParameter );

    BOOST_CHECK_EQUAL( fromFileState6, manualState6 );


    // Semi-latus rectum

    const Eigen::Vector6d fromFileState7 = getCartesianState< double >( parseJSONFile( INPUT( "keplerian7" ) ) );

    Eigen::Vector6d manualState7 = Eigen::Vector6d::Zero( );
    manualState7( semiMajorAxisIndex ) = 3.5 / ( 1 - std::pow( 0.2, 2 ) );
    manualState7( eccentricityIndex ) = 0.2;
    manualState7( inclinationIndex ) = 0.3;
    manualState7( argumentOfPeriapsisIndex ) = 0.4;
    manualState7( longitudeOfAscendingNodeIndex ) = 0.5;
    manualState7( trueAnomalyIndex ) = 0.6;
    manualState7 = convertKeplerianToCartesianElements( manualState7, gravitationalParameter );

    BOOST_CHECK_EQUAL( fromFileState7, manualState7 );


    // Eccentric anomaly

    const Eigen::Vector6d fromFileState8 = getCartesianState< double >( parseJSONFile( INPUT( "keplerian8" ) ) );

    const double eccentricity = 0.1;
    const double eccentricAnomaly = 1.0;
    Eigen::Vector6d manualState8 = Eigen::Vector6d::Zero( );
    manualState8( semiMajorAxisIndex ) = 3.0;
    manualState8( eccentricityIndex ) = eccentricity;
    manualState8( trueAnomalyIndex ) = convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, eccentricity );
    manualState8 = convertKeplerianToCartesianElements( manualState8, gravitationalParameter );

    BOOST_CHECK_EQUAL( fromFileState8, manualState8 );


    // Mean anomaly

    const Eigen::Vector6d fromFileState9 = getCartesianState< double >( parseJSONFile( INPUT( "keplerian9" ) ) );

    const double meanAnomaly = 1.0;
    Eigen::Vector6d manualState9 = Eigen::Vector6d::Zero( );
    manualState9( semiMajorAxisIndex ) = 3.0;
    manualState9( eccentricityIndex ) = eccentricity;
    manualState9( trueAnomalyIndex ) = convertEccentricAnomalyToTrueAnomaly(
                convertMeanAnomalyToEccentricAnomaly( eccentricity, meanAnomaly ), eccentricity );
    manualState9 = convertKeplerianToCartesianElements( manualState9, gravitationalParameter );

    BOOST_CHECK_EQUAL( fromFileState9, manualState9 );

}

// Test 5: Spherical state
BOOST_AUTO_TEST_CASE( test_json_state_spherical )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::json_interface;

    // Central body

    spice_interface::loadStandardSpiceKernels( );
    boost::shared_ptr< Body > centralBody = createBodies( getDefaultBodySettings( { "Earth" } ) ).at( "Earth" );


    // From radius

    const Eigen::Vector6d fromFileState = getCartesianState< double >(
                parseJSONFile( INPUT( "spherical" ) ), KeyPath( ), centralBody );

    Eigen::Vector6d manualState;
    manualState( radiusIndex ) = 2.0;
    manualState( latitudeIndex ) = 0.5;
    manualState( longitudeIndex ) = -1.4;
    manualState( speedIndex ) = 5.0;
    manualState( flightPathIndex ) = 0.08;
    manualState( headingAngleIndex ) = 0.0;

    manualState = tudat::ephemerides::transformStateToGlobalFrame(
                convertSphericalOrbitalToCartesianState( manualState ),
                666.0, centralBody->getRotationalEphemeris( ) );

    BOOST_CHECK_EQUAL( fromFileState, manualState );


    // From altitude

    const Eigen::Vector6d fromFileStateAltitude = getCartesianState< double >(
                parseJSONFile( INPUT( "spherical_altitude" ) ), KeyPath( ), centralBody );

    Eigen::Vector6d manualStateAltitude;
    manualStateAltitude( radiusIndex ) = 2.0;
    manualStateAltitude( latitudeIndex ) = 0.5;
    manualStateAltitude( longitudeIndex ) = -1.4;
    manualStateAltitude( speedIndex ) = 5.0;
    manualStateAltitude( flightPathIndex ) = 0.08;
    manualStateAltitude( headingAngleIndex ) = -0.12;

    manualStateAltitude = tudat::ephemerides::transformStateToGlobalFrame(
                convertSphericalOrbitalToCartesianState( manualStateAltitude ),
                -8.0E+5, centralBody->getRotationalEphemeris( ) );

    BOOST_CHECK_EQUAL( fromFileStateAltitude, manualStateAltitude );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
