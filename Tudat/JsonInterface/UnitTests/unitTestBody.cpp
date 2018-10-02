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
#include "Tudat/JsonInterface/Environment/body.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_body )

// Test 1: body settings
BOOST_AUTO_TEST_CASE( test_json_body_settings )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create BodySettings from JSON file
    const std::shared_ptr< BodySettings > fromFileSettings = createBodySettings( parseJSONFile( INPUT( "body" ) ) );

    // Create BodySettings manually
    std::shared_ptr< BodySettings > manualSettings = std::make_shared< BodySettings >( );
    manualSettings->constantMass = 3000;
    manualSettings->aerodynamicCoefficientSettings = std::make_shared< ConstantAerodynamicCoefficientSettings >(
                5.0, ( Eigen::Vector3d( ) << 1.7, 0.0, 0.0 ).finished( ) );
    manualSettings->atmosphereSettings = std::make_shared< AtmosphereSettings >( nrlmsise00 );
    manualSettings->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 1.0, 0.0, -0.1, 0.0, 0.0 ).finished( ) );
    manualSettings->gravityFieldSettings = std::make_shared< CentralGravityFieldSettings >( 4.0e14 );
    const std::vector< std::string > deformingBodies = { "Moon" };
    const std::vector< std::vector< std::complex< double > > > loveNumbers =
    {
        { std::complex< double >( 1.0, 1.0 ), std::complex< double >( 2.0, 1.0 ), std::complex< double >( 3.0, 1.0 ) },
        { std::complex< double >( 0.0, 0.5 ), std::complex< double >( 0.0, 2.0 ), std::complex< double >( 0.0, 4.0 ) },
        { std::complex< double >( 0.0, 0.0 ), std::complex< double >( 0.0, 0.0 ), std::complex< double >( 0.0, 1.0 ) }
    };
    const double referenceRadius = 5e6;
    manualSettings->gravityFieldVariationSettings = { std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                                                      deformingBodies, loveNumbers, referenceRadius ) };
    manualSettings->radiationPressureSettings =
    { { "Sun", std::make_shared< CannonBallRadiationPressureInterfaceSettings >( "Sun", 5.0, 1.3 ) } };
    manualSettings->rotationModelSettings =
            std::make_shared< RotationModelSettings >( spice_rotation_model, "A", "B" );
    manualSettings->shapeModelSettings = std::make_shared< SphericalBodyShapeSettings >( 5.0e6 );

    Eigen::Vector3d testCartesianPosition( 1917032.190, 6029782.349, -801376.113 );
    Eigen::Vector3d testGeodeticPosition(
                -63.667,  unit_conversions::convertDegreesToRadians( -7.26654999 ),
                unit_conversions::convertDegreesToRadians( 72.36312094 ) );
    Eigen::Vector3d testSphericalPosition = coordinate_conversions::convertCartesianToSpherical(
                testCartesianPosition );
    testSphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - testSphericalPosition( 1 );

    manualSettings->groundStationSettings.push_back(
                std::make_shared< GroundStationSettings >(
                    "Station1", testCartesianPosition, coordinate_conversions::cartesian_position  ) );
    manualSettings->groundStationSettings.push_back(
                std::make_shared< GroundStationSettings >(
                    "Station2", testSphericalPosition, coordinate_conversions::spherical_position ) );
    manualSettings->groundStationSettings.push_back(
                std::make_shared< GroundStationSettings >(
                    "Station3", testGeodeticPosition, coordinate_conversions::geodetic_position ) );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
