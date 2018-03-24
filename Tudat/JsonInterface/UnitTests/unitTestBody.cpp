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
    const boost::shared_ptr< BodySettings > fromFileSettings = createBodySettings( parseJSONFile( INPUT( "body" ) ) );

    // Create BodySettings manually
    boost::shared_ptr< BodySettings > manualSettings = boost::make_shared< BodySettings >( );
    manualSettings->constantMass = 3000;
    manualSettings->aerodynamicCoefficientSettings = boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                5.0, ( Eigen::Vector3d( ) << 1.7, 0.0, 0.0 ).finished( ) );
    manualSettings->atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 );
    manualSettings->ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 1.0, 0.0, -0.1, 0.0, 0.0 ).finished( ) );
    manualSettings->gravityFieldSettings = boost::make_shared< CentralGravityFieldSettings >( 4.0e14 );
    const std::vector< std::string > deformingBodies = { "Moon" };
    const std::vector< std::vector< std::complex< double > > > loveNumbers =
    {
        { std::complex< double >( 1.0, 1.0 ), std::complex< double >( 2.0, 1.0 ), std::complex< double >( 3.0, 1.0 ) },
        { std::complex< double >( 0.0, 0.5 ), std::complex< double >( 0.0, 2.0 ), std::complex< double >( 0.0, 4.0 ) },
        { std::complex< double >( 0.0, 0.0 ), std::complex< double >( 0.0, 0.0 ), std::complex< double >( 0.0, 1.0 ) }
    };
    const double referenceRadius = 5e6;
    manualSettings->gravityFieldVariationSettings = { boost::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                                                      deformingBodies, loveNumbers, referenceRadius ) };
    manualSettings->radiationPressureSettings =
    { { "Sun", boost::make_shared< CannonBallRadiationPressureInterfaceSettings >( "Sun", 5.0, 1.3 ) } };
    manualSettings->rotationModelSettings =
            boost::make_shared< RotationModelSettings >( spice_rotation_model, "A", "B" );
    manualSettings->shapeModelSettings = boost::make_shared< SphericalBodyShapeSettings >( 5.0e6 );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
