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
#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/JsonInterface/Environment/groundStations.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_groundStation )

//// Test 1: gravity field variation types
//BOOST_AUTO_TEST_CASE( test_json_groundStation_elemenType )
//{
//    BOOST_CHECK_EQUAL_ENUM( INPUT( "positionElementTypes" ),
//                            ground_stations::positionElementTypes,
//                            ground_stations::unsupportedPositionElementTypes );
//}

// Test 2: basic solid body gravity field variation
BOOST_AUTO_TEST_CASE( test_json_groundStation )
{
    using namespace simulation_setup;
    using namespace json_interface;
    using namespace unit_conversions;

    // Central body characteristics (WGS84 Earth ellipsoid).
    const double flattening = 1.0 / 298.257223563;
    const double equatorialRadius = 6378137.0;
    std::shared_ptr< basic_astrodynamics::OblateSpheroidBodyShapeModel > oblateSpheroidModel =
            std::make_shared< basic_astrodynamics::OblateSpheroidBodyShapeModel >(
                equatorialRadius, flattening );


    // Expected Cartesian and geodetic state, Montenbruck & Gill (2000) Exercise 5.3.
    const Eigen::Vector3d testCartesianPosition( 1917032.190, 6029782.349, -801376.113 );
    const Eigen::Vector3d testGeodeticPosition( -63.667,
                                                convertDegreesToRadians( -7.26654999 ),
                                                convertDegreesToRadians( 72.36312094 ) );

    // Manually compute associated spherical position
    Eigen::Vector3d testSphericalPosition = coordinate_conversions::convertCartesianToSpherical(
                testCartesianPosition );
    testSphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - testSphericalPosition( 1 );

    // Create GravityFieldVariationSettings from JSON file
    {
        const std::shared_ptr< GroundStationSettings > fromFileSettings =
                parseJSONFile< std::shared_ptr< GroundStationSettings > >( INPUT( "singleGroundStationCartesian" ) );
        const std::shared_ptr< GroundStationSettings > manualSettings =
                std::make_shared< GroundStationSettings >(
                    "station", testCartesianPosition, coordinate_conversions::cartesian_position );
        BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );

    }

    {
        const std::shared_ptr< GroundStationSettings > fromFileSettings =
                parseJSONFile< std::shared_ptr< GroundStationSettings > >( INPUT( "singleGroundStationSpherical" ) );
        const std::shared_ptr< GroundStationSettings > manualSettings =
                std::make_shared< GroundStationSettings >(
                    "station", testSphericalPosition, coordinate_conversions::spherical_position );
        BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
    }

    {
        const std::shared_ptr< GroundStationSettings > fromFileSettings =
                parseJSONFile< std::shared_ptr< GroundStationSettings > >( INPUT( "singleGroundStationGeodetic" ) );
        const std::shared_ptr< GroundStationSettings > manualSettings =
                std::make_shared< GroundStationSettings >(
                    "station", testGeodeticPosition, coordinate_conversions::geodetic_position );
        BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
