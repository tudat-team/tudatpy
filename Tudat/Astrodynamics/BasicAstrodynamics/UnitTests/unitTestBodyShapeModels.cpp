/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/test/unit_test.hpp>


#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_geodetic_coordinate_conversions )

//! Test shape models, for test data, see testGeodeticCoordinateConversions.
BOOST_AUTO_TEST_CASE( testShapeModels )
{
    using namespace tudat::physical_constants;
    using namespace tudat::coordinate_conversions;
    using namespace tudat::unit_conversions;
    using namespace tudat::basic_astrodynamics;

    // Expected Cartesian state, Montenbruck & Gill (2000) Exercise 5.3.
    const Eigen::Vector3d testCartesianPosition( 1917032.190, 6029782.349, -801376.113 );

    // Expected Cartesian state, Montenbruck & Gill (2000) Exercise 5.3.
    const Eigen::Vector3d testGeodeticPosition( -63.667,
                                                convertDegreesToRadians( -7.26654999 ),
                                                convertDegreesToRadians( 72.36312094 ) );

    // Central body characteristics (WGS84 Earth ellipsoid).
    const double flattening = 1.0 / 298.257223563;
    const double equatorialRadius = 6378137.0;

    // Test sphere altitude
    {
        SphericalBodyShapeModel shapeModel = SphericalBodyShapeModel( equatorialRadius );

        // Calculate altitude from shape object.
        const double altitudeFromObject = shapeModel.getAltitude( testCartesianPosition );

        // Calculate object directly.
        const double directAltitude = testCartesianPosition.norm( ) - equatorialRadius;

        // Compare values.
        BOOST_CHECK_EQUAL( altitudeFromObject, directAltitude );
    }

    // Test oblate spheroid
    {
        OblateSpheroidBodyShapeModel shapeModel = OblateSpheroidBodyShapeModel(
                    equatorialRadius, flattening );

        // Calculate altitude from shape object.
        const double altitudeFromObject = shapeModel.getAltitude( testCartesianPosition );

        // Calculate object from free function.
        const double directAltitude = calculateAltitudeOverOblateSpheroid(
                    testCartesianPosition, equatorialRadius, flattening, 1.0E-4 );

        // Compare values.
        BOOST_CHECK_SMALL( altitudeFromObject - testGeodeticPosition.x( ), 1.0E-4 );
        BOOST_CHECK_EQUAL( altitudeFromObject, directAltitude );

        // Test calculation of full geodetic position.
        Eigen::Vector3d calculatedGeodeticPosition = shapeModel.getGeodeticPositionWrtShape(
                    testCartesianPosition );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    calculatedGeodeticPosition, testGeodeticPosition, 1.0E-6 );
    }

    // Test free function altitude calculations
    {
        std::shared_ptr< OblateSpheroidBodyShapeModel > shapeModel =
                std::make_shared< OblateSpheroidBodyShapeModel >(
                    equatorialRadius, flattening );

        Eigen::Vector3d bodyPosition =
                ( Eigen::Vector3d( ) <<
                  ASTRONOMICAL_UNIT / sqrt( 2.0 ),
                  ASTRONOMICAL_UNIT / sqrt( 2.0 ),
                  ASTRONOMICAL_UNIT * 0.01 ).finished( );
        Eigen::Vector3d inertialTestCartesianPosition = testCartesianPosition + bodyPosition;

        double calculatedAltitute = getAltitudeFromNonBodyFixedPosition(
                    shapeModel, testCartesianPosition, Eigen::Vector3d::Zero( ),
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
        BOOST_CHECK_SMALL( calculatedAltitute - testGeodeticPosition.x( ), 1.0E-4 );

        calculatedAltitute = getAltitudeFromNonBodyFixedPosition(
                    shapeModel, inertialTestCartesianPosition, bodyPosition,
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
        BOOST_CHECK_SMALL( calculatedAltitute - testGeodeticPosition.x( ), 1.0E-4 );

        Eigen::Quaterniond dummyTestRotation =
                Eigen::AngleAxisd( 0.4343, Eigen::Vector3d::UnitZ( ) ) *
                Eigen::AngleAxisd( 2.4354, Eigen::Vector3d::UnitX( ) ) *
                Eigen::AngleAxisd( 1.2434, Eigen::Vector3d::UnitY( ) );

        calculatedAltitute = getAltitudeFromNonBodyFixedPosition(
                    shapeModel, dummyTestRotation * inertialTestCartesianPosition,
                    dummyTestRotation * bodyPosition,
                    dummyTestRotation.inverse( ) );
        BOOST_CHECK_SMALL( calculatedAltitute - testGeodeticPosition.x( ), 1.0E-4 );

        calculatedAltitute = getAltitudeFromNonBodyFixedPositionFunctions(
                    shapeModel, dummyTestRotation * inertialTestCartesianPosition,
                    [ & ]( ){ return
                        Eigen::Vector3d( dummyTestRotation * bodyPosition ); },
                    [ & ]( ){ return dummyTestRotation.inverse( ); } );
        BOOST_CHECK_SMALL( calculatedAltitute - testGeodeticPosition.x( ), 1.0E-4 );

    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

