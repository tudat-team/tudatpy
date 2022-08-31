/*    Copyright (c) 2010-2019, Delft University of Technology
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/test/unit_test.hpp>


#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/basic_astro/polyhedronBodyShapeModel.h"
#include "tudat/astro/basic_astro/hybridBodyShapeModel.h"

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

BOOST_AUTO_TEST_CASE( testPolyhedronShapeModel )
{
    using namespace tudat::basic_astrodynamics;

    // Define tolerance
    const double tolerance = 1e-15;

    // Define cuboid polyhedron
    const double w = 10.0; // width
    const double h = 10.0; // height
    const double l = 20.0; // length

    Eigen::MatrixXd verticesCoordinates(8,3);
    Eigen::MatrixXi verticesDefiningEachFacet(12,3);
    verticesCoordinates <<
        0.0, 0.0, 0.0,
        l, 0.0, 0.0,
        0.0, w, 0.0,
        l, w, 0.0,
        0.0, 0.0, h,
        l, 0.0, h,
        0.0, w, h,
        l, w, h;
    verticesDefiningEachFacet <<
        2, 1, 0,
        1, 2, 3,
        4, 2, 0,
        2, 4, 6,
        1, 4, 0,
        4, 1, 5,
        6, 5, 7,
        5, 6, 4,
        3, 6, 7,
        6, 3, 2,
        5, 3, 7,
        3, 5, 1;

    // Test computation of altitude wrt to vertices, without sign
    {
        PolyhedronBodyShapeModel shapeModel = PolyhedronBodyShapeModel (
            verticesCoordinates, verticesDefiningEachFacet, false, true );

        Eigen::Vector3d testCartesianPosition;

        testCartesianPosition << 0.0, 0.0, 0.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );

        testCartesianPosition << 20.0, 10.0, 0.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );

        testCartesianPosition << 10.0, 0.0, 5.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt(100 + 25), tolerance );

        testCartesianPosition << 10.0, 0.0, 0.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 10.0, tolerance );

        testCartesianPosition << 10.0, 5.0, 5.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt(100 + 25 + 25), tolerance );

        testCartesianPosition << 10.0, 5.0, 20.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt(100 + 25 + 100), tolerance );
    }

    // Test computation of altitude wrt to vertices, with sign
//    {
//        PolyhedronBodyShapeModel shapeModel = PolyhedronBodyShapeModel (
//            verticesCoordinates, verticesDefiningEachFacet, true, true );
//
//        Eigen::Vector3d testCartesianPosition;
//
//        testCartesianPosition << 0.0, 0.0, 0.0;
//        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );
//
//        testCartesianPosition << 20.0, 10.0, 0.0;
//        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );
//
//        testCartesianPosition << 10.0, 0.0, 5.0;
//        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt(100 + 25), tolerance );
//
//        testCartesianPosition << 10.0, 0.0, 0.0;
//        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 10.0, tolerance );
//
//        testCartesianPosition << 10.0, 5.0, 5.0;
//        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), -std::sqrt(100 + 25 + 25), tolerance );
//
//        testCartesianPosition << 10.0, 5.0, 20.0;
//        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt(100 + 25 + 100), tolerance );
//    }

    // Test computation of altitude wrt to all polyhedron features, without sign
    {
        PolyhedronBodyShapeModel shapeModel = PolyhedronBodyShapeModel (
            verticesCoordinates, verticesDefiningEachFacet, false, false );

        Eigen::Vector3d testCartesianPosition;

        testCartesianPosition << 10.0, 5.0, 10.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );

        testCartesianPosition << 10.0, 5.0, 10.5;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.5, tolerance );

        testCartesianPosition << 10.0, 5.0, 9.5;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.5, tolerance );

        testCartesianPosition << 10.0, 0.0, 10.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );

        testCartesianPosition << 0.0, -1.0, 11.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt(2), tolerance );

        testCartesianPosition << 20.0, 10.0, 10.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );

        testCartesianPosition << 10.0, -1.0, 11.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt(2), tolerance );

    }

}

BOOST_AUTO_TEST_CASE( testHybridShapeModel )
{
    using namespace tudat::basic_astrodynamics;

    // Define tolerance
    const double tolerance = 1e-15;

    // Define switchover altitude between the two models
    const double switchoverAltitude = 5;

    // Define cuboid polyhedra. High- and low-resolution models are considered to both be cuboids,
    // but with different sizes

    const double wLowRes = 10.0; // width
    const double hLowRes = 10.0; // height
    const double lLowRes = 20.0; // length
    Eigen::MatrixXd verticesCoordinatesLowResModel( 8, 3);
    verticesCoordinatesLowResModel <<
        0.0, 0.0, 0.0,
        lLowRes, 0.0, 0.0,
        0.0, wLowRes, 0.0,
        lLowRes, wLowRes, 0.0,
        0.0, 0.0, hLowRes,
        lLowRes, 0.0, hLowRes,
        0.0, wLowRes, hLowRes,
        lLowRes, wLowRes, hLowRes;

    const double wHighRes = 5.0; // width
    const double hHighRes = 5.0; // height
    const double lHighRes = 10.0; // length
    Eigen::MatrixXd verticesCoordinatesHighResModel( 8, 3);
    verticesCoordinatesHighResModel <<
        0.0, 0.0, 0.0,
        lHighRes, 0.0, 0.0,
        0.0, wHighRes, 0.0,
        lHighRes, wHighRes, 0.0,
        0.0, 0.0, hHighRes,
        lHighRes, 0.0, hHighRes,
        0.0, wHighRes, hHighRes,
        lHighRes, wHighRes, hHighRes;

    Eigen::MatrixXi verticesDefiningEachFacet(12,3);
    verticesDefiningEachFacet <<
        2, 1, 0,
        1, 2, 3,
        4, 2, 0,
        2, 4, 6,
        1, 4, 0,
        4, 1, 5,
        6, 5, 7,
        5, 6, 4,
        3, 6, 7,
        6, 3, 2,
        5, 3, 7,
        3, 5, 1;

    std::shared_ptr< PolyhedronBodyShapeModel > highResShapeModel = std::make_shared< PolyhedronBodyShapeModel >(
    verticesCoordinatesHighResModel, verticesDefiningEachFacet, false, false );

    std::shared_ptr< PolyhedronBodyShapeModel > lowResShapeModel = std::make_shared< PolyhedronBodyShapeModel >(
        verticesCoordinatesLowResModel, verticesDefiningEachFacet, false, false );

    HybridBodyShapeModel shapeModel = HybridBodyShapeModel(
            lowResShapeModel, highResShapeModel, switchoverAltitude);

    Eigen::Vector3d testCartesianPosition;

    // Point above switchover altitude: altitude computed wrt low-resolution shape model
    testCartesianPosition << 5.0, 2.5, 20.0;
    BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 10.0, tolerance );

    // Point below switchover altitude: altitude computed wrt high-resolution shape model
    testCartesianPosition << 5.0, 2.5, 12.0;
    BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 7.0, tolerance );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

