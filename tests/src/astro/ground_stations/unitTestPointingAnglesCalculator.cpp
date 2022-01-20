/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ground_stations/pointingAnglesCalculator.h"
#include "tudat/astro/ground_stations/groundStationState.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/simulation/environment_setup/body.h"

using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::ground_stations;
using namespace tudat;


namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_pointing_angles_calculator )

BOOST_AUTO_TEST_CASE( test_PointingAnglesCalculator )
{
    // Define Earth shape model (sphere)
    std::shared_ptr< SphericalBodyShapeModel > bodyShape = std::make_shared< SphericalBodyShapeModel >( 6.371E6 );

    // Define test ground station point.
    double groundStationDistance = 6371.0E3;
    Eigen::Vector3d groundStationPosition( groundStationDistance, 0.0, 0.0 );

    double degreesToRadians = unit_conversions::convertDegreesToRadians( 1.0 );

    // Test analytically checked azimuth and elevation
    {
        // Create ground station properties
        std::shared_ptr< GroundStationState > stationState = std::make_shared< GroundStationState >(
                    groundStationPosition, coordinate_conversions::cartesian_position, bodyShape );
        std::shared_ptr< PointingAnglesCalculator > pointAnglesCalculator = std::make_shared< PointingAnglesCalculator >(
                    [ & ]( const double ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
                    std::bind( &GroundStationState::getRotationFromBodyFixedToTopocentricFrame, stationState, std::placeholders::_1 ) );

        // Define state of viewed point
        double testLatitude = 30.0 * degreesToRadians;
        double testLongitude = 0.0;
        double testRadius = 8.0E7;
        Eigen::Vector3d testSphericalPoint( testRadius, mathematical_constants::PI / 2.0 - testLatitude, testLongitude );
        Eigen::Vector3d testCartesianPoint = coordinate_conversions::convertSphericalToCartesian( testSphericalPoint );

        // Compute azimuth/elevation angles from PointingAnglesCalculator
        double testAzimuth = pointAnglesCalculator->calculateAzimuthAngle( testCartesianPoint, 0.0 );
        double testElevation = pointAnglesCalculator->calculateElevationAngle( testCartesianPoint, 0.0 );

        double expectedAzimuth = 90.0 * degreesToRadians;
        double expectedElevation = 60.0 * degreesToRadians;

        BOOST_CHECK_CLOSE_FRACTION( expectedAzimuth, testAzimuth, 3.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( expectedElevation, testElevation, 3.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Compare results with data obtained from: http://www.movable-type.co.uk/scripts/latlong.html
    {
        {
            // Create ground station properties
            std::shared_ptr< GroundStationState > stationState = std::make_shared< GroundStationState >(
                        groundStationPosition, coordinate_conversions::cartesian_position, bodyShape );
            std::shared_ptr< PointingAnglesCalculator > pointAnglesCalculator = std::make_shared< PointingAnglesCalculator >(
                        [ & ]( const double ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
                        std::bind( &GroundStationState::getRotationFromBodyFixedToTopocentricFrame, stationState, std::placeholders::_1 ) );

            // Define state of viewed point
            double testLatitude = 21.0 * degreesToRadians;
            double testLongitude = 84.0 * degreesToRadians;
            double testRadius = 8.0E7;
            Eigen::Vector3d testSphericalPoint;
            testSphericalPoint << testRadius, mathematical_constants::PI / 2.0 - testLatitude, testLongitude;
            Eigen::Vector3d testCartesianPoint = coordinate_conversions::convertSphericalToCartesian( testSphericalPoint );

            // Compute azimuth/elevation angles from PointingAnglesCalculator
            double testAzimuth = pointAnglesCalculator->calculateAzimuthAngle( testCartesianPoint, 0.0 );
            double testElevation = pointAnglesCalculator->calculateElevationAngle( testCartesianPoint, 0.0 );

            // Set azimuth/elevation angles retrieved from website.
            double expectedElevation = mathematical_constants::PI / 2.0 - 9385.0 / 6371.0;
            double expectedAzimuth = mathematical_constants::PI / 2.0 - ( 68.0 + 53.0 / 60.0 + 40.0 / 3600.0 ) * degreesToRadians;

            BOOST_CHECK_CLOSE_FRACTION( expectedAzimuth, testAzimuth, 1.0E-5 );
            BOOST_CHECK_CLOSE_FRACTION( expectedElevation, testElevation, 1.0E-3 );
        }

        {
            // Create ground station properties
            std::shared_ptr< GroundStationState > stationState = std::make_shared< GroundStationState >(
                        groundStationPosition, coordinate_conversions::cartesian_position, bodyShape );
            std::shared_ptr< PointingAnglesCalculator > pointAnglesCalculator = std::make_shared< PointingAnglesCalculator >(
                        [ & ]( const double ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
                        std::bind( &GroundStationState::getRotationFromBodyFixedToTopocentricFrame, stationState, std::placeholders::_1 ) );

            // Define state of viewed point
            double testLatitude = -38.0 * degreesToRadians;
            double testLongitude = 234.0 * degreesToRadians;
            double testRadius = 8.0E7;
            Eigen::Vector3d testSphericalPoint;
            testSphericalPoint << testRadius, mathematical_constants::PI / 2.0 - testLatitude, testLongitude;
            Eigen::Vector3d testCartesianPoint = coordinate_conversions::convertSphericalToCartesian( testSphericalPoint );

            // Compute azimuth/elevation angles from PointingAnglesCalculator
            double testAzimuth = pointAnglesCalculator->calculateAzimuthAngle( testCartesianPoint, 0.0 );
            double testElevation = pointAnglesCalculator->calculateElevationAngle( testCartesianPoint, 0.0 );

            // Set azimuth/elevation angles retrieved from website.
            double expectedElevation = mathematical_constants::PI / 2.0 - 13080.0 / 6371.0;
            double expectedAzimuth = mathematical_constants::PI / 2.0 - ( 225.0 + 59.0 / 60.0 + 56.0 / 3600.0 ) * degreesToRadians;

            BOOST_CHECK_CLOSE_FRACTION( expectedAzimuth, testAzimuth, 1.0E-5 );
            BOOST_CHECK_CLOSE_FRACTION( expectedElevation, testElevation, 3.0E-2 );

            std::pair< double, double > pointingAngles = pointAnglesCalculator->calculatePointingAngles( testCartesianPoint, 0.0 );

            BOOST_CHECK_CLOSE_FRACTION( pointingAngles.second, testAzimuth, 1.0E-5 );
            BOOST_CHECK_CLOSE_FRACTION( pointingAngles.first, testElevation, 3.0E-2 );
        }        
    }

    // Check if inertial->topocentric rotation is handled consistently
    {
        // Define inertial->body-fixed frame rotation
        double poleRightAscension = 56.0 * degreesToRadians;
        double poleDeclination = 45.0 * degreesToRadians;
        Eigen::Quaterniond inertialToBodyFixedFrame =
                reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                    poleDeclination, poleRightAscension, 0.0 );

        // Define ground station position
        groundStationPosition = Eigen::Vector3d( 1234.0E3, -4539E3, 4298E3 );
        std::shared_ptr< GroundStationState > stationState = std::make_shared< GroundStationState >(
                    groundStationPosition, coordinate_conversions::cartesian_position, bodyShape );
        std::shared_ptr< PointingAnglesCalculator > pointAnglesCalculator = std::make_shared< PointingAnglesCalculator >(
                    [ & ]( const double ){ return inertialToBodyFixedFrame; },
                    std::bind( &GroundStationState::getRotationFromBodyFixedToTopocentricFrame, stationState, std::placeholders::_1 ) );

        // Define state of viewed point
        double testLatitude = -38.0 * degreesToRadians;
        double testLongitude = 234.0 * degreesToRadians;
        double testRadius = 8.0E7;
        Eigen::Vector3d testSphericalPoint = Eigen::Vector3d(
                    testRadius, mathematical_constants::PI / 2.0 - testLatitude, testLongitude );
        Eigen::Vector3d testCartesianPoint = coordinate_conversions::convertSphericalToCartesian( testSphericalPoint );

        // Retrieve topocentric position of viewed point from GroundStationState and PointingAnglesCalculator and compare.
        Eigen::Vector3d testPointInLocalFrame = pointAnglesCalculator->convertVectorFromInertialToTopocentricFrame(
                    testCartesianPoint, 0.0 );
        Eigen::Vector3d expectedTestPointInLocalFrame =
                stationState->getRotationFromBodyFixedToTopocentricFrame( 0.0 ) * inertialToBodyFixedFrame *
                testCartesianPoint ;
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    testPointInLocalFrame, expectedTestPointInLocalFrame, ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
