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

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateRepresentationConversions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_position_representation_conversions )


BOOST_AUTO_TEST_CASE( testGeneralCoordinateConversions )
{
    using namespace coordinate_conversions;
    using namespace unit_conversions;

    // Expected Cartesian state, Montenbruck & Gill (2000) Exercise 5.3.
    const Eigen::Vector3d testCartesianPosition( 1917032.190, 6029782.349, -801376.113 );

    // Expected Cartesian state, Montenbruck & Gill (2000) Exercise 5.3.
    const Eigen::Vector3d testGeodeticPosition( -63.667,
                                                convertDegreesToRadians( -7.26654999 ),
                                                convertDegreesToRadians( 72.36312094 ) );

    // Manually compute associated spherical position
    Eigen::Vector3d testSphericalPosition = coordinate_conversions::convertCartesianToSpherical(
                testCartesianPosition );
    testSphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - testSphericalPosition( 1 );


    // Central body characteristics (WGS84 Earth ellipsoid).
    const double flattening = 1.0 / 298.257223563;
    const double equatorialRadius = 6378137.0;

    // Create shape model
    boost::shared_ptr< basic_astrodynamics::OblateSpheroidBodyShapeModel > oblateSpheroidModel =
            boost::make_shared< basic_astrodynamics::OblateSpheroidBodyShapeModel >(
                equatorialRadius, flattening );

    // Declare variables for computation
    Eigen::Vector3d computedCartesianPosition, computedSphericalPosition, computedGeodeticPosition;

    // Test computation from spherical position
    {
        computedCartesianPosition = coordinate_conversions::convertPositionElements(
                    testGeodeticPosition, coordinate_conversions::geodetic_position, coordinate_conversions::cartesian_position,
                    oblateSpheroidModel );

        computedSphericalPosition = coordinate_conversions::convertPositionElements(
                    testGeodeticPosition, coordinate_conversions::geodetic_position, coordinate_conversions::spherical_position,
                    oblateSpheroidModel );


        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( testCartesianPosition( i ) - computedCartesianPosition( i ), 1.0E-3 );
        }

        BOOST_CHECK_SMALL( testSphericalPosition( 0 ) - computedSphericalPosition( 0 ), 1.0E-3 );
        BOOST_CHECK_SMALL( testSphericalPosition( 1 ) - computedSphericalPosition( 1 ), 1.0E-3 / equatorialRadius );
        BOOST_CHECK_SMALL( testSphericalPosition( 2 ) - computedSphericalPosition( 2 ), 1.0E-3 / equatorialRadius );
    }

    // Test computation from Cartesian position
    {
        computedSphericalPosition = coordinate_conversions::convertPositionElements(
                    testCartesianPosition, coordinate_conversions::cartesian_position, coordinate_conversions::spherical_position,
                    oblateSpheroidModel );

        computedGeodeticPosition = coordinate_conversions::convertPositionElements(
                    testCartesianPosition, coordinate_conversions::cartesian_position, coordinate_conversions::geodetic_position,
                    oblateSpheroidModel );

        BOOST_CHECK_SMALL( testGeodeticPosition( 0 ) - computedGeodeticPosition( 0 ), 1.0E-3 );
        BOOST_CHECK_SMALL( testGeodeticPosition( 1 ) - computedGeodeticPosition( 1 ), 1.0E-3 / equatorialRadius );
        BOOST_CHECK_SMALL( testGeodeticPosition( 2 ) - computedGeodeticPosition( 2 ), 1.0E-3 / equatorialRadius );

        BOOST_CHECK_SMALL( testSphericalPosition( 0 ) - computedSphericalPosition( 0 ), 1.0E-3 );
        BOOST_CHECK_SMALL( testSphericalPosition( 1 ) - computedSphericalPosition( 1 ), 1.0E-3 / equatorialRadius );
        BOOST_CHECK_SMALL( testSphericalPosition( 2 ) - computedSphericalPosition( 2 ), 1.0E-3 / equatorialRadius );
    }

    // Test computation from geodetic position
    {
        computedCartesianPosition = coordinate_conversions::convertPositionElements(
                    testSphericalPosition, coordinate_conversions::spherical_position, coordinate_conversions::cartesian_position,
                    oblateSpheroidModel );

        computedGeodeticPosition = coordinate_conversions::convertPositionElements(
                    testSphericalPosition, coordinate_conversions::spherical_position, coordinate_conversions::geodetic_position,
                    oblateSpheroidModel );

        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( testCartesianPosition( i ) - computedCartesianPosition( i ), 1.0E-3 );
        }

        BOOST_CHECK_SMALL( testGeodeticPosition( 0 ) - computedGeodeticPosition( 0 ), 1.0E-3 );
        BOOST_CHECK_SMALL( testGeodeticPosition( 1 ) - computedGeodeticPosition( 1 ), 1.0E-3 / equatorialRadius );
        BOOST_CHECK_SMALL( testGeodeticPosition( 2 ) - computedGeodeticPosition( 2 ), 1.0E-3 / equatorialRadius );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
