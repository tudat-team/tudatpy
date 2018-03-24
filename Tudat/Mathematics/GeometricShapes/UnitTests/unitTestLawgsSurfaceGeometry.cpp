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
 *      Craidon, C.B. A Desription of the Langley Wireframe Geometry Standard (LaWGS) format, NASA
 *          TECHNICAL MEMORANDUM 85767.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Mathematics/GeometricShapes/lawgsPartGeometry.h"
#include "Tudat/Mathematics/GeometricShapes/sphereSegment.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_Lawgs_Surface_Geometry )

//! Test implementation of Lawgs surface geometry.
BOOST_AUTO_TEST_CASE( testLawgsSurfaceGeometry )
{
    using namespace tudat;
    using namespace geometric_shapes;

    // Create a full sphere as test geometry, with a radius of 2.0.
    const double sphereRadius = 2.0;
    boost::shared_ptr< SphereSegment > sphere = boost::make_shared< SphereSegment >(
                sphereRadius );

    // Create a Lawgs mesh of the sphere.
    LawgsPartGeometry lawgsSurface;
    const int numberOfLines = 21;
    const int numberOfPoints = 21;
    lawgsSurface.setMesh( sphere, numberOfLines, numberOfPoints );

    // Retrieve the total surface area and check if it is sufficiently close
    // to the expected value.
    using mathematical_constants::PI;
    const double totalArea = lawgsSurface.getTotalArea( );
    BOOST_CHECK_SMALL( std::fabs( totalArea - 4.0 * PI
                                  * ( std::pow( sphereRadius, 2.0 ) ) ), 0.6 );

    // Test if number of lines on mesh is correct.
    BOOST_CHECK_EQUAL( lawgsSurface.getNumberOfLines( ), numberOfLines );

    // Test if number of points per line on mesh is correct.
    BOOST_CHECK_EQUAL( lawgsSurface.getNumberOfPoints( ), numberOfPoints );

    // Set part name.
    std::string partName = "sphere";
    lawgsSurface.setName( partName );

    // Test if part name is properly retrieved.
    BOOST_CHECK_EQUAL( lawgsSurface.getName( ), partName );

    // Retrieve normal and centroid for panel 0, 0.
    Eigen::Vector3d testNormal = lawgsSurface.getPanelSurfaceNormal( 0, 0 );
    Eigen::Vector3d testCentroid = lawgsSurface.getPanelCentroid( 0, 0 );

    // Test whether centroid and normal are collinear for panel 0, 0.
    BOOST_CHECK_SMALL( std::fabs( testCentroid.normalized( ).dot(
                                      testNormal.normalized( ) ) ) - 1.0, 1.0e-5 );

    // Test if the position of the x- and y-coordinate of panel 0, 0 is correct.
    BOOST_CHECK_SMALL( std::fabs( std::atan( testCentroid.y( ) / testCentroid.x( ) ) - PI / 20.0 ),
                       std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
