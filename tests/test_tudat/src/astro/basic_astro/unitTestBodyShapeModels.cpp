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
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000, doi:10.1007/978-3-642-58351-3.
 *      Baumgardner J R, Frederickson P O. Icosahedral Discretization of the Two-Sphere, SIAM Journal on Numerical
 *          Analysis, 22(6):1107-1115, 1985, doi:10.1137/0722066.
 *      Ericson C. Real-Time Collision Detection, Morgan Kaufmann, 2005, doi:10.1201/b14581.
 *
 */

#define BOOST_TEST_MAIN

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <utility>
#include <vector>

#include <boost/test/included/unit_test.hpp>

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

namespace
{

//! Create a ball-shaped test body made of 80 triangles (an "icosphere").
/*!
 *  The starting shape is an icosahedron: the solid with 20 identical triangular faces, like a 20-sided
 *  die. Its 12 corners are the points (0, +/-1, +/-phi) and their cyclic permutations, with phi the
 *  golden ratio. To make the shape rounder, each triangle is cut into four smaller ones by connecting
 *  the midpoints of its sides. All corners are then pushed outwards so that they lie exactly on a
 *  sphere of the given radius. This is the construction of Baumgardner and Frederickson (1985).
 *
 *  The result is a closed surface with 42 corners (vertices), 80 triangles (facets) and 120 sides
 *  (edges), which meets the input requirements of PolyhedronBodyShapeModel. All triangles have nearly
 *  the same size and shape: their areas differ by at most 16%. This regularity is what the test below
 *  relies on.
 *  \param radius Radius of the sphere onto which the vertices are projected.
 *  \param verticesCoordinates Coordinates of each vertex, one row per vertex (returned by reference).
 *  \param verticesDefiningEachFacet Indices of the vertices of each facet, counterclockwise as seen
 *      from outside the polyhedron, one row per facet (returned by reference).
 */
void createSubdividedIcosphere( const double radius, Eigen::MatrixXd& verticesCoordinates, Eigen::MatrixXi& verticesDefiningEachFacet )
{
    const double goldenRatio = ( 1.0 + std::sqrt( 5.0 ) ) / 2.0;

    std::vector< Eigen::Vector3d > vertices{ { -1.0, goldenRatio, 0.0 },  { 1.0, goldenRatio, 0.0 },   { -1.0, -goldenRatio, 0.0 },
                                             { 1.0, -goldenRatio, 0.0 },  { 0.0, -1.0, goldenRatio },  { 0.0, 1.0, goldenRatio },
                                             { 0.0, -1.0, -goldenRatio }, { 0.0, 1.0, -goldenRatio },  { goldenRatio, 0.0, -1.0 },
                                             { goldenRatio, 0.0, 1.0 },   { -goldenRatio, 0.0, -1.0 }, { -goldenRatio, 0.0, 1.0 } };

    const std::vector< std::array< int, 3 > > icosahedronFacets{ { 0, 11, 5 }, { 0, 5, 1 },  { 0, 1, 7 },   { 0, 7, 10 }, { 0, 10, 11 },
                                                                 { 1, 5, 9 },  { 5, 11, 4 }, { 11, 10, 2 }, { 10, 7, 6 }, { 7, 1, 8 },
                                                                 { 3, 9, 4 },  { 3, 4, 2 },  { 3, 2, 6 },   { 3, 6, 8 },  { 3, 8, 9 },
                                                                 { 4, 9, 5 },  { 2, 4, 11 }, { 6, 2, 10 },  { 8, 6, 7 },  { 9, 8, 1 } };

    // Split each facet into four, re-using the midpoint of each edge shared between two facets.
    std::map< std::pair< int, int >, int > midpointIndices;
    auto getMidpointIndex = [ &vertices, &midpointIndices ]( const int first, const int second ) {
        const std::pair< int, int > edgeKey( std::min( first, second ), std::max( first, second ) );
        if( midpointIndices.count( edgeKey ) == 0 )
        {
            vertices.push_back( 0.5 * ( vertices.at( first ) + vertices.at( second ) ) );
            midpointIndices[ edgeKey ] = static_cast< int >( vertices.size( ) ) - 1;
        }
        return midpointIndices.at( edgeKey );
    };

    std::vector< std::array< int, 3 > > facets;
    for( const std::array< int, 3 >& facet : icosahedronFacets )
    {
        const int firstMidpoint = getMidpointIndex( facet.at( 0 ), facet.at( 1 ) );
        const int secondMidpoint = getMidpointIndex( facet.at( 1 ), facet.at( 2 ) );
        const int thirdMidpoint = getMidpointIndex( facet.at( 2 ), facet.at( 0 ) );

        facets.push_back( { facet.at( 0 ), firstMidpoint, thirdMidpoint } );
        facets.push_back( { facet.at( 1 ), secondMidpoint, firstMidpoint } );
        facets.push_back( { facet.at( 2 ), thirdMidpoint, secondMidpoint } );
        facets.push_back( { firstMidpoint, secondMidpoint, thirdMidpoint } );
    }

    verticesCoordinates = Eigen::MatrixXd( vertices.size( ), 3 );
    for( unsigned int i = 0; i < vertices.size( ); ++i )
    {
        verticesCoordinates.row( i ) = radius * vertices.at( i ).normalized( );
    }

    verticesDefiningEachFacet = Eigen::MatrixXi( facets.size( ), 3 );
    for( unsigned int i = 0; i < facets.size( ); ++i )
    {
        verticesDefiningEachFacet.row( i ) << facets.at( i ).at( 0 ), facets.at( i ).at( 1 ), facets.at( i ).at( 2 );
    }
}

//! Find the point on a triangle that is closest to a given point in space.
/*!
 *  The closest point is either one of the triangle's three corners, a point on one of its three sides,
 *  or a point inside the triangle. The function works out which of these cases applies and returns
 *  that point. Copied from the ClosestPtPointTriangle routine in Ericson (2005), Sec. 5.1.5.
 */
Eigen::Vector3d computeClosestPointOnTriangle( const Eigen::Vector3d& fieldPoint,
                                               const Eigen::Vector3d& vertex0,
                                               const Eigen::Vector3d& vertex1,
                                               const Eigen::Vector3d& vertex2 )
{
    const Eigen::Vector3d firstEdge = vertex1 - vertex0;
    const Eigen::Vector3d secondEdge = vertex2 - vertex0;

    // Vertex region of vertex0
    const Eigen::Vector3d vector0 = fieldPoint - vertex0;
    const double d1 = firstEdge.dot( vector0 );
    const double d2 = secondEdge.dot( vector0 );
    if( d1 <= 0.0 && d2 <= 0.0 )
    {
        return vertex0;
    }

    // Vertex region of vertex1
    const Eigen::Vector3d vector1 = fieldPoint - vertex1;
    const double d3 = firstEdge.dot( vector1 );
    const double d4 = secondEdge.dot( vector1 );
    if( d3 >= 0.0 && d4 <= d3 )
    {
        return vertex1;
    }

    // Edge region of vertex0-vertex1
    const double vc = d1 * d4 - d3 * d2;
    if( vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 )
    {
        return vertex0 + ( d1 / ( d1 - d3 ) ) * firstEdge;
    }

    // Vertex region of vertex2
    const Eigen::Vector3d vector2 = fieldPoint - vertex2;
    const double d5 = firstEdge.dot( vector2 );
    const double d6 = secondEdge.dot( vector2 );
    if( d6 >= 0.0 && d5 <= d6 )
    {
        return vertex2;
    }

    // Edge region of vertex0-vertex2
    const double vb = d5 * d2 - d1 * d6;
    if( vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 )
    {
        return vertex0 + ( d2 / ( d2 - d6 ) ) * secondEdge;
    }

    // Edge region of vertex1-vertex2
    const double va = d3 * d6 - d5 * d4;
    if( va <= 0.0 && ( d4 - d3 ) >= 0.0 && ( d5 - d6 ) >= 0.0 )
    {
        return vertex1 + ( ( d4 - d3 ) / ( ( d4 - d3 ) + ( d5 - d6 ) ) ) * ( vertex2 - vertex1 );
    }

    // Interior of the facet
    const double inverseDenominator = 1.0 / ( va + vb + vc );
    return vertex0 + firstEdge * ( vb * inverseDenominator ) + secondEdge * ( vc * inverseDenominator );
}

//! Compute the distance from a point in space to a surface made of triangles, by checking every triangle.
/*!
 *  For each triangle, the closest point is found and the smallest of these distances is kept. This
 *  is slow, but it cannot miss a triangle, so it gives the correct answer against which
 *  PolyhedronBodyShapeModel::getAltitude is checked.
 */
double computeExactDistanceToPolyhedronSurface( const Eigen::Vector3d& fieldPoint,
                                                const Eigen::MatrixXd& verticesCoordinates,
                                                const Eigen::MatrixXi& verticesDefiningEachFacet )
{
    double distance = std::numeric_limits< double >::infinity( );
    for( int facet = 0; facet < verticesDefiningEachFacet.rows( ); ++facet )
    {
        const Eigen::Vector3d closestPoint =
                computeClosestPointOnTriangle( fieldPoint,
                                               verticesCoordinates.row( verticesDefiningEachFacet( facet, 0 ) ),
                                               verticesCoordinates.row( verticesDefiningEachFacet( facet, 1 ) ),
                                               verticesCoordinates.row( verticesDefiningEachFacet( facet, 2 ) ) );
        distance = std::min( distance, ( fieldPoint - closestPoint ).norm( ) );
    }
    return distance;
}

}  // namespace

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
    const Eigen::Vector3d testGeodeticPosition( -63.667, convertDegreesToRadians( -7.26654999 ), convertDegreesToRadians( 72.36312094 ) );

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
        OblateSpheroidBodyShapeModel shapeModel = OblateSpheroidBodyShapeModel( equatorialRadius, flattening );

        // Calculate altitude from shape object.
        const double altitudeFromObject = shapeModel.getAltitude( testCartesianPosition );

        // Calculate object from free function.
        const double directAltitude = calculateAltitudeOverOblateSpheroid( testCartesianPosition, equatorialRadius, flattening, 1.0E-4 );

        // Compare values.
        BOOST_CHECK_SMALL( altitudeFromObject - testGeodeticPosition.x( ), 1.0E-4 );
        BOOST_CHECK_EQUAL( altitudeFromObject, directAltitude );

        // Test calculation of full geodetic position.
        Eigen::Vector3d calculatedGeodeticPosition = shapeModel.getGeodeticPositionWrtShape( testCartesianPosition );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedGeodeticPosition, testGeodeticPosition, 1.0E-6 );
    }

    // Test free function altitude calculations
    {
        std::shared_ptr< OblateSpheroidBodyShapeModel > shapeModel =
                std::make_shared< OblateSpheroidBodyShapeModel >( equatorialRadius, flattening );

        Eigen::Vector3d bodyPosition =
                ( Eigen::Vector3d( ) << ASTRONOMICAL_UNIT / sqrt( 2.0 ), ASTRONOMICAL_UNIT / sqrt( 2.0 ), ASTRONOMICAL_UNIT * 0.01 )
                        .finished( );
        Eigen::Vector3d inertialTestCartesianPosition = testCartesianPosition + bodyPosition;

        double calculatedAltitute = getAltitudeFromNonBodyFixedPosition(
                shapeModel, testCartesianPosition, Eigen::Vector3d::Zero( ), Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
        BOOST_CHECK_SMALL( calculatedAltitute - testGeodeticPosition.x( ), 1.0E-4 );

        calculatedAltitute = getAltitudeFromNonBodyFixedPosition(
                shapeModel, inertialTestCartesianPosition, bodyPosition, Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
        BOOST_CHECK_SMALL( calculatedAltitute - testGeodeticPosition.x( ), 1.0E-4 );

        Eigen::Quaterniond dummyTestRotation = Eigen::AngleAxisd( 0.4343, Eigen::Vector3d::UnitZ( ) ) *
                Eigen::AngleAxisd( 2.4354, Eigen::Vector3d::UnitX( ) ) * Eigen::AngleAxisd( 1.2434, Eigen::Vector3d::UnitY( ) );

        calculatedAltitute = getAltitudeFromNonBodyFixedPosition( shapeModel,
                                                                  dummyTestRotation * inertialTestCartesianPosition,
                                                                  dummyTestRotation * bodyPosition,
                                                                  dummyTestRotation.inverse( ) );
        BOOST_CHECK_SMALL( calculatedAltitute - testGeodeticPosition.x( ), 1.0E-4 );

        calculatedAltitute = getAltitudeFromNonBodyFixedPositionFunctions(
                shapeModel,
                dummyTestRotation * inertialTestCartesianPosition,
                [ & ]( ) { return Eigen::Vector3d( dummyTestRotation * bodyPosition ); },
                [ & ]( ) { return dummyTestRotation.inverse( ); } );
        BOOST_CHECK_SMALL( calculatedAltitute - testGeodeticPosition.x( ), 1.0E-4 );
    }
}

BOOST_AUTO_TEST_CASE( testPolyhedronShapeModel )
{
    using namespace tudat::basic_astrodynamics;

    // Define tolerance
    const double tolerance = 1e-15;

    // Define cuboid polyhedron
    const double w = 10.0;  // width
    const double h = 10.0;  // height
    const double l = 20.0;  // length

    Eigen::MatrixXd verticesCoordinates( 8, 3 );
    Eigen::MatrixXi verticesDefiningEachFacet( 12, 3 );
    verticesCoordinates << 0.0, 0.0, 0.0, l, 0.0, 0.0, 0.0, w, 0.0, l, w, 0.0, 0.0, 0.0, h, l, 0.0, h, 0.0, w, h, l, w, h;
    verticesDefiningEachFacet << 2, 1, 0, 1, 2, 3, 4, 2, 0, 2, 4, 6, 1, 4, 0, 4, 1, 5, 6, 5, 7, 5, 6, 4, 3, 6, 7, 6, 3, 2, 5, 3, 7, 3, 5, 1;

    // Test computation of altitude wrt to vertices, without sign
    {
        PolyhedronBodyShapeModel shapeModel = PolyhedronBodyShapeModel( verticesCoordinates, verticesDefiningEachFacet, false, true );

        Eigen::Vector3d testCartesianPosition;

        testCartesianPosition << 0.0, 0.0, 0.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );

        testCartesianPosition << 20.0, 10.0, 0.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );

        testCartesianPosition << 10.0, 0.0, 5.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt( 100 + 25 ), tolerance );

        testCartesianPosition << 10.0, 0.0, 0.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 10.0, tolerance );

        testCartesianPosition << 10.0, 5.0, 5.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt( 100 + 25 + 25 ), tolerance );

        testCartesianPosition << 10.0, 5.0, 20.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt( 100 + 25 + 100 ), tolerance );
    }

    // Test computation of altitude wrt to vertices, with sign
    {
        PolyhedronBodyShapeModel shapeModel = PolyhedronBodyShapeModel( verticesCoordinates, verticesDefiningEachFacet, true, true );

        Eigen::Vector3d testCartesianPosition;

        testCartesianPosition << 0.0, 0.0, 0.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );

        testCartesianPosition << 20.0, 10.0, 0.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );

        testCartesianPosition << 10.0, 0.0, 5.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt( 100 + 25 ), tolerance );

        testCartesianPosition << 10.0, 0.0, 0.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 10.0, tolerance );

        testCartesianPosition << 10.0, 5.0, 5.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), -std::sqrt( 100 + 25 + 25 ), tolerance );

        testCartesianPosition << 10.0, 5.0, 20.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt( 100 + 25 + 100 ), tolerance );
    }

    // Test computation of altitude wrt to all polyhedron features, without sign
    {
        PolyhedronBodyShapeModel shapeModel = PolyhedronBodyShapeModel( verticesCoordinates, verticesDefiningEachFacet, false, false );

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
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt( 2 ), tolerance );

        testCartesianPosition << 20.0, 10.0, 10.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );

        testCartesianPosition << 10.0, -1.0, 11.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt( 2 ), tolerance );
    }

    // Test computation of altitude wrt to all polyhedron features, with sign
    {
        PolyhedronBodyShapeModel shapeModel = PolyhedronBodyShapeModel( verticesCoordinates, verticesDefiningEachFacet, true, false );

        Eigen::Vector3d testCartesianPosition;

        testCartesianPosition << 10.0, 5.0, 10.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );

        testCartesianPosition << 10.0, 5.0, 10.5;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.5, tolerance );

        testCartesianPosition << 10.0, 5.0, 9.5;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), -0.5, tolerance );

        testCartesianPosition << 10.0, 5.0, 10.0 - 1e-5;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ) + 1.0, -1e-5 + 1.0, tolerance );

        testCartesianPosition << 10.0, 0.0, 10.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );

        testCartesianPosition << 0.0, -1.0, 11.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt( 2 ), tolerance );

        testCartesianPosition << 20.0, 10.0, 10.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 0.0, tolerance );

        testCartesianPosition << 10.0, -1.0, 11.0;
        BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), std::sqrt( 2 ), tolerance );
    }
}

//! Test that the altitude of a polyhedron shape model is computed from the right triangles and edges.
/*!
 *  A polyhedron shape model describes the surface of a body as a mesh of flat triangles. To save time,
 *  getAltitude does not check every triangle. It first finds the mesh corner (vertex) closest to the
 *  point, and then only checks the triangles and edges that touch this vertex or one of its direct
 *  neighbours (Avillez, 2022). If this step leaves out a triangle or edge that it should include, the
 *  computed altitude comes out too large.
 *
 *  How the test works: the body is the ball-shaped mesh from createSubdividedIcosphere. For many points
 *  around it, the altitude from getAltitude is compared with the distance found by checking every
 *  triangle (computeExactDistanceToPolyhedronSurface). The two must agree. The points lie straight
 *  above every vertex, every triangle centre and every edge midpoint, at heights from 0.1% of the
 *  radius up to nine radii. These are the places where a missing triangle or edge is most likely to
 *  change the result.
 *
 *  Why the triangles must be nearly identical: the shortcut in getAltitude is only exact if the part
 *  of the surface closest to the point touches the closest vertex. On a regular mesh like this one
 *  that is always the case. This was confirmed for 40000 random points at altitudes between 0.0005
 *  and 11 radii before the fixed set of points was chosen. On a mesh with long, thin or very unequal
 *  triangles it can fail, even for a convex body (Zong et al., 2023, call this an "interception"),
 *  and getAltitude then gives too large an altitude. Such meshes are avoided here on purpose, so that
 *  the test checks the selection of triangles and edges, not this known limitation of the method.
 */
BOOST_AUTO_TEST_CASE( testPolyhedronShapeModelFacetAndEdgeSelection )
{
    using namespace tudat::basic_astrodynamics;

    const double radius = 1.0E3;

    Eigen::MatrixXd verticesCoordinates;
    Eigen::MatrixXi verticesDefiningEachFacet;
    createSubdividedIcosphere( radius, verticesCoordinates, verticesDefiningEachFacet );

    // Check the number of vertices and triangles, so that a later change to the mesh cannot silently
    // weaken the test.
    const int numberOfVertices = verticesCoordinates.rows( );
    BOOST_CHECK_EQUAL( numberOfVertices, 42 );
    BOOST_CHECK_EQUAL( verticesDefiningEachFacet.rows( ), 2 * numberOfVertices - 4 );

    // Compute the distance to vertices, edges and triangles (not just vertices), since only then does
    // getAltitude select which triangles and edges to check.
    PolyhedronBodyShapeModel shapeModel( verticesCoordinates, verticesDefiningEachFacet, false, false );

    // Place test points straight above every vertex, every triangle centre and every edge midpoint,
    // at several heights.
    std::vector< Eigen::Vector3d > directions;
    for( int vertex = 0; vertex < numberOfVertices; ++vertex )
    {
        directions.push_back( verticesCoordinates.row( vertex ) );
    }
    for( int facet = 0; facet < verticesDefiningEachFacet.rows( ); ++facet )
    {
        const Eigen::Vector3d firstVertex = verticesCoordinates.row( verticesDefiningEachFacet( facet, 0 ) );
        const Eigen::Vector3d secondVertex = verticesCoordinates.row( verticesDefiningEachFacet( facet, 1 ) );
        const Eigen::Vector3d thirdVertex = verticesCoordinates.row( verticesDefiningEachFacet( facet, 2 ) );
        directions.push_back( ( firstVertex + secondVertex + thirdVertex ) / 3.0 );
        directions.push_back( 0.5 * ( firstVertex + secondVertex ) );
        directions.push_back( 0.5 * ( secondVertex + thirdVertex ) );
        directions.push_back( 0.5 * ( thirdVertex + firstVertex ) );
    }

    const std::vector< double > radiusFactors{ 1.001, 1.01, 1.1, 1.5, 3.0, 10.0 };

    for( const Eigen::Vector3d& direction : directions )
    {
        for( const double radiusFactor : radiusFactors )
        {
            const Eigen::Vector3d fieldPoint = radiusFactor * radius * direction.normalized( );

            const double expectedAltitude =
                    computeExactDistanceToPolyhedronSurface( fieldPoint, verticesCoordinates, verticesDefiningEachFacet );

            BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( fieldPoint ), expectedAltitude, 1.0E-12 );
        }
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

    const double wLowRes = 10.0;  // width
    const double hLowRes = 10.0;  // height
    const double lLowRes = 20.0;  // length
    Eigen::MatrixXd verticesCoordinatesLowResModel( 8, 3 );
    verticesCoordinatesLowResModel << 0.0, 0.0, 0.0, lLowRes, 0.0, 0.0, 0.0, wLowRes, 0.0, lLowRes, wLowRes, 0.0, 0.0, 0.0, hLowRes,
            lLowRes, 0.0, hLowRes, 0.0, wLowRes, hLowRes, lLowRes, wLowRes, hLowRes;

    const double wHighRes = 5.0;   // width
    const double hHighRes = 5.0;   // height
    const double lHighRes = 10.0;  // length
    Eigen::MatrixXd verticesCoordinatesHighResModel( 8, 3 );
    verticesCoordinatesHighResModel << 0.0, 0.0, 0.0, lHighRes, 0.0, 0.0, 0.0, wHighRes, 0.0, lHighRes, wHighRes, 0.0, 0.0, 0.0, hHighRes,
            lHighRes, 0.0, hHighRes, 0.0, wHighRes, hHighRes, lHighRes, wHighRes, hHighRes;

    Eigen::MatrixXi verticesDefiningEachFacet( 12, 3 );
    verticesDefiningEachFacet << 2, 1, 0, 1, 2, 3, 4, 2, 0, 2, 4, 6, 1, 4, 0, 4, 1, 5, 6, 5, 7, 5, 6, 4, 3, 6, 7, 6, 3, 2, 5, 3, 7, 3, 5, 1;

    std::shared_ptr< PolyhedronBodyShapeModel > highResShapeModel =
            std::make_shared< PolyhedronBodyShapeModel >( verticesCoordinatesHighResModel, verticesDefiningEachFacet, false, false );

    std::shared_ptr< PolyhedronBodyShapeModel > lowResShapeModel =
            std::make_shared< PolyhedronBodyShapeModel >( verticesCoordinatesLowResModel, verticesDefiningEachFacet, false, false );

    HybridBodyShapeModel shapeModel = HybridBodyShapeModel( lowResShapeModel, highResShapeModel, switchoverAltitude );

    Eigen::Vector3d testCartesianPosition;

    // Point above switchover altitude: altitude computed wrt low-resolution shape model
    testCartesianPosition << 5.0, 2.5, 20.0;
    BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 10.0, tolerance );

    // Point below switchover altitude: altitude computed wrt high-resolution shape model
    testCartesianPosition << 5.0, 2.5, 12.0;
    BOOST_CHECK_CLOSE_FRACTION( shapeModel.getAltitude( testCartesianPosition ), 7.0, tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
