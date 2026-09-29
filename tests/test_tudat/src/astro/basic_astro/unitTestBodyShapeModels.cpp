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

//! Create a once-subdivided icosahedron, projected onto a sphere of the given radius.
/*!
 *  Produces a closed, purely triangular surface with 42 vertices, 80 facets and 120 edges. This
 *  satisfies the Euler relations (F = 2V - 4, E = 3V - 6) that PolyhedronBodyShapeModel requires of
 *  its input. The mesh is convex and near-uniform, its facet areas lying within 16% of one another,
 *  which is what makes it usable as an oracle: for such a mesh the closest primitive to any field
 *  point is incident to the closest vertex, so the localised search performed by getAltitude
 *  reproduces the exhaustive minimum exactly. Uniformity, not convexity, is what carries this.
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

//! Compute the point of a triangle closest to the given field point.
/*!
 *  Voronoi-region algorithm of Ericson, Real-Time Collision Detection (2005), Sec. 5.1.5. The
 *  returned point lies on a vertex, an edge or the interior of the triangle, whichever is closest.
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

//! Compute the exact distance from a field point to the surface of a triangulated polyhedron.
/*!
 *  Exhaustive minimum over every facet of the distance to the closest point of that facet. Because
 *  no facet is skipped, this is independent of the neighbourhood selection performed inside
 *  PolyhedronBodyShapeModel::getAltitude and can therefore be used to verify it.
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

//! Test the selection of facets and edges used to refine the altitude of a polyhedron shape model.
/*!
 *  getAltitude does not evaluate every facet: it locates the closest vertex and then refines the
 *  distance using only the facets and edges around that vertex (Avillez, 2022). This test pins that
 *  selection down by checking it against an exhaustive minimum over the whole mesh.
 *
 *  The body is an icosphere whose facets are all of comparable size and shape. For such a mesh the
 *  closest primitive is always incident to the closest vertex, so the localised search reproduces
 *  the exhaustive minimum and any facet or edge wrongly dropped from the selection shows up as an
 *  altitude that is too large. This was checked over 40000 randomly drawn field points at altitudes
 *  between 0.0005 and 11 body radii before the point set below was fixed.
 *
 *  Convexity alone would not be enough: on a strongly anisotropic convex mesh the closest primitive
 *  need not be incident to the closest vertex, which is an "interception" in the sense of Zong et
 *  al. (2023), and getAltitude then overestimates. The near-uniform mesh used here stays out of that
 *  regime on purpose, so that the test exercises the selection logic and not the known
 *  approximation.
 *
 *  The cuboid used by testPolyhedronShapeModel is too coarse to exercise this: with eight vertices
 *  every facet touches the neighbourhood of every vertex.
 */
BOOST_AUTO_TEST_CASE( testPolyhedronShapeModelFacetAndEdgeSelection )
{
    using namespace tudat::basic_astrodynamics;

    const double radius = 1.0E3;

    Eigen::MatrixXd verticesCoordinates;
    Eigen::MatrixXi verticesDefiningEachFacet;
    createSubdividedIcosphere( radius, verticesCoordinates, verticesDefiningEachFacet );

    // Verify the mesh satisfies the relations required by the shape model, so that a later change to
    // the mesh cannot silently weaken the test.
    const int numberOfVertices = verticesCoordinates.rows( );
    BOOST_CHECK_EQUAL( numberOfVertices, 42 );
    BOOST_CHECK_EQUAL( verticesDefiningEachFacet.rows( ), 2 * numberOfVertices - 4 );

    // Distance is computed with respect to all polyhedron features, which is the path that selects
    // the facets and edges to be evaluated.
    PolyhedronBodyShapeModel shapeModel( verticesCoordinates, verticesDefiningEachFacet, false, false );

    // Build field points along the directions of every vertex, every facet centroid and every edge
    // midpoint of the mesh, at several altitudes. These are the directions along which a dropped
    // facet or edge is most likely to carry the minimum distance.
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
