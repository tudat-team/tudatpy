/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/basic_astro/polyhedronBodyShapeModel.h"
#include "tudat/astro/gravitation/polyhedronGravityField.h"

namespace tudat
{
namespace basic_astrodynamics
{

double PolyhedronBodyShapeModel::getAltitude( const Eigen::Vector3d& bodyFixedPosition )
{
    // Initialize the variable that will hold the altitude
    double altitude;

    // Compute distance to closest vertex and closest vertex id
    unsigned int closestVertex;
    double distanceToVertex = computeDistanceToClosestVertex( bodyFixedPosition, closestVertex );

    // Compute altitude using just the distance to the vertices
    if ( justComputeDistanceToVertices_ )
    {
        altitude = distanceToVertex;
    }
    // Compute altitude using distance to vertices, facets and edges
    else
    {
        const unsigned int numberOfFacets = verticesDefiningEachFacet_.rows();
        const unsigned int numberOfEdges = verticesDefiningEachEdge_.rows();

        std::vector< unsigned int > verticesToTest;
        std::vector< unsigned int > edgesToTest;
        std::vector< unsigned int > facetsToTest;

        // Select vertices connected to the closest vertex
        for (unsigned int edge = 0; edge < numberOfEdges; ++edge)
        {
            if ( closestVertex == static_cast< unsigned int >( verticesDefiningEachEdge_(edge, 0) ) )
            {
                verticesToTest.push_back( verticesDefiningEachEdge_(edge, 1) );
            }
            else if ( closestVertex == static_cast< unsigned int >( verticesDefiningEachEdge_(edge, 1) ) )
            {
                verticesToTest.push_back( verticesDefiningEachEdge_(edge, 0) );
            }
        }

        // Loop over selected vertices
        for ( unsigned int vertex : verticesToTest )
        {
            // Loop over the edges and select the ones that include the selected vertex
            for ( unsigned int edge = 0; edge < numberOfEdges; ++edge )
            {
                // If the edge isn't in the list of edges to test, check whether it should be added to it
                if ( std::count( edgesToTest.begin(), edgesToTest.end(), edge ) == 0)
                {
                    if ( static_cast< unsigned int >( verticesDefiningEachEdge_(edge, 0) ) == vertex ||
                         static_cast< unsigned int >( verticesDefiningEachEdge_(edge, 1) )  == vertex )
                    {
                        edgesToTest.push_back( edge );
                    }
                }
            }

            // Loop over the facets and select the ones that include the selected vertex
            for ( unsigned int facet = 0; facet < numberOfFacets; ++facet )
            {
                // If the edge isn't in the list of edges to test, check whether it should be added to it
                if ( std::count( facetsToTest.begin(), facetsToTest.end(), facet ) == 0)
                {
                    if ( static_cast< unsigned int >( verticesDefiningEachFacet_(facet, 0) ) == vertex ||
                         static_cast< unsigned int >( verticesDefiningEachEdge_(facet, 1) ) == vertex ||
                         static_cast< unsigned int >( verticesDefiningEachFacet_(facet, 2) ) == vertex  )
                    {
                        facetsToTest.push_back( facet );
                    }
                }
            }
        }

        // Create matrix with vertices defining each edge to be tested
        Eigen::MatrixXi verticesDefiningEachEdgeToTest = Eigen::MatrixXi::Constant(
                    edgesToTest.size( ), 2, -1 );
        unsigned int counter = 0;
        for ( unsigned int edge : edgesToTest )
        {
            verticesDefiningEachEdgeToTest(counter,0) = verticesDefiningEachEdge_(edge,0);
            verticesDefiningEachEdgeToTest(counter,1) = verticesDefiningEachEdge_(edge,1);
            ++counter;
        }

        // Create matrix with vertices defining each facet to be tested
        Eigen::MatrixXi verticesDefiningEachFacetToTest = Eigen::MatrixXi::Constant(
                    facetsToTest.size( ), 3, -1 );
        counter = 0;
        for ( unsigned int facet : facetsToTest )
        {
            verticesDefiningEachFacetToTest(counter,0) = verticesDefiningEachFacet_(facet,0);
            verticesDefiningEachFacetToTest(counter,1) = verticesDefiningEachFacet_(facet,1);
            verticesDefiningEachFacetToTest(counter,2) = verticesDefiningEachFacet_(facet,2);
            ++counter;
        }

        // Compute distance to closest edge and facet, using limited set of edges and facets
        double distanceToFacet = computeDistanceToClosestFacet(bodyFixedPosition, verticesDefiningEachFacetToTest);
        double distanceToEdge = computeDistanceToClosestEdge(bodyFixedPosition, verticesDefiningEachEdgeToTest);

        // Altitude is the minimum distance to any of the polyhedrin features
        altitude = std::min({distanceToVertex, distanceToFacet, distanceToEdge});
    }

    // Select the altitude sign if necessary
    if ( computeAltitudeWithSign_ )
    {
        // Compute coordinates of vertices with respect to field point
        Eigen::MatrixXd verticesCoordinatesRelativeToFieldPoint;
        basic_mathematics::calculatePolyhedronVerticesCoordinatesRelativeToFieldPoint(
                verticesCoordinatesRelativeToFieldPoint, bodyFixedPosition, verticesCoordinates_);

        // Compute per-facet factor
        Eigen::VectorXd perFacetFactor;
        basic_mathematics::calculatePolyhedronPerFacetFactor(
                perFacetFactor, verticesCoordinatesRelativeToFieldPoint, verticesDefiningEachFacet_);

        // Compute Laplacian
        double perFacetFactorsSum = - basic_mathematics::calculatePolyhedronLaplacianOfGravitationalPotential(
                1.0, perFacetFactor);

        // If point inside the polyhedron, altitude should be negative
        if ( perFacetFactorsSum > 2.0 * mathematical_constants::PI )
        {
            altitude = - altitude;
        }
    }

    return altitude;
}

double PolyhedronBodyShapeModel::computeDistanceToClosestVertex(
        const Eigen::Vector3d& bodyFixedPosition,
        unsigned int& closestVertexId )
{
    const unsigned int numberOfVertices = verticesCoordinates_.rows();

    // Initialize distance
    double distance = std::numeric_limits< double >::infinity( );

    // Loop over vertices and select the one with smallest distance
    for (unsigned int vertex = 0; vertex < numberOfVertices; ++vertex)
    {
        Eigen::Vector3d vertexCoordinates = verticesCoordinates_.block<1,3>(vertex, 0);
        double distanceToVertex = ( vertexCoordinates - bodyFixedPosition ).norm( );
        if ( distanceToVertex < distance )
        {
            distance = distanceToVertex;
            closestVertexId = vertex;
        }
    }

    return distance;
}

double PolyhedronBodyShapeModel::computeDistanceToClosestFacet (
        const Eigen::Vector3d& bodyFixedPosition,
        const Eigen::MatrixXi& verticesDefiningEachFacetToEvaluate )
{
    const unsigned int numberOfFacets = verticesDefiningEachFacetToEvaluate.rows();

    // Initialize distance: initial value set to NAN
    double distance = std::numeric_limits< double >::max( );

    for ( unsigned int facet = 0; facet < numberOfFacets; ++facet )
    {
        Eigen::Vector3d vertex0 = verticesCoordinates_.block<1,3>(verticesDefiningEachFacetToEvaluate(facet,0),0);
        Eigen::Vector3d vertex1 = verticesCoordinates_.block<1,3>(verticesDefiningEachFacetToEvaluate(facet,1),0);
        Eigen::Vector3d vertex2 = verticesCoordinates_.block<1,3>(verticesDefiningEachFacetToEvaluate(facet,2),0);

        // Compute outward-pointing vector normal to facet
        Eigen::Vector3d facetNormal = ((vertex1 - vertex0).cross(vertex2 - vertex1)).normalized();

        double d = (bodyFixedPosition - vertex0).dot( facetNormal );

        Eigen::Vector3d pointInPlane = bodyFixedPosition - d * facetNormal;

        Eigen::Vector3d u0 = pointInPlane - vertex0;
        Eigen::Vector3d u1 = pointInPlane - vertex1;
        Eigen::Vector3d u2 = pointInPlane - vertex2;

        Eigen::Vector3d v0 = u2.cross(u0);
        Eigen::Vector3d v1 = u0.cross(u1);
        Eigen::Vector3d v2 = u1.cross(u2);

        // If computed distance is valid
        if ( ( v0.dot(v1) > 0 ) && ( v1.dot(v2) > 0 ) )
        {
            // Compute absolute value of distance
            d = std::abs( d );
            // If d is smaller than the previous distance value: save d
            if ( d < distance )
            {
                distance = d;
            }
        }
    }

    return distance;
}

double PolyhedronBodyShapeModel::computeDistanceToClosestEdge (
        const Eigen::Vector3d& bodyFixedPosition,
        const Eigen::MatrixXi& verticesDefiningEachEdgeToEvaluate )
{
    const unsigned int numberOfEdges = verticesDefiningEachEdgeToEvaluate.rows();

    // Initialize distance: initial value set to NAN
    double distance = TUDAT_NAN;

    for (unsigned int edge = 0; edge < numberOfEdges; ++edge)
    {
        Eigen::Vector3d vertex0 = verticesCoordinates_.block<1,3>(verticesDefiningEachEdgeToEvaluate(edge,0),0);
        Eigen::Vector3d vertex1 = verticesCoordinates_.block<1,3>(verticesDefiningEachEdgeToEvaluate(edge,1),0);

        Eigen::Vector3d r_v0_p = bodyFixedPosition - vertex0;
        Eigen::Vector3d r_v0_v1 = vertex1 - vertex0;

        double d_v0_h = r_v0_p.dot( r_v0_v1.normalized() );
        double t = d_v0_h / r_v0_v1.norm();

        if ( 0 < t && t < 1)
        {
            Eigen::Vector3d h = vertex0 + t * r_v0_v1;
            double d = ( bodyFixedPosition - h ).norm();

            // If there was no previous valid distance value: save d
            if ( std::isnan( distance ) )
            {
                distance = d;
            }
            // If d is smaller than the previous distance value: save d
            else if ( d < distance )
            {
                distance = d;
            }
        }
    }

    return distance;
}

void PolyhedronBodyShapeModel::computeVerticesDefiningEachEdge( )
{
    const unsigned int numberOfVertices = verticesCoordinates_.rows();
    const unsigned int numberOfFacets = verticesDefiningEachFacet_.rows();
    const unsigned int numberOfEdges = 3 * ( numberOfVertices - 2 );

    verticesDefiningEachEdge_ = Eigen::MatrixXi::Constant( numberOfEdges, 2, -1 );

    unsigned int numberOfInsertedEdges = 0;
    for ( unsigned int facet = 0; facet < numberOfFacets; ++facet )
    {
        // Extract edges from face
        std::vector< std::vector< int > > edgesToInsert;
        const int vertex0 = verticesDefiningEachFacet_(facet,0);
        const int vertex1 = verticesDefiningEachFacet_(facet,1);
        const int vertex2 = verticesDefiningEachFacet_(facet,2);
        edgesToInsert.push_back( { vertex0, vertex1 } );
        edgesToInsert.push_back( { vertex1, vertex2 } );
        edgesToInsert.push_back( { vertex2, vertex0 } );

        // Check if all edges of the facet have been included in verticesDefiningEachEdge. If so, remove the edge from
        // the edgesToInsert vector and add the facet to facetsDefiningEachEdge
        for ( unsigned int edge = 0; ( edge < numberOfInsertedEdges && !edgesToInsert.empty() ); ++edge )
        {
            // Loop over edges of current facet still to be inserted
            for ( unsigned int i = 0; i < edgesToInsert.size(); ++i)
            {
                if ( ( edgesToInsert.at(i).at(0) == verticesDefiningEachEdge_(edge, 0) && edgesToInsert.at(i).at(1) == verticesDefiningEachEdge_(edge, 1) ) ||
                        ( edgesToInsert.at(i).at(0) == verticesDefiningEachEdge_(edge, 1) && edgesToInsert.at(i).at(1) == verticesDefiningEachEdge_(edge, 0) ) )
                {
                    edgesToInsert.erase(edgesToInsert.begin() + i);
                }
            }
        }

        // Check if any of the facet's edges still needs to be added to verticesDefiningEachEdge, and add it/them if so
        for ( unsigned int i = 0; i < edgesToInsert.size(); ++i)
        {
            verticesDefiningEachEdge_(numberOfInsertedEdges,0) = edgesToInsert.at(i).at(0);
            verticesDefiningEachEdge_(numberOfInsertedEdges,1) = edgesToInsert.at(i).at(1);
            ++numberOfInsertedEdges;
        }
    }

    // Sanity checks
    if ( numberOfInsertedEdges != numberOfEdges )
    {
        throw std::runtime_error( "Extracted number of polyhedron edges not correct." );
    }
    for ( unsigned int i = 0; i < numberOfEdges; ++i )
    {
        for (unsigned int j : {0,1} )
        {
            if ( verticesDefiningEachEdge_(i,j) != verticesDefiningEachEdge_(i,j) )
            {
                throw std::runtime_error( "The vertices defining some edge were not selected." );
            }
        }
    }
}

} // namespace basic_astrodynamics
} // namespace tudat

