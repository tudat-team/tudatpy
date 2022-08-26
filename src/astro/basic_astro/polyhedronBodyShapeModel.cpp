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

namespace tudat
{
namespace basic_astrodynamics
{

double PolyhedronBodyShapeModel::getAltitude( const Eigen::Vector3d& bodyFixedPosition )
{
    double altitude = std::numeric_limits< double >::infinity( );

    // Compute altitude using just the distance to the vertices
    if ( justComputeDistanceToVertices_ )
    {
        altitude = computeDistanceToClosestVertex( bodyFixedPosition );
    }

    // Compute altitude using distance to vertices, facets and edges
    else
    {

    }

    // Select the altitude sign if necessary
    if ( computeAltitudeWithSign_ )
    {
        // Compute coordinates of vertices with respect to field point
        Eigen::MatrixXd verticesCoordinatesRelativeToFieldPoint;
        gravitation::calculatePolyhedronVerticesCoordinatesRelativeToFieldPoint(
                verticesCoordinatesRelativeToFieldPoint, bodyFixedPosition, verticesCoordinates_);

        // Compute per-facet factor
        Eigen::VectorXd perFacetFactor;
        gravitation::calculatePolyhedronPerFacetFactor(
                perFacetFactor, verticesCoordinatesRelativeToFieldPoint, verticesDefiningEachFacet_);

        // Compute Laplacian
        double perFacetFactorsSum = - gravitation::calculatePolyhedronLaplacianOfGravitationalPotential(
                1.0, perFacetFactor);

        // If point inside the polyhedron, altitude should be negative
        if ( perFacetFactorsSum > 2.0 * mathematical_constants::PI )
        {
            altitude = - altitude;
        }
    }

    return altitude;
}

double PolyhedronBodyShapeModel::computeDistanceToClosestVertex( const Eigen::Vector3d& bodyFixedPosition )
{
    const unsigned int numberOfVertices = verticesCoordinates_.rows();

    // Initialize distance
    double distance = std::numeric_limits< double >::infinity( );

    // Loop over vertices and select the one with smallest distance
    for (unsigned int vertex = 0; vertex < numberOfVertices; ++vertex)
    {
        double distanceToVertex = (
                ( Eigen::Vector3d() << verticesCoordinates_.block<1,3>(vertex, 0) ).finished() -
                bodyFixedPosition ).norm( );
        if ( distanceToVertex < distance )
        {
            distance = distanceToVertex;
        }
    }

    return distance;
}

double PolyhedronBodyShapeModel::computeDistanceToClosestFacet (
        const Eigen::Vector3d& bodyFixedPosition,
        const Eigen::MatrixXd& verticesDefiningEachFacetToEvaluate )
{
    const unsigned int numberOfFacets = verticesDefiningEachFacetToEvaluate.rows();

    // Initialize distance: initial value set to NAN
    double distance = TUDAT_NAN;

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
        if ( ( v0.dot(v1) >= 0 ) && ( v1.dot(v2) >= 0 ) )
        {
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

double PolyhedronBodyShapeModel::computeDistanceToClosestEdge (
        const Eigen::Vector3d& bodyFixedPosition,
        const Eigen::MatrixXd& verticesDefiningEachEdgeToEvaluate )
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

        if ( 0 < t < 1)
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

} // namespace basic_astrodynamics
} // namespace tudat

