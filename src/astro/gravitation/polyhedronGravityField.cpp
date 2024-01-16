/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/gravitation/polyhedronGravityField.h"

namespace tudat
{

namespace gravitation
{

void PolyhedronGravityCache::update (const Eigen::Vector3d& currentBodyFixedPosition)
{
//    if ( currentBodyFixedPosition_.hasNaN( ) || currentBodyFixedPosition != currentBodyFixedPosition_ )
    {
        currentBodyFixedPosition_ = currentBodyFixedPosition;

        // Compute coordinates of vertices with respect to field point
        basic_mathematics::calculatePolyhedronVerticesCoordinatesRelativeToFieldPoint(
                currentVerticesCoordinatesRelativeToFieldPoint_, currentBodyFixedPosition_, verticesCoordinates_);

        // Compute per-facet factor
        basic_mathematics::calculatePolyhedronPerFacetFactor(
                currentPerFacetFactor_, currentVerticesCoordinatesRelativeToFieldPoint_, verticesDefiningEachFacet_);

        // Compute per-edge factor
        basic_mathematics::calculatePolyhedronPerEdgeFactor(
                currentPerEdgeFactor_, currentVerticesCoordinatesRelativeToFieldPoint_, verticesDefiningEachEdge_);
    }
}

void PolyhedronGravityField::computeVerticesAndFacetsDefiningEachEdge ( )
{
    const unsigned int numberOfVertices = verticesCoordinates_.rows();
    const unsigned int numberOfFacets = verticesDefiningEachFacet_.rows();
    const unsigned int numberOfEdges = 3 * ( numberOfVertices - 2 );

    verticesDefiningEachEdge_ = Eigen::MatrixXi::Constant( numberOfEdges, 2, -1 );
    facetsDefiningEachEdge_ = Eigen::MatrixXi::Constant( numberOfEdges, 2, -1 );

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
        for ( unsigned int edge = 0; edge < numberOfInsertedEdges && !edgesToInsert.empty(); ++edge )
        {
            // Loop over edges of current facet still to be inserted
            for ( unsigned int i = 0; i < edgesToInsert.size(); ++i)
            {
                if ( ( edgesToInsert.at(i).at(0) == verticesDefiningEachEdge_(edge, 0) && edgesToInsert.at(i).at(1) == verticesDefiningEachEdge_(edge, 1) ) ||
                        ( edgesToInsert.at(i).at(0) == verticesDefiningEachEdge_(edge, 1) && edgesToInsert.at(i).at(1) == verticesDefiningEachEdge_(edge, 0) ) )
                {
                    edgesToInsert.erase(edgesToInsert.begin() + i);
                    facetsDefiningEachEdge_(edge, 1) = facet;
                }
            }
        }

        // Check if any of the facet's edges still needs to be added to verticesDefiningEachEdge, and add it/them if so
        for ( unsigned int i = 0; i < edgesToInsert.size(); ++i)
        {
            verticesDefiningEachEdge_(numberOfInsertedEdges,0) = edgesToInsert.at(i).at(0);
            verticesDefiningEachEdge_(numberOfInsertedEdges,1) = edgesToInsert.at(i).at(1);
            facetsDefiningEachEdge_(numberOfInsertedEdges, 0) = facet;
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
            if ( verticesDefiningEachEdge_(i,j) < 0 )
            {
                throw std::runtime_error( "The vertices defining some edge were not selected." );
            }
            else if ( facetsDefiningEachEdge_(i,j) < 0 )
            {
                throw std::runtime_error( "The facets defining some edge were not selected." );
            }
        }
    }
}

void PolyhedronGravityField::computeFacetNormalsAndDyads ( )
{
    const unsigned int numberOfFacets = verticesDefiningEachFacet_.rows();

    facetNormalVectors_.resize(numberOfFacets);
    facetDyads_.resize(numberOfFacets);

    // Loop over facets, and for each facet compute the facet normal and the facet dyad
    for (unsigned int facet = 0; facet < numberOfFacets; ++facet)
    {
        Eigen::Vector3d vertex0 = verticesCoordinates_.block<1,3>(verticesDefiningEachFacet_(facet,0),0);
        Eigen::Vector3d vertex1 = verticesCoordinates_.block<1,3>(verticesDefiningEachFacet_(facet,1),0);
        Eigen::Vector3d vertex2 = verticesCoordinates_.block<1,3>(verticesDefiningEachFacet_(facet,2),0);

        // Compute outward-pointing vector normal to facet
        facetNormalVectors_.at(facet) = (vertex1 - vertex0).cross(vertex2 - vertex1);
        facetNormalVectors_.at(facet).normalize();

        // Compute facet dyad (outer product)
        facetDyads_.at(facet) = facetNormalVectors_.at(facet) * facetNormalVectors_.at(facet).transpose();
    }

}

void PolyhedronGravityField::computeEdgeDyads ( )
{
    const unsigned int numberOfEdges = verticesDefiningEachEdge_.rows();

    edgeDyads_.resize(numberOfEdges);

    // Loop over edges, and for each edge compute the edge dyad
    for (unsigned int edge = 0; edge < numberOfEdges; ++edge)
    {
        Eigen::Vector3d vertex0 = verticesCoordinates_.block<1,3>(verticesDefiningEachEdge_(edge,0),0);
        Eigen::Vector3d vertex1 = verticesCoordinates_.block<1,3>(verticesDefiningEachEdge_(edge,1),0);

        // Compute edge normals, with arbitrary direction
        Eigen::Vector3d facetNormalFacetA = facetNormalVectors_.at(facetsDefiningEachEdge_(edge, 0) );
        Eigen::Vector3d edgeNormalFacetA = ( vertex0 - vertex1).cross(facetNormalFacetA );
        edgeNormalFacetA.normalize();
        Eigen::Vector3d facetNormalFacetB = facetNormalVectors_.at(facetsDefiningEachEdge_(edge, 1) );
        Eigen::Vector3d edgeNormalFacetB = ( vertex0 - vertex1).cross(facetNormalFacetB );
        edgeNormalFacetB.normalize();

        // Check order of vertices for the facet associated with edgeNormalA, and correct the edge normal directions based on that
        int vertex0FacetAIndex = -1, vertex1FacetAIndex = -1;
        for (int facetVertex = 0; facetVertex < 3; ++facetVertex )
        {
            if ( verticesDefiningEachFacet_(facetsDefiningEachEdge_(edge,0), facetVertex) == verticesDefiningEachEdge_(edge, 0) )
            {
                vertex0FacetAIndex = facetVertex;
            }
            else if ( verticesDefiningEachFacet_(facetsDefiningEachEdge_(edge,0), facetVertex) == verticesDefiningEachEdge_(edge, 1) )
            {
                vertex1FacetAIndex = facetVertex;
            }
        }

        if( vertex0FacetAIndex < 0 || vertex1FacetAIndex < 0 )
        {
            throw std::runtime_error( "Error when computing edge dyads in polyhedron gravity field; could not identify vertices" );
        }


        if (( vertex0FacetAIndex == vertex1FacetAIndex + 1) || ( vertex0FacetAIndex == 0 && vertex1FacetAIndex == 2 ) )
        {
            edgeNormalFacetB = - edgeNormalFacetB;
        }
        else
        {
            edgeNormalFacetA = - edgeNormalFacetA;
        }

        // Compute edge dyads: outer product
        edgeDyads_.at(edge) = facetNormalFacetA * edgeNormalFacetA.transpose() + facetNormalFacetB * edgeNormalFacetB.transpose();

    }
}

} // namespace gravitation

} // namespace tudat
