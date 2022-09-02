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

void calculatePolyhedronVerticesCoordinatesRelativeToFieldPoint (
        Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::Vector3d& bodyFixedPosition,
        const Eigen::MatrixXd& verticesCoordinates)
{
    const unsigned int numberOfVertices = verticesCoordinates.rows();
    verticesCoordinatesRelativeToFieldPoint.resize(numberOfVertices, 3);

    const Eigen::MatrixXd bodyFixedPositionMatrix = bodyFixedPosition.transpose();

    for ( unsigned int vertex = 0; vertex < numberOfVertices; ++vertex )
    {
        verticesCoordinatesRelativeToFieldPoint.block<1,3>(vertex,0) = verticesCoordinates.block<1,3>(vertex,0) -
                bodyFixedPositionMatrix;
    }
}

void calculatePolyhedronPerFacetFactor (
        Eigen::VectorXd& perFacetFactor,
        const Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::MatrixXi& verticesDefiningEachFacet)
{
    const unsigned int numberOfFacets = verticesDefiningEachFacet.rows();
    perFacetFactor.resize(numberOfFacets);

    for ( unsigned int facet = 0; facet < numberOfFacets; ++facet )
    {
        // Rename position vectors of each facet's vertices relative to field point
        const Eigen::Vector3d relPosI = verticesCoordinatesRelativeToFieldPoint.block<1,3>(verticesDefiningEachFacet(facet,0),0);
        const Eigen::Vector3d relPosJ = verticesCoordinatesRelativeToFieldPoint.block<1,3>(verticesDefiningEachFacet(facet,1),0);
        const Eigen::Vector3d relPosK = verticesCoordinatesRelativeToFieldPoint.block<1,3>(verticesDefiningEachFacet(facet,2),0);

        const double numerator = relPosI.dot( relPosJ.cross(relPosK) );
        if ( numerator == 0 )
        {
            perFacetFactor(facet) = 0;
        }
        else
        {
            perFacetFactor(facet) = 2.0 * atan2( numerator,
                ( relPosI.norm() * relPosJ.norm() * relPosK.norm() + relPosI.norm() * relPosJ.dot(relPosK) +
                  relPosJ.norm() * relPosK.dot(relPosI) + relPosK.norm() * relPosI.dot(relPosJ) ) );
        }

    }

}

void calculatePolyhedronPerEdgeFactor (
        Eigen::VectorXd& perEdgeFactor,
        const Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::MatrixXi& verticesDefiningEachEdge)
{
    const unsigned int numberOfEdges = verticesDefiningEachEdge.rows();
    perEdgeFactor.resize(numberOfEdges);

    for ( unsigned int edge = 0; edge < numberOfEdges; ++edge )
    {
        // Rename position vectors of each edge's vertices relative to field point
        const Eigen::Vector3d relPosI = verticesCoordinatesRelativeToFieldPoint.block<1,3>(verticesDefiningEachEdge(edge,0),0);
        const Eigen::Vector3d relPosJ = verticesCoordinatesRelativeToFieldPoint.block<1,3>(verticesDefiningEachEdge(edge,1),0);

        // Compute edge vector
        const Eigen::Vector3d eIJ = relPosI - relPosJ;

        // Selection of the edgeFactor to be 0 is only valid when computing the potential and the derivative of the
        // potential, not when computing the 2nd derivative! See "The solid angle hidden in polyhedron gravitation
        // formulations", Werner (2017), appendix C1
        const double denominator = relPosI.norm() + relPosJ.norm() - eIJ.norm();
        if ( std::abs(denominator) <  1e-18 )
        {
            perEdgeFactor(edge) = 0;
        }
        else
        {
            perEdgeFactor(edge) = log( (relPosI.norm() + relPosJ.norm() + eIJ.norm()) / denominator );
        }

    }
}

double calculatePolyhedronGravitationalPotential(
        const double gravitationalConstantTimesDensity,
        const Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::MatrixXi& verticesDefiningEachFacet,
        const Eigen::MatrixXi& verticesDefiningEachEdge,
        const std::vector< Eigen::MatrixXd >& facetDyads,
        const std::vector< Eigen::MatrixXd >& edgeDyads,
        const Eigen::VectorXd& perFacetFactor,
        const Eigen::VectorXd& perEdgeFactor )
{
    const unsigned int numberOfFacets = verticesDefiningEachFacet.rows();
    const unsigned int numberOfEdges = verticesDefiningEachEdge.rows();
    double perEdgeSum=0, perFacetSum=0;

    // Loop over edges
    for ( unsigned int edge = 0; edge < numberOfEdges; ++edge )
    {
        const Eigen::Vector3d toEdgeVector =
                verticesCoordinatesRelativeToFieldPoint.block<1,3>(verticesDefiningEachEdge(edge,0),0);

        perEdgeSum += (toEdgeVector.transpose() * edgeDyads.at(edge) * toEdgeVector * perEdgeFactor(edge) )(0,0);
    }

    // Loop over facets
    for ( unsigned int facet = 0; facet < numberOfFacets; ++facet )
    {
        const Eigen::Vector3d toFacetVector =
                verticesCoordinatesRelativeToFieldPoint.block<1,3>(verticesDefiningEachFacet(facet,0),0);

        perFacetSum += ( toFacetVector.transpose() * facetDyads.at(facet) * toFacetVector * perFacetFactor(facet) )(0,0);
    }

    return 0.5 * gravitationalConstantTimesDensity * ( perEdgeSum - perFacetSum);
}

Eigen::Vector3d calculatePolyhedronGradientOfGravitationalPotential(
        const double gravitationalConstantTimesDensity,
        const Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::MatrixXi& verticesDefiningEachFacet,
        const Eigen::MatrixXi& verticesDefiningEachEdge,
        const std::vector< Eigen::MatrixXd >& facetDyads,
        const std::vector< Eigen::MatrixXd >& edgeDyads,
        const Eigen::VectorXd& perFacetFactor,
        const Eigen::VectorXd& perEdgeFactor)
{
    const unsigned int numberOfFacets = verticesDefiningEachFacet.rows();
    const unsigned int numberOfEdges = verticesDefiningEachEdge.rows();
    Eigen::Vector3d perEdgeSum, perFacetSum;
    perEdgeSum.setZero();
    perFacetSum.setZero();

    // Loop over edges
    for ( unsigned int edge = 0; edge < numberOfEdges; ++edge )
    {
        const Eigen::Vector3d toEdgeVector = (
             verticesCoordinatesRelativeToFieldPoint.block<1,3>(verticesDefiningEachEdge(edge,0),0) +
             verticesCoordinatesRelativeToFieldPoint.block<1,3>(verticesDefiningEachEdge(edge,1),0) ) / 2;

        perEdgeSum += edgeDyads.at(edge) * toEdgeVector * perEdgeFactor(edge);
    }

    // Loop over facets
    for ( unsigned int facet = 0; facet < numberOfFacets; ++facet )
    {
        const Eigen::Vector3d toFacetVector = (
          verticesCoordinatesRelativeToFieldPoint.block<1,3>(verticesDefiningEachFacet(facet,0),0) +
          verticesCoordinatesRelativeToFieldPoint.block<1,3>(verticesDefiningEachFacet(facet,1),0) +
          verticesCoordinatesRelativeToFieldPoint.block<1,3>(verticesDefiningEachFacet(facet,2),0)
          ) / 3.0;

        perFacetSum += facetDyads.at(facet) * toFacetVector * perFacetFactor(facet);
    }

    return - gravitationalConstantTimesDensity * (perEdgeSum - perFacetSum);
}

Eigen::Matrix3d calculatePolyhedronHessianOfGravitationalPotential(
        const double gravitationalConstantTimesDensity,
        const std::vector< Eigen::MatrixXd >& facetDyads,
        const std::vector< Eigen::MatrixXd >& edgeDyads,
        const Eigen::VectorXd& perFacetFactor,
        const Eigen::VectorXd& perEdgeFactor)
{
    const unsigned int numberOfFacets = facetDyads.size();
    const unsigned int numberOfEdges = edgeDyads.size();
    Eigen::Matrix3d perEdgeSum, perFacetSum;
    perEdgeSum.setZero();
    perFacetSum.setZero();

    // Loop over edges
    for ( unsigned int edge = 0; edge < numberOfEdges; ++edge )
    {
        if ( perEdgeFactor(edge) == 0 )
        {
            // When computing the per edge factor, it is taken to be 0 at edges singularities (see function
            // calculatePolyhedronPerEdgeFactor, and reference within). This is not valid when computing the hessian matrix!
            throw std::runtime_error( "Computation of hessian matrix of singular for points at edges." );
        }
        perEdgeSum += edgeDyads.at(edge) * perEdgeFactor(edge);
    }

    // Loop over facets
    for ( unsigned int facet = 0; facet < numberOfFacets; ++facet )
    {
        perFacetSum += facetDyads.at(facet) * perFacetFactor(facet);
    }

    return gravitationalConstantTimesDensity * (perEdgeSum - perFacetSum);
}

double calculatePolyhedronLaplacianOfGravitationalPotential(
        const double gravitationalConstantTimesDensity,
        const Eigen::VectorXd& perFacetFactor)
{
    const unsigned int numberOfFacets = perFacetFactor.size();
    double perEdgeFactorSum = 0;

    for ( unsigned int facet = 0; facet < numberOfFacets; ++facet )
    {
        perEdgeFactorSum += perFacetFactor(facet);
    }

    return - gravitationalConstantTimesDensity * perEdgeFactorSum;
}

void PolyhedronGravityCache::update (const Eigen::Vector3d& currentBodyFixedPosition)
{
    if ( currentBodyFixedPosition != currentBodyFixedPosition_ )
    {
        currentBodyFixedPosition_ = currentBodyFixedPosition;

        // Compute coordinates of vertices with respect to field point
        calculatePolyhedronVerticesCoordinatesRelativeToFieldPoint(
                currentVerticesCoordinatesRelativeToFieldPoint_, currentBodyFixedPosition_, verticesCoordinates_);

        // Compute per-facet factor
        calculatePolyhedronPerFacetFactor(
                currentPerFacetFactor_, currentVerticesCoordinatesRelativeToFieldPoint_, verticesDefiningEachFacet_);

        // Compute per-edge factor
        calculatePolyhedronPerEdgeFactor(
                currentPerEdgeFactor_, currentVerticesCoordinatesRelativeToFieldPoint_, verticesDefiningEachEdge_);
    }
}

void PolyhedronGravityField::computeVerticesAndFacetsDefiningEachEdge ( )
{
    const unsigned int numberOfVertices = verticesCoordinates_.rows();
    const unsigned int numberOfFacets = verticesDefiningEachFacet_.rows();
    const unsigned int numberOfEdges = 3 * ( numberOfVertices - 2 );

    verticesDefiningEachEdge_ = Eigen::MatrixXi::Constant( numberOfEdges, 2, TUDAT_NAN );
    facetsDefiningEachEdge_ = Eigen::MatrixXi::Constant( numberOfEdges, 2, TUDAT_NAN );

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
        for ( unsigned int edge = 0; edge < numberOfInsertedEdges and !edgesToInsert.empty(); ++edge )
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
            if ( verticesDefiningEachEdge_(i,j) != verticesDefiningEachEdge_(i,j) )
            {
                throw std::runtime_error( "The vertices defining some edge were not selected." );
            }
            else if ( facetsDefiningEachEdge_(i,j) != facetsDefiningEachEdge_(i,j) )
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
        unsigned int vertex0FacetAIndex = TUDAT_NAN, vertex1FacetAIndex = TUDAT_NAN;
        for (unsigned int facetVertex = 0; facetVertex < 3; ++facetVertex )
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