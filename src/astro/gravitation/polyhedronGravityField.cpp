/*    Copyright (c) 2010-2019, Delft University of Technology
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
                verticesCoordinatesRelativeToFieldPoint_, currentBodyFixedPosition_, verticesCoordinates_);

        // Compute per-facet factor
        calculatePolyhedronPerFacetFactor(
                perFacetFactor_, verticesCoordinatesRelativeToFieldPoint_, verticesDefiningEachFacet_);

        // Compute per-edge factor
        calculatePolyhedronPerEdgeFactor(
                perEdgeFactor_, verticesCoordinatesRelativeToFieldPoint_, verticesDefiningEachEdge_);
    }
}


} // namespace gravitation

} // namespace tudat