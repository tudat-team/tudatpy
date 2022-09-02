/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *       "EXTERIOR GRAVITATION OF A POLYHEDRON DERIVED AND COMPARED WITH HARMONIC AND MASCON GRAVITATION REPRESENTATIONS
 *          OF ASTEROID 4769 CASTALIA", Werner and Scheeres (1997), Celestial Mechanics and Dynamical Astronomy
 *       "The solid angle hidden in polyhedron gravitation formulations", Werner (2017), Journal of Geodesy
 */

#ifndef TUDATBUNDLE_POLYHEDRON_H
#define TUDATBUNDLE_POLYHEDRON_H

#include <memory>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

namespace tudat
{
namespace basic_mathematics
{

/*! Checks if the provided polyhedron settings are valid.
 *
 * Checks if the provided polyhedron settings are valid. It verifies that the provided matrices with the vertices coordinates
 * and the vertices defining each facet have valid dimensions. Throws an error if invalid dimensions.
 *
 * @param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex, 3 columns).
 * @param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet, 3 columns).
 */
void checkValidityOfPolyhedronSettings( const Eigen::MatrixXd& verticesCoordinates,
                                        const Eigen::MatrixXi& verticesDefiningEachFacet);

/*! Computes the coordinates of the polyhedron vertices with respect a the field point.
 *
 * Computes the coordinates of the polyhedron vertices with respect to the field point.
 * @param verticesCoordinatesRelativeToFieldPoint Matrix with coordinates of each vertex wrt field point (output).
 * @param bodyFixedPosition Body fixed position of field point (input).
 * @param verticesCoordinates Coordinates of polyehdron vertices (input).
 */
void calculatePolyhedronVerticesCoordinatesRelativeToFieldPoint (
        Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::Vector3d& bodyFixedPosition,
        const Eigen::MatrixXd& verticesCoordinates);

/*! Calculates the per-facet factor of each polyhedron facet.
 *
 * Calculates the per-facet factor of each polyhedron facet, according to Eq. 27 of Werner and Scheeres (1997).
 * @param perFacetFactor Vector with the per-facet factor of each facet (output).
 * @param verticesCoordinatesRelativeToFieldPoint Matrix with coordinates of each vertex wrt field point (input).
 * @param verticesDefiningEachFacet Identification of the vertices defining each facet (0 indexed) (input)
 */
void calculatePolyhedronPerFacetFactor (
        Eigen::VectorXd& perFacetFactor,
        const Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::MatrixXi& verticesDefiningEachFacet);

/*! Calculates the per-edge factor of each polyhedron edge.
 *
 * Calculates the per-edge factor of each polyhedron edge, according to Eq. 7 of Werner and Scheeres (1997).
 * @param perEdgeFactor Vector with the per-edge factor of each edge (output).
 * @param verticesCoordinatesRelativeToFieldPoint Matrix with coordinates of each vertex wrt field point (input).
 * @param verticesDefiningEachEdge Identification of the vertices defining each facet (0 indexed) (input)
 */
void calculatePolyhedronPerEdgeFactor (
        Eigen::VectorXd& perEdgeFactor,
        const Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::MatrixXi& verticesDefiningEachEdge);

/*! Calculates the gravitational potential of a constant-density polyhedron.
 *
 * Calculates the gravitational potential of a constant-density polyhedron, according to Eq. 10 of Werner and Scheeres
 * (1997).
 * @param gravitationalConstantTimesDensity Product of the gravitational constant and density.
 * @param verticesCoordinatesRelativeToFieldPoint Matrix with coordinates of each vertex wrt field point.
 * @param verticesDefiningEachFacet Identification of the vertices defining each facet (0 indexed).
 * @param verticesDefiningEachEdge Identification of the vertices defining each edge (0 indexed).
 * @param facetDyads Vector containing facet dyads.
 * @param edgeDyads Vector containing edge dyads.
 * @param perFacetFactor Vector containing per-facet factors.
 * @param perEdgeFactor Vector containing per-edge factors.
 * @return Gravitational potential.
 */
double calculatePolyhedronGravitationalPotential(
        const double gravitationalConstantTimesDensity,
        const Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::MatrixXi& verticesDefiningEachFacet,
        const Eigen::MatrixXi& verticesDefiningEachEdge,
        const std::vector< Eigen::MatrixXd >& facetDyads,
        const std::vector< Eigen::MatrixXd >& edgeDyads,
        const Eigen::VectorXd& perFacetFactor,
        const Eigen::VectorXd& perEdgeFactor );

/*! Calculates the gradient of the potential of a constant-density polyhedron.
 *
 * Calculates the gradient of the potential of a constant-density polyhedron, according to Eq. 15 of Werner and Scheeres
 * (1997).
 * @param gravitationalConstantTimesDensity Product of the gravitational constant and density.
 * @param verticesCoordinatesRelativeToFieldPoint Matrix with coordinates of each vertex wrt field point.
 * @param verticesDefiningEachFacet Identification of the vertices defining each facet (0 indexed).
 * @param verticesDefiningEachEdge Identification of the vertices defining each edge (0 indexed).
 * @param facetDyads Vector containing facet dyads.
 * @param edgeDyads Vector containing edge dyads.
 * @param perFacetFactor Vector containing per-facet factors.
 * @param perEdgeFactor Vector containing per-edge factors.
 * @return Gradient of gravitational potential.
 */
Eigen::Vector3d calculatePolyhedronGradientOfGravitationalPotential(
        const double gravitationalConstantTimesDensity,
        const Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::MatrixXi& verticesDefiningEachFacet,
        const Eigen::MatrixXi& verticesDefiningEachEdge,
        const std::vector< Eigen::MatrixXd >& facetDyads,
        const std::vector< Eigen::MatrixXd >& edgeDyads,
        const Eigen::VectorXd& perFacetFactor,
        const Eigen::VectorXd& perEdgeFactor);

/*! Calculates the Hessian matrix of the potential of a constant-density polyhedron.
 *
 * Calculates the Hessian matrix of the potential of a constant-density polyhedron, according to Eq. 16 of Werner and
 * Scheeres (1997).
 * @param gravitationalConstantTimesDensity Product of the gravitational constant and density.
 * @param facetDyads Vector containing facet dyads.
 * @param edgeDyads Vector containing edge dyads.
 * @param perFacetFactor Vector containing per-facet factors.
 * @param perEdgeFactor Vector containing per-edge factors.
 * @return Hessian matrix of the potential.
 */
Eigen::Matrix3d calculatePolyhedronHessianOfGravitationalPotential(
        const double gravitationalConstantTimesDensity,
        const std::vector< Eigen::MatrixXd >& facetDyads,
        const std::vector< Eigen::MatrixXd >& edgeDyads,
        const Eigen::VectorXd& perFacetFactor,
        const Eigen::VectorXd& perEdgeFactor);

/*! Calculate the laplacian of the gravitational potential of a constant-density polyhedron.
 *
 * Calculate the laplacian of the gravitational potential of a constant-density polyhedron, according to Eq. 17 of
 * Werner and Scheeres (1997).
 * @param gravitationalConstantTimesDensity Product of the gravitational constant and density.
 * @param perFacetFactor Vector containing per-facet factors.
 * @return Laplacian of the gravitational potential
 */
double calculatePolyhedronLaplacianOfGravitationalPotential(
        const double gravitationalConstantTimesDensity,
        const Eigen::VectorXd& perFacetFactor);

} // namespace basic_mathematics
} // namespace tudat


#endif //TUDATBUNDLE_POLYHEDRON_H
