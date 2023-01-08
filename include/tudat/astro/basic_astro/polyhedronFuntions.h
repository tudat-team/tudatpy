/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      A. Dobrovolskis (1996), "Inertia of Any Polyhedron", Icarus, 124 (243), 698-704
 *      D.J. Scheeres (2012), "Orbital Motion in Strongly Perturbed Environments: Applications to Asteroid, Comet and
 *          Planetary Satellite Orbiters", Springer-Praxis.
 */

#ifndef TUDAT_POLYHEDRONFUNTIONS_H
#define TUDAT_POLYHEDRONFUNTIONS_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>

#include <iostream>

namespace tudat
{
namespace basic_astrodynamics
{

/*! Computes the surface area of a polyhedron.
 *
 * Computes the surface area of a polyhedron, according to Dobrovolskis (1996), section 2.
 *
 * @param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex, 3 columns).
 * @param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet, 3 columns).
 * @return Volume.
 */
double computePolyhedronSurfaceArea (const Eigen::MatrixXd& verticesCoordinates,
                                     const Eigen::MatrixXi& verticesDefiningEachFacet );

/*! Computes the volume of a polyhedron.
 *
 * Computes the volume of a polyhedron, according to Dobrovolskis (1996), section 3.
 *
 * @param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex, 3 columns).
 * @param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet, 3 columns).
 * @return Volume.
 */
double computePolyhedronVolume (const Eigen::MatrixXd& verticesCoordinates,
                                const Eigen::MatrixXi& verticesDefiningEachFacet );

/*! Computes the mean radius of the polyhedron.
 *
 * Computes the mean radius of the polyhedron, according to Scheeres (2012), section 2.2.
 *
 * @param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex, 3 columns).
 * @param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet, 3 columns).
 * @return Mean radius.
 */
double computePolyhedronMeanRadius( const Eigen::MatrixXd& verticesCoordinates,
                                    const Eigen::MatrixXi& verticesDefiningEachFacet );

/*! Computes the centroid a polyhedron.
 *
 * Computes the volume of a polyhedron, according to Dobrovolskis (1996), section 4. When using the polyhedron as a
 * constant density gravity model, the centroid coincides with the center of mass.
 *
 * @param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex, 3 columns).
 * @param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet, 3 columns).
 * @return Position of centroid.
 */
Eigen::Vector3d computePolyhedronCentroidPosition (const Eigen::MatrixXd& verticesCoordinates,
                                                   const Eigen::MatrixXi& verticesDefiningEachFacet );

/*! Adjusts the centroid of the polyhedron to the desired value.
 *
 * Adjusts the centroid of the polyhedron to the desired value. When using the polyhedron as a constant density gravity
 * model, the centroid coincides with the center of mass; this function might be useful e.g. to ensure that the center
 * of mass coincides with the origin of some body fixed frame.
 *
 * @param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex, 3 columns).
 * @param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet, 3 columns).
 * @param desiredCentroid Desired position of the centroid wrt the frame in which the vertices were defined.
 * @return Corrected coordinates of vertices.
 */
Eigen::MatrixXd modifyPolyhedronCentroidPosition (const Eigen::MatrixXd& verticesCoordinates,
                                                  const Eigen::MatrixXi& verticesDefiningEachFacet,
                                                  const Eigen::Vector3d desiredCentroid );


/*! Computes the inertia tensor a polyhedron.
 *
 * Computes the inertia tensor of a polyhedron, according to Dobrovolskis (1996), section 5.
 *
 * @param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex, 3 columns).
 * @param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet, 3 columns).
 * @param density Density of polyhedron
 * @return Inertia tensor.
 */
Eigen::Matrix3d computePolyhedronInertiaTensor (const Eigen::MatrixXd& verticesCoordinates,
                                                const Eigen::MatrixXi& verticesDefiningEachFacet,
                                                const double density );

/*! Computes the inertia tensor a polyhedron.
 *
 * Computes the inertia tensor of a polyhedron, according to Dobrovolskis (1996), section 5.
 *
 * @param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex, 3 columns).
 * @param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet, 3 columns).
 * @param gravitationalParameter Gravitational parameter of the polyhedron.
 * @param gravitationalConstant Gravitational constant
 * @return Inertia tensor.
 */
Eigen::Matrix3d computePolyhedronInertiaTensor (const Eigen::MatrixXd& verticesCoordinates,
                                                const Eigen::MatrixXi& verticesDefiningEachFacet,
                                                const double gravitationalParameter,
                                                const double gravitationalConstant );

} // namespace basic_astrodynamics
} // namespace tudat

#endif //TUDAT_POLYHEDRONFUNTIONS_H
