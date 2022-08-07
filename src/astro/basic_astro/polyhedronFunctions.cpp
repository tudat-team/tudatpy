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

#include "tudat/astro/basic_astro/polyhedronFuntions.h"

namespace tudat
{
namespace polyhedron_utilities
{

void checkValidityOfPolyhedronSettings( const Eigen::MatrixXd& verticesCoordinates,
                                        const Eigen::MatrixXi& verticesDefiningEachFacet )
{
    const unsigned int numberOfFacets = verticesDefiningEachFacet.rows();
    const unsigned int numberOfVertices = verticesCoordinates.rows();

    if ( numberOfFacets != 2 * ( numberOfVertices - 2) )
    {
        throw std::runtime_error( "Error when processing polyhedron: number of facets (" + std::to_string( numberOfFacets ) +
            ") and vertices (" + std::to_string( numberOfVertices ) + ") are not consistent." );
    }
    else if ( verticesCoordinates.cols() != 3 )
    {
        throw std::runtime_error( "Error when processing polyhedron: table with vertices coordinates has invalid "
                                  "number of columns (" + std::to_string( verticesCoordinates.cols() ) + ")." );
    }
    else if ( verticesDefiningEachFacet.cols() != 3 )
    {
        throw std::runtime_error( "Error when processing polyhedron: table with vertices defining each facet has invalid "
                                  "number of columns (" + std::to_string( verticesCoordinates.cols() ) + ")." );
    }

}

double computeSurfaceArea ( const Eigen::MatrixXd& verticesCoordinates,
                            const Eigen::MatrixXi& verticesDefiningEachFacet )
{
    // Check if inputs are valid
    checkValidityOfPolyhedronSettings ( verticesCoordinates, verticesDefiningEachFacet );

    const unsigned int numberOfFacets = verticesDefiningEachFacet.rows();

    // Initialize area
    double area = 0;

    // Loop over facets, computing the area of each one
    for (unsigned int facet = 0; facet < numberOfFacets; ++facet)
    {
        Eigen::Vector3d vertex0 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,0),0);
        Eigen::Vector3d vertex1 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,1),0);
        Eigen::Vector3d vertex2 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,2),0);

        Eigen::Vector3d facetNormal = ( vertex1 - vertex0 ).cross( vertex2 - vertex0 );

        area += facetNormal.norm() / 2.0;
    }

    return area;
}

double computeVolume ( const Eigen::MatrixXd& verticesCoordinates,
                      const Eigen::MatrixXi& verticesDefiningEachFacet )
{
    // Check if inputs are valid
    checkValidityOfPolyhedronSettings ( verticesCoordinates, verticesDefiningEachFacet );

    const unsigned int numberOfFacets = verticesDefiningEachFacet.rows();

    // Initialize volume
    double volume = 0;

    // Loop over tetrahedra, computing the volume of each one
    for (unsigned int facet = 0; facet < numberOfFacets; ++facet)
    {
        Eigen::Vector3d vertex0 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,0),0);
        Eigen::Vector3d vertex1 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,1),0);
        Eigen::Vector3d vertex2 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,2),0);

        volume += 1.0/6.0 * vertex0.dot( vertex1.cross( vertex2 ) );
    }

    return volume;
}

Eigen::Vector3d computeCentroidPosition ( const Eigen::MatrixXd& verticesCoordinates,
                                          const Eigen::MatrixXi& verticesDefiningEachFacet )
{
    // Check if inputs are valid
    checkValidityOfPolyhedronSettings ( verticesCoordinates, verticesDefiningEachFacet );

    const unsigned int numberOfFacets = verticesDefiningEachFacet.rows();

    // Initialize centroid
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();

    // Loop over tetrahedra, and compute the centroid of the full polyhedron
    for (unsigned int facet = 0; facet < numberOfFacets; ++facet)
    {
        Eigen::Vector3d vertex0 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,0),0);
        Eigen::Vector3d vertex1 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,1),0);
        Eigen::Vector3d vertex2 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,2),0);

        double volumeOfTetrahedron = 1.0/6.0 * vertex0.dot( vertex1.cross( vertex2 ) );
        Eigen::Vector3d centerOfMassOfTetrahedron = (vertex0 + vertex1 + vertex2) / 4.0;
        centroid += volumeOfTetrahedron * centerOfMassOfTetrahedron;
    }

    centroid = centroid / computeVolume( verticesCoordinates, verticesDefiningEachFacet );

    return centroid;
}

Eigen::MatrixXd modifyCentroidPosition ( Eigen::MatrixXd verticesCoordinates,
                                         const Eigen::MatrixXi& verticesDefiningEachFacet,
                                         const Eigen::Vector3d desiredCentroid )
{
    Eigen::Vector3d initialCentroid = computeCentroidPosition( verticesCoordinates, verticesDefiningEachFacet );
    Eigen::Vector3d centroidCorrection = desiredCentroid - initialCentroid;

    // Correct the centroid position if necessary
    if ( desiredCentroid != initialCentroid )
    {
        const unsigned int numberOfVertices = verticesCoordinates.rows();

        for (unsigned int vertex = 0; vertex < numberOfVertices; ++vertex)
        {
            verticesCoordinates.block<1,3>(vertex,0) += centroidCorrection;
        }
    }

    return verticesCoordinates;
}

Eigen::Matrix3d computeInertiaTensor ( const Eigen::MatrixXd& verticesCoordinates,
                                       const Eigen::MatrixXi& verticesDefiningEachFacet,
                                       const double density )
{
    // Check if inputs are valid
    checkValidityOfPolyhedronSettings ( verticesCoordinates, verticesDefiningEachFacet );

    const unsigned int numberOfFacets = verticesDefiningEachFacet.rows();

    // Initialize matrix with products of inertia. Only values on and above the diagonal will be filled
    Eigen::Matrix3d inertiaProducts = Eigen::Matrix3d::Zero();

    // Loop over tetrahedra, and compute the centroid of the full polyhedron
    for (unsigned int facet = 0; facet < numberOfFacets; ++facet)
    {
        Eigen::Vector3d vertex0 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,0),0);
        Eigen::Vector3d vertex1 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,1),0);
        Eigen::Vector3d vertex2 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,2),0);

        Eigen::Vector3d edge0 = vertex1 - vertex0;
        Eigen::Vector3d edge1 = vertex2 - vertex0;

        double volumeOfTetrahedron = 1.0/6.0 * vertex0.dot( vertex1.cross( vertex2 ) );

        for (unsigned int row = 0; row < 3; ++row)
        {
            for (unsigned int col = row; col <= 3; ++col)
            {

                inertiaProducts(row, col) += density * volumeOfTetrahedron / 20.0 * (
                        12 * vertex0(row) * vertex0(col) + 2 * edge0(row) * edge0(col) + 2 * edge1(row) * edge1(col) +
                        4 * ( vertex0(row) * edge0(col) + vertex0(col) * edge0(row) ) +
                        4 * ( vertex0(row) * edge1(col) + vertex0(col) * edge1(row) ) +
                        edge0(col) * edge1(row) + edge0(row) * edge1(col) );
            }
        }
    }

     // Compute inertia tensor
    Eigen::Matrix3d inertiaTensor;
    inertiaTensor(0,0) = inertiaProducts(1,1) + inertiaProducts(2,2);
    inertiaTensor(1,1) = inertiaProducts(0,0) + inertiaProducts(2,2);
    inertiaTensor(2,2) = inertiaProducts(0,0) + inertiaProducts(1,1);
    inertiaTensor(0,1) = - inertiaProducts(0,1);
    inertiaTensor(1,0) = inertiaTensor(0,1);
    inertiaTensor(0,2) = - inertiaProducts(0,2);
    inertiaTensor(2,0) = inertiaTensor(0,2);
    inertiaTensor(1,2) = - inertiaProducts(1,2);
    inertiaTensor(2,1) = inertiaTensor(1,2);

    return inertiaTensor;
}

Eigen::Matrix3d computeInertiaTensor ( const Eigen::MatrixXd& verticesCoordinates,
                                       const Eigen::MatrixXi& verticesDefiningEachFacet,
                                       const double gravitationalParameter,
                                       const double gravitationalConstant )
{
    double density = gravitationalParameter / ( computeVolume(verticesCoordinates, verticesDefiningEachFacet) *
            gravitationalConstant );

    return computeInertiaTensor(verticesCoordinates, verticesDefiningEachFacet, density);
}

} // namespace polyhedron_utilities
} // namespace tudat