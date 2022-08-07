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

double computeVolume ( const Eigen::MatrixXd& verticesCoordinates,
                      const Eigen::MatrixXi& verticesDefiningEachFacet )
{
    // Check if inputs are valid
    checkValidityOfPolyhedronSettings ( verticesCoordinates, verticesDefiningEachFacet );

    const unsigned int numberOfFacets = verticesDefiningEachFacet.rows();

    // Initialize volume
    double volume = 0;

    // Loop of tetrahedra, computing the volume of each one
    for (unsigned int facet = 0; facet < numberOfFacets; ++facet)
    {
        Eigen::Vector3d vertex0 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,0),0);
        Eigen::Vector3d vertex1 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,1),0);
        Eigen::Vector3d vertex2 = verticesCoordinates.block<1,3>(verticesDefiningEachFacet(facet,2),0);

        volume += 1.0/6.0 * vertex0.dot( vertex1.cross( vertex2 ) );
    }

    return volume;
}


Eigen::Vector3d computeCentroid ( const Eigen::MatrixXd& verticesCoordinates,
                                  const Eigen::MatrixXi& verticesDefiningEachFacet )
{
    // Check if inputs are valid
    checkValidityOfPolyhedronSettings ( verticesCoordinates, verticesDefiningEachFacet );

    const unsigned int numberOfFacets = verticesDefiningEachFacet.rows();

    // Initialize centroid
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();

    // Loop of tetrahedra, and compute the centroid of the full polyhedron
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
    Eigen::Vector3d initialCentroid = computeCentroid( verticesCoordinates, verticesDefiningEachFacet );
    Eigen::Vector3d centroidCorrection = desiredCentroid - initialCentroid;

    const unsigned int numberOfVertices = verticesCoordinates.rows();

    for (unsigned int vertex = 0; vertex < numberOfVertices; ++vertex)
    {
        verticesCoordinates.block<1,3>(vertex,0) += centroidCorrection;
    }

    return verticesCoordinates;
}