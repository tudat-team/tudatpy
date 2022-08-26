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
 *       Avillez (2022), MSc thesis (TU Delft) - TODO: add proper reference
 */

#ifndef TUDATBUNDLE_POLYHEDRONBODYSHAPEMODEL_H
#define TUDATBUNDLE_POLYHEDRONBODYSHAPEMODEL_H

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/bodyShapeModel.h"
#include "tudat/astro/basic_astro/polyhedronFuntions.h"
#include "tudat/astro/gravitation/polyhedronGravityField.h"
#include <iostream>

namespace tudat
{
namespace basic_astrodynamics
{

// Body shape model defined by a polyhedron, usually used for approximating small bodies (e.g. small moons and asteroids)
class PolyhedronBodyShapeModel: public BodyShapeModel
{
public:

    /*! Constructor.
     *
     * Constructor.
     * @param verticesCoordinates Matrix with coordinates of the polyhedron vertices. Each row represents the (x,y,z)
     * coordinates of one vertex.
     * @param verticesDefiningEachFacet Matrix with the indices (0 indexed) of the vertices defining each facet. Each
     * row contains 3 indices, which must be provided in counterclockwise order when seen from outside the polyhedron.
     * @param computeAltitudeWithSign Flag indicating whether the altitude should be computed with sign (i.e. >0 if
     * above surface, <0 otherwise) or having always a positive value.
     * @param justComputeDistanceToVertices Flag indicating whether the distance should be computed wrt to all the
     * polyhedron features or wrt to just the vertices.
     */
    PolyhedronBodyShapeModel(
            const Eigen::MatrixXd& verticesCoordinates,
            const Eigen::MatrixXi& verticesDefiningEachFacet,
            const bool computeAltitudeWithSign,
            const bool justComputeDistanceToVertices ):
        verticesCoordinates_( verticesCoordinates ),
        verticesDefiningEachFacet_( verticesDefiningEachFacet ),
        computeAltitudeWithSign_( computeAltitudeWithSign ),
        justComputeDistanceToVertices_( justComputeDistanceToVertices ),
        averageRadius_( TUDAT_NAN )
    {
        // Check if provided settings are valid
        polyhedron_utilities::checkValidityOfPolyhedronSettings( verticesCoordinates, verticesDefiningEachFacet );

        // If necessary, get list with vertices defining each edge
        if ( !justComputeDistanceToVertices_ )
        {
            computeVerticesDefiningEachEdge();
        }
    }

    //! Destructor
    ~PolyhedronBodyShapeModel( ){ }

    //! Calculates the altitude above the polyhedron
    /*!
     *  Function to calculate the altitude above the polyhedron from a body fixed position.
     *  Function computes the minimum distance to each of the polyhedron features (vertices, edges and facets); the
     *  distance is only computed wrt to the edges and facets around the closest vertex. See Avillez (2022).
     *  \param bodyFixedPosition Cartesian, body-fixed position of the point at which the altitude
     *  is to be determined.
     *  \return Altitude above the oblate spheroid.
     */
    double getAltitude( const Eigen::Vector3d& bodyFixedPosition );

    /*! Function to set the average radius of the polyhedron.
     *
     * Function to set the average radius of the polyhedron. This function is necessary as currently no function to
     * compute the average radius is implemented.
     * @param averageRadius Average radius of the polyhedron.
     */
    void setAverageRadius( const double averageRadius )
    {
        averageRadius_ = averageRadius;
    }

    //! Function to return the average radius of the polyhedron.
    /*!
     *  Function to return the average radius of the polyhedron.
     *  \return Average radius of the polyhedron.
     */
    double getAverageRadius( )
    {
        if ( std::isnan( averageRadius_ ) )
        {
            throw std::runtime_error( "The average radius of the polyhedron shape model was not defined." );
        }

        std::cerr << "Warning: the returned polyhedron average radius was manually defined and does not corresponds to"
                     "the true average radius." << std::endl;

        return averageRadius_;
    }

    // Function to return the vertices coordinates.
    const Eigen::MatrixXd& getVerticesCoordinates( )
    {
        return verticesCoordinates_;
    }

    // Function to return the vertices defining each facet.
    const Eigen::MatrixXi& getVerticesDefiningEachFacet( )
    {
        return verticesDefiningEachFacet_;
    }

    // Function to return the computeAltitudeWithSign flag.
    bool getComputeAltitudeWithSign( )
    {
        return computeAltitudeWithSign_;
    }

    // Function to return the justComputeDistanceToVertices flag.
    bool getJustComputeDistanceToVertices( )
    {
        return justComputeDistanceToVertices_;
    }

private:

    /*! Computes the distance to the vertex closest to the field point.
     *
     * Computes the distance to the vertex closest to the field point.
     * @param bodyFixedPosition Cartesian, body-fixed position of the point at which the altitude
     *  is to be determined.
     * @return Distance to the vertex closest to the field point.
     */
    double computeDistanceToClosestVertex( const Eigen::Vector3d& bodyFixedPosition,
                                           unsigned int& closestVertexId);

    /*! Computes the distance to the facet closest to the field point.
     *
     * Computes the distance to the facet closest to the field point, according to Avillez (2022). The distance is
     * computed only wrt to the limited set of facets provided as argument. If no valid distance is computed, the
     * function returns NAN. The returned distance is unsigned.
     * @param bodyFixedPosition Cartesian, body-fixed position of the point at which the altitude
     *  is to be determined.
     * @param verticesDefiningEachFacetToEvaluate Facets wrt which the distance is to be computed.
     * @return Distance to closest facet.
     */
    double computeDistanceToClosestFacet ( const Eigen::Vector3d& bodyFixedPosition,
                                           const Eigen::MatrixXi& verticesDefiningEachFacetToEvaluate );

    /*! Computes the distance to the edge closest to the field point.
     *
     * Computes the distance to the edge closest to the field point, according to Avillez (2022). The distance is
     * computed only wrt to the limited set of edges provided as argument. If no valid distance is computed, the
     * function returns NAN. The returned distance is unsigned.
     * @param bodyFixedPosition Cartesian, body-fixed position of the point at which the altitude
     *  is to be determined.
     * @param verticesDefiningEachEdgeToEvaluate Edges wrt which the distance is to be computed.
     * @return Distance to closest edge.
     */
    double computeDistanceToClosestEdge ( const Eigen::Vector3d& bodyFixedPosition,
                                          const Eigen::MatrixXi& verticesDefiningEachEdgeToEvaluate );

    /*! Computes the matrix with the indices of the vertices defining each edge.
     *
     * Computes the matrix with the indices of the vertices defining each edge; saves the list to verticesDefiningEachEdge_.
     */
    void computeVerticesDefiningEachEdge( );


    // Matrix with coordinates of the polyhedron vertices.
    Eigen::MatrixXd verticesCoordinates_;

    // Matrix with the indices (0 indexed) of the vertices defining each facet.
    Eigen::MatrixXi verticesDefiningEachFacet_;

    // Matrix with the indices (0 indexed) of the vertices defining each edge.
    Eigen::MatrixXi verticesDefiningEachEdge_;

    // Flag indicating whether the altitude should be computed with sign (i.e. >0 if above surface, <0 otherwise) or
    // having always a positive value
    bool computeAltitudeWithSign_;

    // Flag indicating whether the distance should be computed wrt to all the polyhedron features or wrt to just the
    // vertices.
    bool justComputeDistanceToVertices_;

    // Average radius of the polyhedron
    double averageRadius_;

};

} // namespace basic_astrodynamics
} // namespace tudat

#endif //TUDATBUNDLE_POLYHEDRONBODYSHAPEMODEL_H
