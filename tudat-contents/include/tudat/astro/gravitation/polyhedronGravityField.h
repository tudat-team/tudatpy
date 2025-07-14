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

#ifndef TUDAT_POLYHEDRONGRAVITYFIELD_H
#define TUDAT_POLYHEDRONGRAVITYFIELD_H

#include <memory>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <iostream>

#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/polyhedronFuntions.h"
#include "tudat/math/basic/polyhedron.h"

namespace tudat
{

namespace gravitation
{

//! Cache object in which variables that are required for the computation of polyhedron gravity field are stored.
class PolyhedronGravityCache
{
public:
    /*! Constructor.
     *
     * Constructor.
     * @param verticesCoordinates Matrix with coordinates of the polyhedron vertices. Each row represents the (x,y,z)
     * coordinates of one vertex.
     * @param verticesDefiningEachFacet Matrix with the indices (0 indexed) of the vertices defining each facet. Each
     * row contains 3 indices, which must be provided in counterclockwise order when seen from outise the polyhedron.
     * @param verticesDefiningEachEdge Matrix with the indices (0 indexed) of the vertices defining each facet. Each
     * row contains 2 indices.
     */
    PolyhedronGravityCache( const Eigen::MatrixXd& verticesCoordinates,
                            const Eigen::MatrixXi& verticesDefiningEachFacet,
                            const Eigen::MatrixXi& verticesDefiningEachEdge ):
        verticesCoordinates_( verticesCoordinates ), verticesDefiningEachFacet_( verticesDefiningEachFacet ),
        verticesDefiningEachEdge_( verticesDefiningEachEdge )
    {
        currentBodyFixedPosition_ = ( Eigen::Vector3d( ) << TUDAT_NAN, TUDAT_NAN, TUDAT_NAN ).finished( );
    }

    /*! Update cached variables to current state.
     *
     * Update cached variables to current state.
     * @param currentBodyFixedPosition Current body fixed position.
     */
    void update( const Eigen::Vector3d& currentBodyFixedPosition );

    /*! Function to retrieve the coordinates of the polyhedron vertices wrt field point.
     *
     * Function to retrieve the coordinates of the polyhedron vertices wrt field point.
     * @return Coordinates of the polyhedron vertices wrt field point.
     */
    Eigen::MatrixXd& getVerticesCoordinatesRelativeToFieldPoint( )
    {
        return currentVerticesCoordinatesRelativeToFieldPoint_;
    }

    /*! Function to retrieve the vector of per-facet factors.
     *
     * Function to retrieve the vector of per-facet factors.
     * @return Per-facet factors.
     */
    Eigen::VectorXd& getPerFacetFactor( )
    {
        return currentPerFacetFactor_;
    }

    /*! Function to retrieve the vector of per-edge factors.
     *
     * Function to retrieve the vector of per-edge factors.
     * @return Per-edge factors.
     */
    Eigen::VectorXd& getPerEdgeFactor( )
    {
        return currentPerEdgeFactor_;
    }

protected:
private:
    // Current body fixed position.
    Eigen::Vector3d currentBodyFixedPosition_;

    // Matrix with coordinates of the polyhedron vertices.
    const Eigen::MatrixXd verticesCoordinates_;

    // Matrix with the indices (0 indexed) of the vertices defining each facet.
    const Eigen::MatrixXi verticesDefiningEachFacet_;

    // Matrix with the indices (0 indexed) of the vertices defining each facet.
    const Eigen::MatrixXi verticesDefiningEachEdge_;

    // Current vertices coordinates wrt body fixed position.
    Eigen::MatrixXd currentVerticesCoordinatesRelativeToFieldPoint_;

    // Current value of the per-facet factors.
    Eigen::VectorXd currentPerFacetFactor_;

    // Current value of the per-edge factors.
    Eigen::VectorXd currentPerEdgeFactor_;
};

//! Class to represent the gravity field of a constant density polyhedron.
class PolyhedronGravityField : public GravityFieldModel
{
public:
    /*! Constructor.
     *
     * Constructor.
     * @param gravitationalParameter Gravitational parameter of the polyhedron.
     * @param volume Volume of the polyhedron.
     * @param verticesCoordinates Cartesian coordinates of each vertex (one row per vertex, 3 columns).
     * @param verticesDefiningEachFacet Index (0 based) of the vertices constituting each facet (one row per facet, 3 columns).
     * @param fixedReferenceFrame Identifier for body-fixed reference frame to which the field is fixed (optional).
     * @param updateInertiaTensor Function that is to be called to update the inertia tensor (typicaly in Body class;
     * default empty)
     */
    PolyhedronGravityField( const double gravitationalParameter,
                            const Eigen::MatrixXd& verticesCoordinates,
                            const Eigen::MatrixXi& verticesDefiningEachFacet,
                            const std::string& fixedReferenceFrame = "",
                            const std::function< void( ) > updateInertiaTensor = std::function< void( ) >( ) ):
        GravityFieldModel( gravitationalParameter, updateInertiaTensor ), gravitationalParameter_( gravitationalParameter ),
        verticesCoordinates_( verticesCoordinates ), verticesDefiningEachFacet_( verticesDefiningEachFacet ),
        fixedReferenceFrame_( fixedReferenceFrame )
    {
        // Check if provided arguments are valid
        basic_mathematics::checkValidityOfPolyhedronSettings( verticesCoordinates, verticesDefiningEachFacet );

        // Compute volume
        volume_ = basic_astrodynamics::computePolyhedronVolume( verticesCoordinates, verticesDefiningEachFacet );
        density_ = gravitationalParameter_ / physical_constants::GRAVITATIONAL_CONSTANT / volume_;

        // Compute edges in polyhedron
        computeVerticesAndFacetsDefiningEachEdge( );

        // Compute facet dyads
        computeFacetNormalsAndDyads( );

        // Compute edge dyads
        computeEdgeDyads( );

        // Create cache object
        polyhedronGravityCache_ =
                std::make_shared< PolyhedronGravityCache >( verticesCoordinates_, verticesDefiningEachFacet_, verticesDefiningEachEdge_ );

        inertiaTensor_ = basic_astrodynamics::computePolyhedronInertiaTensor( verticesCoordinates_, verticesDefiningEachFacet_, density_ );
    }

    /*! Function to calculate the gravitational potential.
     *
     * Function to calculate the gravitational potential.
     * @param bodyFixedPosition Position of point at which potential is to be calculated, in body-fixed frame.
     * @return Gravitational potential.
     */
    virtual double getGravitationalPotential( const Eigen::Vector3d& bodyFixedPosition )
    {
        polyhedronGravityCache_->update( bodyFixedPosition );

        return basic_mathematics::calculatePolyhedronGravitationalPotential(
                gravitationalParameter_ / volume_,
                polyhedronGravityCache_->getVerticesCoordinatesRelativeToFieldPoint( ),
                verticesDefiningEachFacet_,
                verticesDefiningEachEdge_,
                facetDyads_,
                edgeDyads_,
                polyhedronGravityCache_->getPerFacetFactor( ),
                polyhedronGravityCache_->getPerEdgeFactor( ) );
    }

    /*! Function to calculate the gradient of the gravitational potential (i.e. the acceleration).
     *
     * Function to calculate the gradient of the gravitational potential (i.e. the acceleration).
     * @param bodyFixedPosition Position of point at which potential is to be calculated, in body-fixed frame.
     * @return Gradient of the gravitational potential.
     */
    virtual Eigen::Vector3d getGradientOfPotential( const Eigen::Vector3d& bodyFixedPosition )
    {
        polyhedronGravityCache_->update( bodyFixedPosition );

        return basic_mathematics::calculatePolyhedronGradientOfGravitationalPotential(
                gravitationalParameter_ / volume_,
                polyhedronGravityCache_->getVerticesCoordinatesRelativeToFieldPoint( ),
                verticesDefiningEachFacet_,
                verticesDefiningEachEdge_,
                facetDyads_,
                edgeDyads_,
                polyhedronGravityCache_->getPerFacetFactor( ),
                polyhedronGravityCache_->getPerEdgeFactor( ) );
    }

    /*! Function to calculate the hessian matrix of the gravitational potential.
     *
     * Function to calculate the hessian matrix of the gravitational potential.
     * @param bodyFixedPosition Position of point at which potential is to be calculated, in body-fixed frame.
     * @return Hessian of the gravitational potential.
     */
    Eigen::Matrix3d getHessianOfPotential( const Eigen::Vector3d& bodyFixedPosition )
    {
        polyhedronGravityCache_->update( bodyFixedPosition );

        return basic_mathematics::calculatePolyhedronHessianOfGravitationalPotential( gravitationalParameter_ / volume_,
                                                                                      facetDyads_,
                                                                                      edgeDyads_,
                                                                                      polyhedronGravityCache_->getPerFacetFactor( ),
                                                                                      polyhedronGravityCache_->getPerEdgeFactor( ) );
    }

    /*! Function to calculate the laplacian of the gravitational potential.
     *
     * Function to calculate the laplacian of the gravitational potential.
     * @param bodyFixedPosition Position of point at which potential is to be calculated, in body-fixed frame.
     * @return Laplacian of the gravitational potential.
     */
    virtual double getLaplacianOfPotential( const Eigen::Vector3d& bodyFixedPosition )
    {
        polyhedronGravityCache_->update( bodyFixedPosition );

        return basic_mathematics::calculatePolyhedronLaplacianOfGravitationalPotential( gravitationalParameter_ / volume_,
                                                                                        polyhedronGravityCache_->getPerFacetFactor( ) );
    }

    //! Function to retrieve the identifier for the body-fixed reference frame.
    std::string getFixedReferenceFrame( )
    {
        return fixedReferenceFrame_;
    }

    //! Function to return the volume
    double getVolume( )
    {
        return volume_;
    }

    //! Function to return the vertices coordinates.
    const Eigen::MatrixXd& getVerticesCoordinates( )
    {
        return verticesCoordinates_;
    }

    //! Function to return the vertices defining each facet.
    const Eigen::MatrixXi& getVerticesDefiningEachFacet( )
    {
        return verticesDefiningEachFacet_;
    }

    //! Function to return the vertices defining each edge.
    const Eigen::MatrixXi& getVerticesDefiningEachEdge( )
    {
        return verticesDefiningEachEdge_;
    }

    //! Function to return the facet dyads.
    const std::vector< Eigen::MatrixXd >& getFacetDyads( )
    {
        return facetDyads_;
    }

    //! Function to return the edge dyads.
    const std::vector< Eigen::MatrixXd >& getEdgeDyads( )
    {
        return edgeDyads_;
    }

    virtual Eigen::Matrix3d getInertiaTensor( )
    {
        return inertiaTensor_;
    }

protected:
private:
    /*! Function to compute the vertices and facets defining each edge.
     *
     * Function to compute the vertices and facets defining each edge, according to Werner and Scheeres (1997).
     * VerticesDefiningEachEdge is a matrix with
     * the index (0 based) of the vertices constituting each edge (one row per edge, 2 columns). FacetsDefiningEachEdge
     * is a matrix with the index (0 based) of the facets defining each edge (one row per edge, 2 columns). The function
     * saves the two computed matrices into member variables of the class.
     */
    void computeVerticesAndFacetsDefiningEachEdge( );

    /*! Function to compute the facet normals and facet dyads.
     *
     * Function to compute the facet normals and facet dyads, according to Werner and Scheeres (1997).
     * Computes a vector of facet normals, each member of the
     * vector being the normal to one facet. Computes a vector of facet dyads, each member of the vector being the ´
     * dyad associated with one facet. The function saves the two computed vectors into member variables of the class.
     */
    void computeFacetNormalsAndDyads( );

    /*! Function to compute the edge dyads.
     *
     * Function to compute the edge dyads, according to Werner and Scheeres (1997). Computes a vector of edge dyads,
     * each member of the vector being the dyad associated with one edge. The function saves the computed vector
     * into a member variable of the class.
     */
    void computeEdgeDyads( );

    //! Gravitational parameter of the polyhedron.
    double gravitationalParameter_;

    //! Volume of the polyhedron.
    double volume_;

    double density_;

    Eigen::Matrix3d inertiaTensor_;

    //! Cartesian coordinates of each vertex (one row per vertex, 3 columns).
    Eigen::MatrixXd verticesCoordinates_;

    //! Index (0 based) of the vertices constituting each facet (one row per facet, 3 columns).
    Eigen::MatrixXi verticesDefiningEachFacet_;

    //! Index (0 based) of the vertices constituting each edge (one row per edge, 2 columns).
    Eigen::MatrixXi verticesDefiningEachEdge_;

    //! Index (0 based) of the facets defining each edge (one row per edge, 2 columns).
    Eigen::MatrixXi facetsDefiningEachEdge_;

    //! Vector with the outward-pointing normal vector of each facet.
    std::vector< Eigen::Vector3d > facetNormalVectors_;

    //! Vector with the facet dyad of each facet.
    std::vector< Eigen::MatrixXd > facetDyads_;

    //! Vector with the edge dyad of each edge.
    std::vector< Eigen::MatrixXd > edgeDyads_;

    //! Polyhedron cache.
    std::shared_ptr< PolyhedronGravityCache > polyhedronGravityCache_;

    //! Identifier for body-fixed reference frame
    std::string fixedReferenceFrame_;
};

}  // namespace gravitation

}  // namespace tudat

#endif  // TUDAT_POLYHEDRONGRAVITYFIELD_H
