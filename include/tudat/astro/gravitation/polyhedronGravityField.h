/*    Copyright (c) 2010-2019, Delft University of Technology
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
 */

#ifndef TUDATBUNDLE_POLYHEDRONGRAVITYFIELD_H
#define TUDATBUNDLE_POLYHEDRONGRAVITYFIELD_H

#include <memory>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <iostream>

#include "tudat/astro/gravitation/gravityFieldModel.h"

namespace tudat
{

namespace gravitation
{

// Function to calculate the position of the vertices of a polyhedron relative to a field point
void calculatePolyhedronVerticesCoordinatesRelativeToFieldPoint (
        Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::Vector3d& bodyFixedPosition,
        const Eigen::MatrixXd& verticesCoordinates);

// Function to calculate the per-facet factor of a polyhedron
void calculatePolyhedronPerFacetFactor (
        Eigen::VectorXd& perFacetFactor,
        const Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::MatrixXi& verticesDefiningEachFacet);

// Function to calculate the per-edge factor of a polyhedron
void calculatePolyhedronPerEdgeFactor (
        Eigen::VectorXd& perEdgeFactor,
        const Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::MatrixXi& verticesDefiningEachEdge);

// Function to calculate the gravitational potential of a constant-density polyhedron
double calculatePolyhedronGravitationalPotential(
        const double gravitationalConstantTimesDensity,
        const Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::MatrixXi& verticesDefiningEachFacet,
        const Eigen::MatrixXi& verticesDefiningEachEdge,
        const std::vector< Eigen::MatrixXd >& facetDyads,
        const std::vector< Eigen::MatrixXd >& edgeDyads,
        const Eigen::VectorXd& perFacetFactor,
        const Eigen::VectorXd& perEdgeFactor );

// Function to calculate the gradient of the potential of a constant-density polyhedron
Eigen::Vector3d calculatePolyhedronGradientOfGravitationalPotential(
        const double gravitationalConstantTimesDensity,
        const Eigen::MatrixXd& verticesCoordinatesRelativeToFieldPoint,
        const Eigen::MatrixXi& verticesDefiningEachFacet,
        const Eigen::MatrixXi& verticesDefiningEachEdge,
        const std::vector< Eigen::MatrixXd >& facetDyads,
        const std::vector< Eigen::MatrixXd >& edgeDyads,
        const Eigen::VectorXd& perFacetFactor,
        const Eigen::VectorXd& perEdgeFactor);

// Function to calculate the Hessian matrix of the potential of a constant-density polyhedron
Eigen::Matrix3d calculatePolyhedronHessianOfGravitationalPotential(
        const double gravitationalConstantTimesDensity,
        const std::vector< Eigen::MatrixXd >& facetDyads,
        const std::vector< Eigen::MatrixXd >& edgeDyads,
        const Eigen::VectorXd& perFacetFactor,
        const Eigen::VectorXd& perEdgeFactor);

// Function to calculate the laplacian of the gravitational potential of a constant-density polyhedron
double calculatePolyhedronLaplacianOfGravitationalPotential(
        const double gravitationalConstantTimesDensity,
        const Eigen::VectorXd& perFacetFactor);


//! Cache object in which variables that are required for the computation of polyhedron gravity field are stored.
class PolyhedronGravityCache
{
public:
    PolyhedronGravityCache(
            const Eigen::MatrixXd& verticesCoordinates,
            const Eigen::MatrixXi& verticesDefiningEachFacet,
            const Eigen::MatrixXi& verticesDefiningEachEdge)
            : verticesCoordinates_( verticesCoordinates ),
              verticesDefiningEachFacet_( verticesDefiningEachFacet ),
              verticesDefiningEachEdge_( verticesDefiningEachEdge )
    {
        currentBodyFixedPosition_ = (Eigen::Vector3d() << TUDAT_NAN, TUDAT_NAN, TUDAT_NAN).finished();
    }

    void update( const Eigen::Vector3d& currentBodyFixedPosition );

    Eigen::MatrixXd& getVerticesCoordinatesRelativeToFieldPoint ( )
    { return verticesCoordinatesRelativeToFieldPoint_; }

    Eigen::VectorXd& getPerFacetFactor ( )
    { return perFacetFactor_; }

    Eigen::VectorXd& getPerEdgeFactor ( )
    { return perEdgeFactor_; }


protected:

private:

    Eigen::Vector3d currentBodyFixedPosition_;

    const Eigen::MatrixXd verticesCoordinates_;

    const Eigen::MatrixXi verticesDefiningEachFacet_;

    const Eigen::MatrixXi verticesDefiningEachEdge_;

    Eigen::MatrixXd verticesCoordinatesRelativeToFieldPoint_;

    Eigen::VectorXd perFacetFactor_;

    Eigen::VectorXd perEdgeFactor_;
};


//! Class to represent the gravity field of a constant density polyhedron
class PolyhedronGravityField: public GravityFieldModel
{
public:

    PolyhedronGravityField(
            const double gravitationalParameter,
            const double volume,
            const Eigen::MatrixXd& verticesCoordinates,
            const Eigen::MatrixXi& verticesDefiningEachFacet,
            const Eigen::MatrixXi& verticesDefiningEachEdge,
            std::vector< Eigen::MatrixXd >& facetDyads,
            std::vector< Eigen::MatrixXd >& edgeDyads,
            const std::string& fixedReferenceFrame = "",
            const std::function< void( ) > updateInertiaTensor = std::function< void( ) > ( ) )
        : GravityFieldModel(gravitationalParameter, updateInertiaTensor),
        gravitationalParameter_( gravitationalParameter ),
        volume_( volume ),
        verticesCoordinates_( verticesCoordinates ),
        verticesDefiningEachFacet_( verticesDefiningEachFacet ),
        verticesDefiningEachEdge_( verticesDefiningEachEdge ),
        facetDyads_( facetDyads ),
        edgeDyads_( edgeDyads ),
        fixedReferenceFrame_( fixedReferenceFrame )
    {
        polyhedronGravityCache_ = std::make_shared< PolyhedronGravityCache >(
                verticesCoordinates_, verticesDefiningEachFacet_, verticesDefiningEachEdge_);
    }

    virtual double getGravitationalPotential( const Eigen::Vector3d& bodyFixedPosition )
    {
        polyhedronGravityCache_->update(bodyFixedPosition);

        return calculatePolyhedronGravitationalPotential(
                gravitationalParameter_ / volume_,
                polyhedronGravityCache_->getVerticesCoordinatesRelativeToFieldPoint(),
                verticesDefiningEachFacet_,
                verticesDefiningEachEdge_,
                facetDyads_,
                edgeDyads_,
                polyhedronGravityCache_->getPerFacetFactor(),
                polyhedronGravityCache_->getPerEdgeFactor() );
    }

    virtual Eigen::Vector3d getGradientOfPotential( const Eigen::Vector3d& bodyFixedPosition )
    {
        polyhedronGravityCache_->update(bodyFixedPosition);

        return calculatePolyhedronGradientOfGravitationalPotential(
                gravitationalParameter_ / volume_,
                polyhedronGravityCache_->getVerticesCoordinatesRelativeToFieldPoint(),
                verticesDefiningEachFacet_,
                verticesDefiningEachEdge_,
                facetDyads_,
                edgeDyads_,
                polyhedronGravityCache_->getPerFacetFactor(),
                polyhedronGravityCache_->getPerEdgeFactor() );
    }

    Eigen::Matrix3d getHessianOfPotential( const Eigen::Vector3d& bodyFixedPosition )
    {
        polyhedronGravityCache_->update(bodyFixedPosition);

        return calculatePolyhedronHessianOfGravitationalPotential(
                gravitationalParameter_ / volume_,
                facetDyads_,
                edgeDyads_,
                polyhedronGravityCache_->getPerFacetFactor(),
                polyhedronGravityCache_->getPerEdgeFactor() );
    }

    virtual double getLaplacianOfPotential( const Eigen::Vector3d& bodyFixedPosition )
    {
        polyhedronGravityCache_->update(bodyFixedPosition);

        return calculatePolyhedronLaplacianOfGravitationalPotential(
                gravitationalParameter_ / volume_,
                polyhedronGravityCache_->getPerFacetFactor() );
    }

    //! Function to retrieve the identifier for body-fixed reference frame
    std::string getFixedReferenceFrame( )
    {
        return fixedReferenceFrame_;
    }

    //! Function to retrieve the volume
    double getVolume( )
    {
        return volume_;
    }

    const Eigen::MatrixXd& getVerticesCoordinates( )
    {
        return verticesCoordinates_;
    }

    const Eigen::MatrixXi& getVerticesDefiningEachFacet( )
    {
        return verticesDefiningEachFacet_;
    }

    const Eigen::MatrixXi& getVerticesDefiningEachEdge( )
    {
        return verticesDefiningEachEdge_;
    }

    const std::vector< Eigen::MatrixXd >& getFacetDyads( )
    {
        return facetDyads_;
    }

    const std::vector< Eigen::MatrixXd >& getEdgeDyads( )
    {
        return edgeDyads_;
    }

protected:

private:

    const double gravitationalParameter_;

    const double volume_;

    const Eigen::MatrixXd verticesCoordinates_;

    const Eigen::MatrixXi verticesDefiningEachFacet_;

    const Eigen::MatrixXi verticesDefiningEachEdge_;

    std::vector< Eigen::MatrixXd > facetDyads_;

    std::vector< Eigen::MatrixXd > edgeDyads_;

    std::shared_ptr< PolyhedronGravityCache > polyhedronGravityCache_;

    //! Identifier for body-fixed reference frame
    std::string fixedReferenceFrame_;

};

} // namespace gravitation

} // namespace tudat

#endif //TUDATBUNDLE_POLYHEDRONGRAVITYFIELD_H
