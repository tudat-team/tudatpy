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

#ifndef TUDATBUNDLE_POLYHEDRONGRAVITYMODEL_H
#define TUDATBUNDLE_POLYHEDRONGRAVITYMODEL_H

#include <memory>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <iostream>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/gravitation/polyhedronGravityField.h"

namespace tudat
{
namespace gravitation
{

template< typename StateMatrix >
class PolyhedronGravitationalAccelerationModel: public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{

protected:
    //! Typedef for a position-returning function.
    typedef std::function< void( StateMatrix& ) > StateFunction;

public:

    //! Constructor taking position-functions for bodies, and constant parameters of polyhedron paramers.
    /*!
     * Constructor taking a pointer to a function returning the position of the body subject to
     * gravitational acceleration, constant gravitational parameter and volume of the
     * body exerting the acceleration, and a pointer to a function returning the position of the body exerting the
     * gravitational acceleration (typically the central body). This constructor uses the
     * Boost::lambda library to create a function on-the-fly that returns the constant
     * gravitational parameter and volume. The
     * constructor also updates all the internal members. The position of the body exerting the
     * gravitational acceleration is an optional parameter; the default position is the origin.
     * \param positionOfBodySubjectToAccelerationFunction Pointer to function returning position of
     *          body subject to gravitational acceleration.
     * \param aGravitationalParameter A (constant) gravitational parameter [m^2 s^-3].
     * \param aVolume A (constant) volume [m^3].
     * \param positionOfBodyExertingAccelerationFunction Pointer to function returning position of
     *          body exerting gravitational acceleration (default = (0,0,0)).
     * \param rotationFromBodyFixedToIntegrationFrameFunction Function providing the rotation from
     * body-fixes from to the frame in which the numerical integration is performed.
     * \param isMutualAttractionUsed Variable denoting whether attraction from body undergoing acceleration on
     * body exerting acceleration is included (i.e. whether aGravitationalParameter refers to the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     * \param polyhedronCache Cache object for computing/retrieving repeated terms in polyhedron potential
     *          gradient calculation.
     */
    PolyhedronGravitationalAccelerationModel(
            const StateFunction positionOfBodySubjectToAccelerationFunction,
            const double aGravitationalParameter,
            const double aVolume,
            const Eigen::MatrixXd& aVerticesCoordinatesMatrix,
            const Eigen::MatrixXi& aVerticesDefiningEachFacetMatrix,
            const Eigen::MatrixXi& aVerticesDefiningEachEdgeMatrix,
            std::vector< Eigen::MatrixXd >& aFacetDyadsVector,
            std::vector< Eigen::MatrixXd >& aEdgeDyadsVector,
            const StateFunction positionOfBodyExertingAccelerationFunction =
                [ ]( Eigen::Vector3d& input ){ input = Eigen::Vector3d::Zero( ); },
            const std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction =
                [ ]( ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); },
            const bool isMutualAttractionUsed = 0,
            std::shared_ptr< PolyhedronGravityCache > polyhedronCache =
                std::make_shared< PolyhedronGravityCache >( ) )
        : subjectPositionFunction_( positionOfBodySubjectToAccelerationFunction ),
        gravitationalParameterFunction_( [ = ]( ){ return aGravitationalParameter; } ),
        volumeFunction_( [ = ]( ){ return aVolume; } ),
        getVerticesCoordinates_( [ = ]( ){ return aVerticesCoordinatesMatrix; } ),
        getVerticesDefiningEachFacet_( [ = ]( ){ return aVerticesDefiningEachFacetMatrix; } ),
        getVerticesDefiningEachEdge_( [ = ]( ){ return aVerticesDefiningEachEdgeMatrix; } ),
        getFacetDyads_( [ = ]( ){ return aFacetDyadsVector; } ),
        getEdgeDyads_( [ = ]( ){ return aEdgeDyadsVector; } ),
        sourcePositionFunction_( positionOfBodyExertingAccelerationFunction ),
        rotationFromBodyFixedToIntegrationFrameFunction_( rotationFromBodyFixedToIntegrationFrameFunction ),
        isMutualAttractionUsed_( isMutualAttractionUsed ),
        polyhedronCache_( polyhedronCache )
    { }

    //! Update class members.
    /*!
     * Updates all the base class members to their current values and also updates the class
     * members of this class.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {

            rotationToIntegrationFrame_ = rotationFromBodyFixedToIntegrationFrameFunction_( );

            currentInertialRelativePosition_ =
                    this->positionOfBodySubjectToAcceleration - this->positionOfBodyExertingAcceleration ;

            currentRelativePosition_ = rotationToIntegrationFrame_.inverse( ) * currentInertialRelativePosition_;

            polyhedronCache_->update( currentRelativePosition_ );

            currentAcceleration_ = calculatePolyhedronGradientOfGravitationalPotential(
                    gravitationalParameterFunction_() / volumeFunction_(),
                    polyhedronCache_->getVerticesCoordinatesRelativeToFieldPoint(),
                    getVerticesDefiningEachFacet_(),
                    getVerticesDefiningEachEdge_(),
                    getFacetDyads_(),
                    getEdgeDyads_(),
                    polyhedronCache_->getPerFacetFactor(),
                    polyhedronCache_->getPerEdgeFactor() );

            currentAccelerationInBodyFixedFrame_ = rotationToIntegrationFrame_.inverse( ) * currentAcceleration_;
        }
    }

    //! Function to return current position vector from body exerting acceleration to body undergoing acceleration, in frame
    //! fixed to body undergoing acceleration
    Eigen::Vector3d getCurrentRelativePosition( )
    {
        return currentRelativePosition_;
    }

    //! Function to return current position vector from body exerting acceleration to body undergoing acceleration, in inertial
    //! frame
    Eigen::Vector3d getCurrentInertialRelativePosition( )
    {
        return currentInertialRelativePosition_;
    }

    //! Function to retrieve the current rotation from body-fixed frame to integration frame, in the form of a quaternion.
    Eigen::Quaterniond getCurrentRotationToIntegrationFrame( )
    {
        return rotationToIntegrationFrame_;
    }

    //! Function to retrieve the current rotation from body-fixed frame to integration frame, as a rotation matrix.
    Eigen::Matrix3d getCurrentRotationToIntegrationFrameMatrix( )
    {
        return rotationToIntegrationFrame_.toRotationMatrix( );
    }

private:

    //! Pointer to function returning position of body subject to acceleration.
    const StateFunction subjectPositionFunction_;

    //! Function returning a gravitational parameter [m^3 s^-2].
    const std::function< double( ) > gravitationalParameterFunction_;

    //! Function returning a volume [m^3].
    const std::function< double( ) > volumeFunction_;

    //! Pointer to function returning the vertices coordinates
    const std::function< Eigen::MatrixXd( ) > getVerticesCoordinates_;

    //! Pointer to function returning the vertices defining each facet
    const std::function< Eigen::MatrixXi( ) > getVerticesDefiningEachFacet_;

    //! Pointer to function returning the vertices defining each edge
    const std::function< Eigen::MatrixXi( ) > getVerticesDefiningEachEdge_;

    //! Pointer to function returning the facet dyads
    const std::function< std::vector< Eigen::MatrixXd >( ) > getFacetDyads_;

    //! Pointer to function returning the edge dyads
    const std::function< std::vector< Eigen::MatrixXd >( ) > getEdgeDyads_;

    //! Pointer to function returning position of body exerting acceleration.
    const StateFunction sourcePositionFunction_;

    //! Function returning the current rotation from body-fixed frame to integration frame.
    std::function< Eigen::Quaterniond( ) > rotationFromBodyFixedToIntegrationFrameFunction_;

    //! Variable denoting whether mutual acceleration between bodies is included.
    bool isMutualAttractionUsed_;

    //!  Polyhedron for this acceleration
    std::shared_ptr< PolyhedronGravityCache > polyhedronCache_;

    //! Current rotation from body-fixed frame to integration frame.
    Eigen::Quaterniond rotationToIntegrationFrame_;

    //! Current position vector from body exerting acceleration to body undergoing acceleration, in inertial frame
    Eigen::Vector3d currentInertialRelativePosition_;

    //! Current position vector from body exerting acceleration to body undergoing acceleration, in frame fixed to body
    //! undergoing acceleration
    Eigen::Vector3d currentRelativePosition_;

    //! Current acceleration in frame fixed to body undergoing acceleration, as computed by last call to updateMembers function
    Eigen::Vector3d currentAccelerationInBodyFixedFrame_;

};


} // namespace gravitation

} // namespace tudat

#endif //TUDATBUNDLE_POLYHEDRONGRAVITYMODEL_H
