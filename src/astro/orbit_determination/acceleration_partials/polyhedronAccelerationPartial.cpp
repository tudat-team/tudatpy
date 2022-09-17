/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/polyhedronAccelerationPartial.h"
#include "tudat/math/basic/polyhedron.h"
#include "tudat/astro/orbit_determination/acceleration_partials/centralGravityAccelerationPartial.h"

namespace tudat
{
namespace acceleration_partials
{

//! Contructor.
PolyhedronGravityPartial::PolyhedronGravityPartial (
        const std::string& acceleratedBody,
        const std::string& acceleratingBody,
        const std::shared_ptr< gravitation::PolyhedronGravitationalAccelerationModel > accelerationModel,
        const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartials ):
    AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::polyhedron_gravity ),
    gravitationalParameterFunction_( accelerationModel->getGravitationalParameterFunction( ) ),
    volumeFunction_( accelerationModel->getVolumeFunction( ) ),
    polyhedronCache_( accelerationModel->getPolyhedronCache() ),
    facetDyads_( accelerationModel->getFacetDyadsFunction( )( ) ),
    edgeDyads_( accelerationModel->getEdgeDyadsFunction( )( ) ),
    positionFunctionOfAcceleratedBody_( std::bind( &gravitation::PolyhedronGravitationalAccelerationModel::
                                                   getCurrentPositionOfBodySubjectToAcceleration, accelerationModel ) ),
    positionFunctionOfAcceleratingBody_( std::bind( &gravitation::PolyhedronGravitationalAccelerationModel::
                                                    getCurrentPositionOfBodyExertingAcceleration, accelerationModel ) ),
    fromBodyFixedToIntegrationFrameRotation_( std::bind( &gravitation::PolyhedronGravitationalAccelerationModel::
                                                         getCurrentRotationToIntegrationFrameMatrix, accelerationModel ) ),
    accelerationFunction_( std::bind( &gravitation::PolyhedronGravitationalAccelerationModel::getAcceleration,
                                      accelerationModel ) ),
    updateFunction_( std::bind( &gravitation::PolyhedronGravitationalAccelerationModel::updateMembers,
                                accelerationModel, std::placeholders::_1 ) ),
    rotationMatrixPartials_( rotationMatrixPartials )
{

}

void PolyhedronGravityPartial::update( const double currentTime )
{
    using namespace tudat::coordinate_conversions;

    if( !( currentTime_ == currentTime ) )
    {
        // Update acceleration model
        updateFunction_( currentTime );

        // Calculate Cartesian position in frame fixed to body exerting acceleration
        Eigen::Matrix3d currentRotationToBodyFixedFrame_ = fromBodyFixedToIntegrationFrameRotation_( ).inverse( );
        bodyFixedPosition_ = currentRotationToBodyFixedFrame_ *
                ( positionFunctionOfAcceleratedBody_( ) - positionFunctionOfAcceleratingBody_( ) );

        // Calculate spherical position in frame fixed to body exerting acceleration
        bodyFixedSphericalPosition_ = convertCartesianToSpherical( bodyFixedPosition_ );
        bodyFixedSphericalPosition_( 1 ) = mathematical_constants::PI / 2.0 - bodyFixedSphericalPosition_( 1 );

        // Calculate partial of acceleration wrt position of body undergoing acceleration.
        currentBodyFixedPartialWrtPosition_ = basic_mathematics::calculatePolyhedronHessianOfGravitationalPotential(
                gravitationalParameterFunction_() / volumeFunction_(),
                facetDyads_, edgeDyads_,
                polyhedronCache_->getPerFacetFactor(),
                polyhedronCache_->getPerEdgeFactor() );

        currentPartialWrtVelocity_.setZero( );
        currentPartialWrtPosition_.setZero( );

        // Compute partial w.r.t. position in inertial frame
        currentPartialWrtPosition_ +=
                currentRotationToBodyFixedFrame_.inverse( ) * currentBodyFixedPartialWrtPosition_ * currentRotationToBodyFixedFrame_;

        // If rotation matrix depends on translational state, add correction partials
        if( rotationMatrixPartials_.count(
                    std::make_pair( estimatable_parameters::initial_body_state, "" ) ) > 0 )
        {
            // Compute the acceleration and body-fixed position partial, without the central term (to avoid numerical errors)
            Eigen::Vector3d nonCentralAcceleration = accelerationFunction_( );
            Eigen::Matrix3d nonCentralBodyFixedPartial = currentBodyFixedPartialWrtPosition_;

            nonCentralAcceleration -= gravitation::computeGravitationalAcceleration(
                        positionFunctionOfAcceleratedBody_( ), gravitationalParameterFunction_( ), positionFunctionOfAcceleratingBody_( ) );

            nonCentralBodyFixedPartial -=
                    currentRotationToBodyFixedFrame_ * calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody(
                        positionFunctionOfAcceleratedBody_( ), positionFunctionOfAcceleratingBody_( ), gravitationalParameterFunction_( ) ) *
                    currentRotationToBodyFixedFrame_.inverse( );

            // Compute rotation matrix partials
            std::vector< Eigen::Matrix3d > rotationPositionPartials =
                    rotationMatrixPartials_.at(
                        std::make_pair( estimatable_parameters::initial_body_state, "" ) )->
                    calculatePartialOfRotationMatrixToBaseFrameWrParameter( currentTime );

            // Add correction terms to position and velocity partials
            for( unsigned int i = 0; i < 3; i++ )
            {
                currentPartialWrtPosition_.block( 0, i, 3, 1 ) -=
                        rotationPositionPartials.at( i ) * ( currentRotationToBodyFixedFrame_ * nonCentralAcceleration );

                currentPartialWrtPosition_.block( 0, i, 3, 1 ) -=
                        currentRotationToBodyFixedFrame_.inverse( ) * nonCentralBodyFixedPartial *
                        ( rotationPositionPartials.at( i ).transpose( ) *
                          ( positionFunctionOfAcceleratedBody_( ) - positionFunctionOfAcceleratingBody_( ) ) );
                currentPartialWrtVelocity_.block( 0, i, 3, 1 ) -=
                        rotationPositionPartials.at( i + 3 ) * ( currentRotationToBodyFixedFrame_ * nonCentralAcceleration );
            }
        }

        currentTime_ = currentTime;
    }

}

} // namespace acceleration_partials

} // namespace tudat
