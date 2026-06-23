/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/observation_partials/pixelCoordinatesPointingPartial.h"

#include <cmath>

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace observation_partials
{

namespace
{

Eigen::Matrix3d getLeftJacobianOfSO3( const Eigen::Vector3d& rotationVector )
{
    const double angle = rotationVector.norm( );
    const Eigen::Matrix3d rotationVectorCrossProduct = linear_algebra::getCrossProductMatrix( rotationVector );
    const Eigen::Matrix3d rotationVectorCrossProductSquared = rotationVectorCrossProduct * rotationVectorCrossProduct;

    if( angle < 1.0E-8 )
    {
        return Eigen::Matrix3d::Identity( ) + 0.5 * rotationVectorCrossProduct + ( 1.0 / 6.0 ) * rotationVectorCrossProductSquared;
    }

    return Eigen::Matrix3d::Identity( ) + ( 1.0 - std::cos( angle ) ) / ( angle * angle ) * rotationVectorCrossProduct +
            ( angle - std::sin( angle ) ) / ( angle * angle * angle ) * rotationVectorCrossProductSquared;
}

}  // namespace

PixelCoordinatesPointingPartial::PixelCoordinatesPointingPartial(
        const std::shared_ptr< PixelCoordinatesScaling > pixelCoordinatesScaling,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > pointingCorrectionParameter ):
    ObservationPartial< 2 >( pointingCorrectionParameter->getParameterName( ) ), pixelCoordinatesScaling_( pixelCoordinatesScaling ),
    pointingCorrectionParameter_( pointingCorrectionParameter )
{
    if( pixelCoordinatesScaling_ == nullptr )
    {
        throw std::runtime_error( "Error when creating pixel-coordinate pointing partial: pixel-coordinate scaling is null." );
    }
    if( pointingCorrectionParameter_ == nullptr )
    {
        throw std::runtime_error( "Error when creating pixel-coordinate pointing partial: pointing parameter is null." );
    }
    receiverIndex_ =
            observation_models::getSingleLinkStateEntryIndices( observation_models::pixel_coordinates ).at( observation_models::receiver );
}

std::vector< std::pair< Eigen::Matrix< double, 2, Eigen::Dynamic >, double > > PixelCoordinatesPointingPartial::calculatePartial(
        const std::vector< Eigen::Vector6d >&,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings >,
        const Eigen::Matrix< double, 2, 1 >& )
{
    if( linkEndOfFixedTime != observation_models::receiver )
    {
        throw std::runtime_error( "Error when computing pixel-coordinate pointing partial: fixed link end must be receiver." );
    }

    const Eigen::Vector3d relativeRangeVectorCameraFrame = pixelCoordinatesScaling_->getCurrentRelativeRangeVectorCameraFrame( );
    const Eigen::Matrix< double, 2, 3 > pixelLinePartialWrtCameraFramePosition =
            pixelCoordinatesScaling_->getCurrentPixelLinePartialWrtCameraFramePosition( );

    const Eigen::Vector3d pointingCorrection = pointingCorrectionParameter_->getParameterValue( ).segment( 0, 3 );

    // For r_cam = Exp(delta_theta) * r_nominal and additive perturbations of delta_theta:
    // d r_cam / d delta_theta = -skew(r_cam) * J_left(delta_theta).
    // At zero correction, J_left = I and this reduces to -skew(r_cam).
    const Eigen::Matrix3d relativeRangePartialWrtPointingCorrection =
            -linear_algebra::getCrossProductMatrix( relativeRangeVectorCameraFrame ) * getLeftJacobianOfSO3( pointingCorrection );

    Eigen::Matrix< double, 2, Eigen::Dynamic > partial = Eigen::Matrix< double, 2, Eigen::Dynamic >::Zero( 2, 3 );
    partial = pixelLinePartialWrtCameraFramePosition * relativeRangePartialWrtPointingCorrection;

    return { std::make_pair( partial, times.at( receiverIndex_ ) ) };
}

}  // namespace observation_partials

}  // namespace tudat
