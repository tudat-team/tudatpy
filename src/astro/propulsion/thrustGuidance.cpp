/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/propulsion/thrustGuidance.h"

namespace tudat
{

namespace propulsion
{

//! Function to get the unit vector colinear with velocity segment of a translational state.
Eigen::Vector3d getDirectionColinearWithVelocity(
        const std::function< void( Eigen::Vector6d& ) > currentStateFunction,
        const double currentTime, const bool putVectorInOppositeDirection )
{
    static Eigen::Vector6d currentState;
    currentStateFunction( currentState );
    return ( ( putVectorInOppositeDirection == 1 ) ? -1.0 : 1.0 ) * ( currentState.segment( 3, 3 ) ).normalized( );
}

//! Function to get the unit vector colinear with position segment of a translational state.
Eigen::Vector3d getDirectionColinearWithPosition(
        const std::function< void( Eigen::Vector6d& ) > currentStateFunction,
        const double currentTime, const bool putVectorInOppositeDirection )
{
    static Eigen::Vector6d currentState;
    currentStateFunction( currentState );
    return ( ( putVectorInOppositeDirection == 1 ) ? -1.0 : 1.0 ) * ( currentState.segment( 0, 3 ) ).normalized( );
}

//! Function to get the force direction from a time-only function.
Eigen::Vector3d getForceDirectionFromTimeOnlyFunction(
        const double currentTime,
        const std::function< Eigen::Vector3d( const double ) > timeOnlyFunction )
{
    return timeOnlyFunction( currentTime ).normalized( );
}

void DirectThrustDirectionCalculator::updateQuaternion( const double time )
{
    if( time != currentQuaterionTime_  && time == time )
    {
        currentRotationToBaseFrame_ = directionBasedRotationModel_->getRotationToBaseFrame( time );
    }
    currentQuaterionTime_ = time;
}

Eigen::Vector3d DirectThrustDirectionCalculator::getInertialThrustDirection(
        const std::shared_ptr< system_models::EngineModel > engineModel )
{
    if( engineModel->getBodyFixedThrustDirection( ) == directionBasedRotationModel_->getAssociatedBodyFixedDirection( ) )
    {
        return currentInertialDirection_;
    }
    else if( engineModel->getBodyFixedThrustDirection( ) == -directionBasedRotationModel_->getAssociatedBodyFixedDirection( ) )
    {
        return currentInertialDirection_;
    }
    else
    {
        updateQuaternion( currentTime_ );
        return ( currentRotationToBaseFrame_ * engineModel->getBodyFixedThrustDirection( ) ).normalized( );
    }
}


} // namespace propulsion

} // namespace tudat

