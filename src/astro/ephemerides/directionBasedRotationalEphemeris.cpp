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

#include "tudat/astro/ephemerides/directionBasedRotationalEphemeris.h"

namespace tudat
{


namespace ephemerides
{

void CustomBodyFixedDirectionCalculator::update( const double time )
{
    if( time != currentTime_ )
    {
        currentBodyAxisDirection_ = inertialBodyAxisDirectionFunction_( time );
        currentTime_ =  time;
    }
}

Eigen::Quaterniond DirectionBasedRotationalEphemeris::getRotationToBaseFrame(
        const double currentTime )
{
    Eigen::Vector3d eulerAngles = getEulerAngles( currentTime );
    return Eigen::Quaterniond(
                  Eigen::AngleAxisd( eulerAngles( 0 ), Eigen::Vector3d::UnitZ( ) )  *
                  Eigen::AngleAxisd( eulerAngles( 1 ), Eigen::Vector3d::UnitY( ) )  *
                  Eigen::AngleAxisd( eulerAngles( 2 ), Eigen::Vector3d::UnitX( ) ) );
}

void DirectionBasedRotationalEphemeris::update( const double currentTime )
{
    if( !( currentTime == currentTime_ ) )
    {
        if( currentTime == currentTime )
        {
            currentTime_ = currentTime;
            if( directionCalculator_ != nullptr )
            {
                currentInertialDirection_ = directionCalculator_->getInertialDirection( currentTime_ ).normalized( );
            }
            else
            {
                throw std::runtime_error( "Error when using DirectionBasedRotationalEphemeris, inerial body axis direction function no set" );
            }
        }
        else
        {
            resetCurrentTime( );
        }
    }
}


void DirectionBasedRotationalEphemeris::resetCurrentTime(  )
{
    currentTime_ = TUDAT_NAN;
    currentInertialDirection_.setConstant( TUDAT_NAN );

    currentEulerAnglesTime_ = TUDAT_NAN;
    eulerAngles_.setConstant( TUDAT_NAN );

    directionCalculator_->resetCurrentTime( );
}

void DirectionBasedRotationalEphemeris::calculateEulerAngles( )
{
    eulerAngles_( 0 ) = std::atan2( currentInertialDirection_.y( ), currentInertialDirection_.x( ) );
    eulerAngles_( 1 ) = -std::atan2( currentInertialDirection_.z( ), currentInertialDirection_.segment( 0, 2 ).norm( ) );
    if( freeRotationAngleFunction_ == nullptr )
    {
        eulerAngles_( 2 ) = 0.0;
    }
    else
    {
        eulerAngles_( 2 ) = freeRotationAngleFunction_( currentTime_ );
    }

    currentEulerAnglesTime_ = currentTime_;
}

} // namespace tudat
} // namespace ephemerides
