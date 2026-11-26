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

#include "tudat/astro/ephemerides/aeordynamicAngleRotationalEphemeris.h"
#include "tudat/math/basic/rotationRepresentations.h"

namespace tudat
{

namespace reference_frames
{

Eigen::Vector3d computeBodyFixedAeroAngles( const Eigen::Matrix3d& inertialToBodyFixedFrame,
                                            const Eigen::Matrix3d& trajectoryToInertialFrame )
{
    // Retrieve rotation matrix that is to be converted to orientation angles.
    Eigen::Matrix3d currentRotationFromBodyToTrajectoryFrame_ = ( inertialToBodyFixedFrame * trajectoryToInertialFrame ).transpose( );

    // Compute associated Euler angles and set as orientation angles.
    Eigen::Vector3d eulerAngles = basic_mathematics::get132EulerAnglesFromRotationMatrix( currentRotationFromBodyToTrajectoryFrame_ );
    return ( Eigen::Vector3d( ) << -eulerAngles( 2 ), eulerAngles( 1 ), eulerAngles( 0 ) ).finished( );
}
Eigen::Vector3d FromGenericEphemerisAerodynamicAngleInterface::getAngles( const double time,
                                                                          const Eigen::Matrix3d& trajectoryToInertialFrame )
{
    return computeBodyFixedAeroAngles( ephemeris_->getRotationMatrixToTargetFrame( time ), trajectoryToInertialFrame );
}

void FromGenericEphemerisAerodynamicAngleInterface::resetCurrentTime( )
{
    ephemeris_->resetCurrentTime( );
}

Eigen::Vector3d FromAeroEphemerisAerodynamicAngleInterface::getAngles( const double time, const Eigen::Matrix3d& trajectoryToInertialFrame )
{
    ephemeris_->update( time );
    return ephemeris_->getBodyAngles( time );
}

void FromAeroEphemerisAerodynamicAngleInterface::resetCurrentTime( )
{
    ephemeris_->resetCurrentTime( );
}

}  // namespace reference_frames

namespace ephemerides
{

AerodynamicAngleRotationalEphemeris::AerodynamicAngleRotationalEphemeris(
        const std::string& centralBodyName,
        const std::string& baseFrameOrientation,
        const std::string& targetFrameOrientation,
        const std::function< Eigen::Vector3d( const double ) > aerodynamicAngleFunction ):
    RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ), centralBodyName_( centralBodyName ),
    aerodynamicAngleFunction_( aerodynamicAngleFunction ), currentTime_( TUDAT_NAN ), isBodyInPropagation_( false )
{}

void AerodynamicAngleRotationalEphemeris::applyEnvironmentModelClosure(
        const std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator )
{
    aerodynamicAngleCalculator_ = aerodynamicAngleCalculator;
    if( aerodynamicAngleCalculator_ != nullptr )
    {
        aerodynamicAngleCalculator_->setAerodynamicAngleClosureIsIncomplete( );
    }
}

void AerodynamicAngleRotationalEphemeris::setIsPitchTrimmed(
        const bool isPitchTrimmed,
        const std::function< Eigen::Vector2d( const double ) > sideslipAndBankAngleFunction )
{
    isPitchTrimmed_ = isPitchTrimmed;
    trimSideslipAndBankAngleFunction_ = sideslipAndBankAngleFunction;
}

void AerodynamicAngleRotationalEphemeris::updateBodyAngles( )
{
    if( aerodynamicAngleFunction_ != nullptr )
    {
        currentBodyAngles_ = aerodynamicAngleFunction_( currentTime_ );
    }
    else
    {
        currentBodyAngles_.setZero( );
    }
}

void AerodynamicAngleRotationalEphemeris::update( const double currentTime )
{
    if( !isBodyInPropagation_ )
    {
        throw std::runtime_error( "Error, can only use AerodynamicAngleRotationalEphemeris during propagation" );
    }
    else
    {
        if( !( currentTime == currentTime_ ) )
        {
            if( currentTime == currentTime )
            {
                if( aerodynamicAngleCalculator_ == nullptr )
                {
                    throw std::runtime_error(
                            "Error in aerodynamic angle based rotation model, model is not closed (no aerodynamic angle calculator)" );
                }
                if( isPitchTrimmed_ && aerodynamicAngleFunction_ == nullptr )
                {
                    throw std::runtime_error(
                            "Error in trimmed aerodynamic angle based rotation model, model is not closed (trim conditions not set, these "
                            "require closed atmospheric flight conditions of the body)" );
                }
                currentTime_ = currentTime;
                aerodynamicAngleCalculator_->update( currentTime, false );
                updateBodyAngles( );
                aerodynamicAngleCalculator_->update( currentTime, true );
            }
            else
            {
                resetCurrentTime( );
            }
        }
    }
}

void AerodynamicAngleRotationalEphemeris::resetCurrentTime( )
{
    if( currentTime_ == currentTime_ )
    {
        currentTime_ = TUDAT_NAN;
        if( aerodynamicAngleCalculator_ != nullptr )
        {
            aerodynamicAngleCalculator_->resetCurrentTime( );
        }
        if( aerodynamicAngleFunction_ != nullptr )
        {
            aerodynamicAngleFunction_( TUDAT_NAN );
        }
    }
    else
    {
        currentTime_ = TUDAT_NAN;
    }
}

}  // namespace ephemerides
}  // namespace tudat
