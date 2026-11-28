/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/ephemerides/iauRotationModel.h"

namespace tudat
{

namespace ephemerides
{

Eigen::Quaterniond IauRotationModel::getRotationToBaseFrame( const double secondsSinceEpoch )
{
    return getRotationToTargetFrame( secondsSinceEpoch ).inverse( );
}

Eigen::Quaterniond IauRotationModel::getRotationToTargetFrame( const double currentTime )
{
    updateRotationMatrices( currentTime );
    return Eigen::Quaterniond( meridianRotationAboutZAxis_ * declinationRotationAboutXAxis_ * rightAscensionRotationAboutZAxis_ *
                               anglesFrameToBaseFrameRotation_ );
}

Eigen::Matrix3d IauRotationModel::getDerivativeOfRotationToTargetFrame( const double secondsSinceEpoch )
{
    using namespace reference_frames;

    updateRotationMatrices( secondsSinceEpoch );
    updateAngleDerivatives( secondsSinceEpoch );

    return currentMeridianDerivative_ * getDerivativeOfZAxisRotationWrtAngle( meridianRotationAboutZAxis_ ) *
            declinationRotationAboutXAxis_ * rightAscensionRotationAboutZAxis_ * anglesFrameToBaseFrameRotation_ -
            currentPoleDerivative_( 1 ) * meridianRotationAboutZAxis_ *
            getDerivativeOfXAxisRotationWrtAngle( declinationRotationAboutXAxis_ ) * rightAscensionRotationAboutZAxis_ *
            anglesFrameToBaseFrameRotation_ +
            currentPoleDerivative_( 0 ) * meridianRotationAboutZAxis_ * declinationRotationAboutXAxis_ *
            getDerivativeOfZAxisRotationWrtAngle( rightAscensionRotationAboutZAxis_ ) * anglesFrameToBaseFrameRotation_;
}

Eigen::Matrix3d IauRotationModel::getDerivativeOfRotationToBaseFrame( const double secondsSinceEpoch )
{
    return getDerivativeOfRotationToTargetFrame( secondsSinceEpoch ).transpose( );
}

void IauRotationModel::updateRotationMatrices( const double currentTime )
{
    if( currentTime == currentTime )
    {
        updateAngles( currentTime );
        meridianRotationAboutZAxis_ = Eigen::AngleAxisd( -currentMeridian_, Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );

        declinationRotationAboutXAxis_ =
                Eigen::AngleAxisd( -( mathematical_constants::PI / 2.0 - currentPolePosition_( 1 ) ), Eigen::Vector3d::UnitX( ) )
                        .toRotationMatrix( );

        rightAscensionRotationAboutZAxis_ =
                Eigen::AngleAxisd( -( currentPolePosition_( 0 ) + mathematical_constants::PI / 2.0 ), Eigen::Vector3d::UnitZ( ) )
                        .toRotationMatrix( );
    }
}

void IauRotationModel::updateAngles( const double currentTime )
{
    if( currentTime == currentTime )
    {
        epochsSinceReference_ = currentTime - referenceEpochJ2000_;
        currentMeridian_ = nominalMeridian_ + rotationRate_ * epochsSinceReference_;
        currentPolePosition_ = nominalPole_ + polePrecession_ * epochsSinceReference_;

        for( auto it : meridianPeriodicTerms_ )
        {
            currentMeridian_ += it.second.first * std::sin( it.first * epochsSinceReference_ + it.second.second );
        }

        for( auto it : polePeriodicTerms_ )
        {
            double currentAngle = it.first * epochsSinceReference_ + it.second.second;

            currentPolePosition_( 0 ) += it.second.first( 0 ) * std::sin( currentAngle );
            currentPolePosition_( 1 ) += it.second.first( 1 ) * std::cos( currentAngle );
        }
    }
}

void IauRotationModel::updateAngleDerivatives( const double currentTime )
{
    if( currentTime == currentTime )
    {
        epochsSinceReference_ = currentTime - referenceEpochJ2000_;

        currentMeridianDerivative_ = rotationRate_;
        currentPoleDerivative_ = polePrecession_;

        for( auto it : meridianPeriodicTerms_ )
        {
            currentMeridianDerivative_ += it.second.first * it.first * std::cos( it.first * epochsSinceReference_ + it.second.second );
        }

        for( auto it : polePeriodicTerms_ )
        {
            double currentAngle = it.first * epochsSinceReference_ + it.second.second;

            currentPoleDerivative_( 0 ) += it.first * it.second.first( 0 ) * std::cos( currentAngle );
            currentPoleDerivative_( 1 ) -= it.first * it.second.first( 1 ) * std::sin( currentAngle );
        }
    }
}

}  // namespace ephemerides

}  // namespace tudat
