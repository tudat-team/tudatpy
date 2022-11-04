/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_IAUROTATIONMODEL_H
#define TUDAT_IAUROTATIONMODEL_H

#include <vector>
#include <map>

#include <boost/function.hpp>

#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/astro/reference_frames.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"

namespace tudat
{

namespace ephemerides
{

class IauRotationModel: public RotationalEphemeris
{
public:
    IauRotationModel(  const std::string& baseFrameOrientation,
                       const std::string& targetFrameOrientation,
                       const double nominalMeridian,
                       const Eigen::Vector2d& nominalPole,
                       const double rotationRate,
                       const Eigen::Vector2d& polePrecession,
                       const std::map< double, std::pair< double, double > >& meridianPeriodicTerms,
                       const std::map< double, std::pair< Eigen::Vector2d, double > >& polePeriodicTerms ):
        RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
        nominalMeridian_( nominalMeridian ), nominalPole_( nominalPole ),
        rotationRate_( rotationRate ), polePrecession_( polePrecession ),
        meridianPeriodicTerms_( meridianPeriodicTerms ),
        polePeriodicTerms_( polePeriodicTerms )

    { }

    Eigen::Quaterniond getRotationToBaseFrame(
            const double secondsSinceEpoch )
    {
        return getRotationToTargetFrame( secondsSinceEpoch ).inverse( );
    }

    Eigen::Quaterniond getRotationToTargetFrame(
            const double currentTime )
    {
        updateRotationMatrices( currentTime );
        return Eigen::Quaterniond(
                    meridianRotationAboutZAxis_ * declinationRotationAboutXAxis_ * rightAscensionRotationAboutZAxis_ );
    }

    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame(
            const double secondsSinceEpoch )
    {
        using namespace reference_frames;

        updateRotationMatrices( secondsSinceEpoch );
        updateAngleDerivatives( secondsSinceEpoch );

        return currentMeridianDerivative_ * getDerivativeOfZAxisRotationWrtAngle( meridianRotationAboutZAxis_ ) * declinationRotationAboutXAxis_ * rightAscensionRotationAboutZAxis_ -
                currentPoleDerivative_( 1 ) * meridianRotationAboutZAxis_ * getDerivativeOfXAxisRotationWrtAngle( declinationRotationAboutXAxis_ ) * rightAscensionRotationAboutZAxis_ +
                currentPoleDerivative_( 0 ) * meridianRotationAboutZAxis_ * declinationRotationAboutXAxis_ * getDerivativeOfZAxisRotationWrtAngle( rightAscensionRotationAboutZAxis_ );
    }

    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame(
            const double secondsSinceEpoch )
    {
        return getDerivativeOfRotationToTargetFrame( secondsSinceEpoch ).transpose( );
    }


    Eigen::Vector2d getNominalPole( )
    {
        return nominalPole_;
    }

    void setNominalPole( const Eigen::Vector2d& nominalPole )
    {
        nominalPole_ = nominalPole;
    }


    Eigen::Vector2d getPolePrecession( )
    {
        return polePrecession_;
    }

    void setPolePrecession( const Eigen::Vector2d& polePrecession )
    {
        polePrecession_ = polePrecession;
    }

    Eigen::Matrix3d getCurrentMeridianRotationAboutZAxis( )
    {
        return meridianRotationAboutZAxis_;
    }

    Eigen::Matrix3d getCurrentDeclinationRotationAboutXAxis( )
    {
        return declinationRotationAboutXAxis_;
    }

    Eigen::Matrix3d getCurrentRightAscensionRotationAboutZAxis( )
    {
        return rightAscensionRotationAboutZAxis_;
    }

    void updateRotationMatrices( const double currentTime )
    {
        if( currentTime == currentTime )
        {
            updateAngles( currentTime );
            meridianRotationAboutZAxis_ =
                    Eigen::AngleAxisd( -currentMeridian_,
                                       Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );

            declinationRotationAboutXAxis_ =
                    Eigen::AngleAxisd( -( mathematical_constants::PI / 2.0 - currentPolePosition_( 1 ) ),
                                       Eigen::Vector3d::UnitX( ) ).toRotationMatrix( );

            rightAscensionRotationAboutZAxis_ =
                    Eigen::AngleAxisd( -( currentPolePosition_( 0 ) + mathematical_constants::PI / 2.0 ),
                                       Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );
        }
    }

private:

    void updateAngles( const double currentTime )
    {
        if( currentTime == currentTime )
        {
            currentMeridian_ = nominalMeridian_ + rotationRate_ * currentTime;
            currentPolePosition_ = nominalPole_ + polePrecession_ * currentTime;

            for( auto it : meridianPeriodicTerms_ )
            {
                currentMeridian_ += it.second.first * std::sin(
                            it.first * currentTime +
                            it.second.second );
            }

            for( auto it : polePeriodicTerms_ )
            {
                double currentAngle = it.first * currentTime +
                        it.second.second;

                currentPolePosition_( 0 ) += it.second.first( 0 ) * std::sin( currentAngle );
                currentPolePosition_( 1 ) += it.second.first( 1 ) * std::cos( currentAngle );

            }
        }
    }


    void updateAngleDerivatives( const double currentTime )
    {
        if( currentTime == currentTime )
        {
            currentMeridianDerivative_ = rotationRate_;
            currentPoleDerivative_ = polePrecession_;

            for( auto it : meridianPeriodicTerms_ )
            {
                currentMeridianDerivative_ += it.second.first * it.first * std::cos(
                            it.first * currentTime +
                            it.second.second );
            }

            for( auto it : polePeriodicTerms_ )
            {
                double currentAngle = it.first * currentTime +
                        it.second.second;

                currentPoleDerivative_( 0 ) += it.first * it.second.first( 0 ) * std::cos( currentAngle );
                currentPoleDerivative_( 1 ) -= it.first * it.second.first( 1 ) * std::sin( currentAngle );

            }
        }
    }


    double nominalMeridian_;

    Eigen::Vector2d nominalPole_;

    double rotationRate_;

    Eigen::Vector2d polePrecession_;

    std::map< double, std::pair< double, double > > meridianPeriodicTerms_;

    std::map< double, std::pair< Eigen::Vector2d, double > > polePeriodicTerms_;

    double currentMeridian_;

    double currentMeridianDerivative_;

    Eigen::Vector2d currentPolePosition_;

    Eigen::Vector2d currentPoleDerivative_;

    Eigen::Matrix3d meridianRotationAboutZAxis_;

    Eigen::Matrix3d declinationRotationAboutXAxis_;

    Eigen::Matrix3d rightAscensionRotationAboutZAxis_;

};

}

}

#endif // TUDAT_IAUROTATIONMODEL_H
