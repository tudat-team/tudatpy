/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_DIRECTTIDALDISSIPATIONACCELERATION_H
#define TUDAT_DIRECTTIDALDISSIPATIONACCELERATION_H

#include <iostream>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Basics/basicTypedefs.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

namespace tudat
{

namespace gravitation
{

Eigen::Vector3d computeDirectTidalAccelerationDueToTideOnPlanet(
        const Eigen::Vector6d relativeStateOfBodyExertingTide, const Eigen::Vector3d planetAngularVelocityVector,
        const double currentTidalAccelerationMultiplier, const double timeLag,
        const bool includeDirectRadialComponent = true );

Eigen::Vector3d computeDirectTidalAccelerationDueToTideOnSatellite(
        const Eigen::Vector6d relativeStateOfBodyExertingTide,
        const double currentTidalAccelerationMultiplier,
        const double timeLag,const bool includeDirectRadialComponent );

class DirectTidalDissipationAcceleration: public basic_astrodynamics::AccelerationModel3d
{
public:
    DirectTidalDissipationAcceleration(
            const boost::function< Eigen::Vector6d( ) > stateFunctionOfBodyExertingTide,
            const boost::function< Eigen::Vector6d( ) > stateFunctionOfBodyUndergoingTide_,
            const boost::function< double( ) > massFunctionOfBodyExertingTide,
            const boost::function< Eigen::Vector3d( ) > angularVelocityVectorOfBodyUndergoingTide,
            const double k2LoveNumber,
            const double timeLag,
            const double equatorialRadiusOfBodyUndergoingTide,
            const bool includeDirectRadialComponent ):
        stateFunctionOfBodyExertingTide_( stateFunctionOfBodyExertingTide ), stateFunctionOfBodyUndergoingTide_( stateFunctionOfBodyUndergoingTide_ ),
        massFunctionOfBodyExertingTide_( massFunctionOfBodyExertingTide ),
        angularVelocityVectorOfBodyUndergoingTide_( angularVelocityVectorOfBodyUndergoingTide ),
        k2LoveNumber_( k2LoveNumber ), timeLag_( timeLag ), equatorialRadiusOfBodyUndergoingTide_( equatorialRadiusOfBodyUndergoingTide ),
        includeDirectRadialComponent_( includeDirectRadialComponent ), modelTideOnPlanet_( true )
    {
        equatorialRadiusToFifthPower_ =
                equatorialRadiusOfBodyUndergoingTide_ * equatorialRadiusOfBodyUndergoingTide_ *
                equatorialRadiusOfBodyUndergoingTide_ * equatorialRadiusOfBodyUndergoingTide_ *
                equatorialRadiusOfBodyUndergoingTide_;
    }

    DirectTidalDissipationAcceleration(
            const boost::function< Eigen::Vector6d( ) > stateFunctionOfBodyExertingTide,
            const boost::function< Eigen::Vector6d( ) > stateFunctionOfBodyUndergoingTide_,
            const boost::function< double( ) > massFunctionOfBodyExertingTide,
            const boost::function< double( ) > massFunctionOfBodyUndergoingTide,
            const double k2LoveNumber,
            const double timeLag,
            const double equatorialRadiusOfBodyUndergoingTide,
            const bool includeDirectRadialComponent ):
        stateFunctionOfBodyExertingTide_( stateFunctionOfBodyExertingTide ), stateFunctionOfBodyUndergoingTide_( stateFunctionOfBodyUndergoingTide_ ),
        massFunctionOfBodyExertingTide_( massFunctionOfBodyExertingTide ), massFunctionOfBodyUndergoingTide_( massFunctionOfBodyUndergoingTide ),
        angularVelocityVectorOfBodyUndergoingTide_( boost::lambda::constant( Eigen::Vector3d::Constant( TUDAT_NAN ) ) ),
        k2LoveNumber_( k2LoveNumber ), timeLag_( timeLag ), equatorialRadiusOfBodyUndergoingTide_( equatorialRadiusOfBodyUndergoingTide ),
        includeDirectRadialComponent_( includeDirectRadialComponent ), modelTideOnPlanet_( false )
    {
        equatorialRadiusToFifthPower_ =
                equatorialRadiusOfBodyUndergoingTide_ * equatorialRadiusOfBodyUndergoingTide_ *
                equatorialRadiusOfBodyUndergoingTide_ * equatorialRadiusOfBodyUndergoingTide_ *
                equatorialRadiusOfBodyUndergoingTide_;
    }

    ~DirectTidalDissipationAcceleration( ){ }

    Eigen::Vector3d getAcceleration( )
    {
        return currentAcceleration_;
    }

    //! Update member variables used by the aerodynamic acceleration model.
    /*!
    * Updates member variables used by the aerodynamic acceleration model.
    * Function pointers to retrieve the current values of quantities from which the
    * acceleration is to be calculated are set by constructor. This function calls
    * them to update the associated variables to their current state.
    * \param currentTime Time at which acceleration model is to be updated.
    */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            currentRelativeState_ = stateFunctionOfBodyExertingTide_( ) - stateFunctionOfBodyUndergoingTide_( );
            double distance = currentRelativeState_.segment( 0, 3 ).norm( );
            double distanceSquared = distance * distance;
            double distanceToEighthPower = distanceSquared * distanceSquared * distanceSquared * distanceSquared;
            currentTidalAccelerationMultiplier_ =
                    - 3.0 * massFunctionOfBodyExertingTide_( ) * equatorialRadiusToFifthPower_ / distanceToEighthPower * k2LoveNumber_;

            if( modelTideOnPlanet_ )
            {
                currentAngularVelocityVectorOfBodyUndergoingTide_ = angularVelocityVectorOfBodyUndergoingTide_( );


                currentAcceleration_ = computeDirectTidalAccelerationDueToTideOnPlanet(
                            currentRelativeState_, currentAngularVelocityVectorOfBodyUndergoingTide_,
                            currentTidalAccelerationMultiplier_, timeLag_, includeDirectRadialComponent_ );
            }
            else
            {
                currentTidalAccelerationMultiplier_ *= massFunctionOfBodyExertingTide_( ) / massFunctionOfBodyUndergoingTide_( );

                currentAcceleration_ = computeDirectTidalAccelerationDueToTideOnSatellite(
                            currentRelativeState_, currentTidalAccelerationMultiplier_,
                            timeLag_, includeDirectRadialComponent_ );
            }

            currentTime_ = currentTime;
        }
    }


    Eigen::Vector6d getCurrentRelativeState( )
    {
        return currentRelativeState_;
    }

    Eigen::Vector3d getCurrentAngularVelocityVectorOfBodyUndergoingTide( )
    {
        return currentAngularVelocityVectorOfBodyUndergoingTide_;
    }



    boost::function< double( ) > getMassFunctionOfBodyExertingTide( )
    {
        return massFunctionOfBodyExertingTide_;
    }

    boost::function< double( ) > getMassFunctionOfBodyUndergoingTide( )
    {
        return massFunctionOfBodyUndergoingTide_;
    }




    double getK2LoveNumber( )
    {
        return k2LoveNumber_;
    }

    double getTimeLag( )
    {
        return timeLag_;
    }

    double getEquatorialRadiusOfBodyUndergoingTide( )
    {
        return equatorialRadiusOfBodyUndergoingTide_;
    }

    bool getIncludeDirectRadialComponent( )
    {
        return includeDirectRadialComponent_;
    }

    bool getModelTideOnPlanet( )
    {
        return modelTideOnPlanet_;
    }


    Eigen::Vector3d getCurrentAcceleration( )
    {
        return currentAcceleration_;
    }

    double getCurrentTidalAccelerationMultiplier( )
    {
        return currentTidalAccelerationMultiplier_;
    }

private:

    Eigen::Vector6d currentRelativeState_;

    Eigen::Vector3d currentAngularVelocityVectorOfBodyUndergoingTide_;


    boost::function< Eigen::Vector6d( ) > stateFunctionOfBodyExertingTide_;

    boost::function< Eigen::Vector6d( ) > stateFunctionOfBodyUndergoingTide_;

    boost::function< double( ) > massFunctionOfBodyExertingTide_;

    boost::function< double( ) > massFunctionOfBodyUndergoingTide_;

    boost::function< Eigen::Vector3d( ) > angularVelocityVectorOfBodyUndergoingTide_;


    double k2LoveNumber_;

    double timeLag_;

    double equatorialRadiusOfBodyUndergoingTide_;

    double equatorialRadiusToFifthPower_;

    bool includeDirectRadialComponent_;

    bool modelTideOnPlanet_;


    Eigen::Vector3d currentAcceleration_;

    double currentTidalAccelerationMultiplier_;
};

} // namespace gravitation

} // namespace tudat

#endif // TUDAT_DIRECTTIDALDISSIPATIONACCELERATION_H
