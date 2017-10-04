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
        const Eigen::Vector6d relativeStateOfBodyExertingTide, const Eigen::Vector3d angularVelocityVectorOfBodyUndergoingTide,
        const double massOfBodyExertingTide, const double k2LoveNumber, const double timeLag, const double referenceRadius,
        const bool includeDirectRadialComponent = true );

Eigen::Vector3d computeDirectTidalAccelerationDueToTideOnSatellite(
        const Eigen::Vector6d relativeStateOfBodyExertingTide, const double massOfBodyExertingTide, const double k2LoveNumber,
        const double timeLag, const double referenceRadius, const bool includeDirectRadialComponent = true );

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
        stateFunctionOfBodyExertingTide_( stateFunctionOfBodyExertingTide ), stateFunctionOfBodyUndergoingTide__( stateFunctionOfBodyUndergoingTide_ ),
        massFunctionOfBodyExertingTide_( massFunctionOfBodyExertingTide ),
        angularVelocityVectorOfBodyUndergoingTide_( angularVelocityVectorOfBodyUndergoingTide ),
        k2LoveNumber_( k2LoveNumber ), timeLag_( timeLag ), equatorialRadiusOfBodyUndergoingTide_( equatorialRadiusOfBodyUndergoingTide ),
        includeDirectRadialComponent_( includeDirectRadialComponent ), useBodyRotationTerm_( true )
    { }

    DirectTidalDissipationAcceleration(
            const boost::function< Eigen::Vector6d( ) > stateFunctionOfBodyExertingTide,
            const boost::function< Eigen::Vector6d( ) > stateFunctionOfBodyUndergoingTide_,
            const boost::function< double( ) > massFunctionOfBodyExertingTide,
            const double k2LoveNumber,
            const double timeLag,
            const double equatorialRadiusOfBodyUndergoingTide,
            const bool includeDirectRadialComponent ):
        stateFunctionOfBodyExertingTide_( stateFunctionOfBodyExertingTide ), stateFunctionOfBodyUndergoingTide__( stateFunctionOfBodyUndergoingTide_ ),
        massFunctionOfBodyExertingTide_( massFunctionOfBodyExertingTide ),
        angularVelocityVectorOfBodyUndergoingTide_( boost::lambda::constant( Eigen::Vector3d::Constant( TUDAT_NAN ) ) ),
        k2LoveNumber_( k2LoveNumber ), timeLag_( timeLag ), equatorialRadiusOfBodyUndergoingTide_( equatorialRadiusOfBodyUndergoingTide ),
        includeDirectRadialComponent_( includeDirectRadialComponent ), useBodyRotationTerm_( false )
    { }

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
            if( useBodyRotationTerm_ )
            {
            currentAcceleration_ = computeDirectTidalAccelerationDueToTideOnPlanet(
                        stateFunctionOfBodyExertingTide_( ) - stateFunctionOfBodyUndergoingTide__( ), angularVelocityVectorOfBodyUndergoingTide_( ),
                        massFunctionOfBodyExertingTide_( ),
                        k2LoveNumber_, timeLag_, equatorialRadiusOfBodyUndergoingTide_, includeDirectRadialComponent_ );
            }
            else
            {
                currentAcceleration_ = computeDirectTidalAccelerationDueToTideOnSatellite(
                            stateFunctionOfBodyExertingTide_( ) - stateFunctionOfBodyUndergoingTide__( ),
                            massFunctionOfBodyExertingTide_( ),
                            k2LoveNumber_, timeLag_, equatorialRadiusOfBodyUndergoingTide_, includeDirectRadialComponent_ );
            }

            currentTime_ = currentTime;
        }
    }

private:

    boost::function< Eigen::Vector6d( ) > stateFunctionOfBodyExertingTide_;

    boost::function< Eigen::Vector6d( ) > stateFunctionOfBodyUndergoingTide__;

    boost::function< double( ) > massFunctionOfBodyExertingTide_;

    boost::function< Eigen::Vector3d( ) > angularVelocityVectorOfBodyUndergoingTide_;

    double k2LoveNumber_;

    double timeLag_;

    double equatorialRadiusOfBodyUndergoingTide_;

    bool includeDirectRadialComponent_;

    bool useBodyRotationTerm_;

    Eigen::Vector3d currentAcceleration_;
};

} // namespace gravitation

} // namespace tudat

#endif // TUDAT_DIRECTTIDALDISSIPATIONACCELERATION_H
