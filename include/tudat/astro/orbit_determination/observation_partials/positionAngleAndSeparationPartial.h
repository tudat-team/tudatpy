/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_POSITIONANGLEANDSEPARATIONPARTIAL_H
#define TUDAT_POSITIONANGLEANDSEPARATIONPARTIAL_H

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/positionPartials.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/observation_partials/lightTimeCorrectionPartial.h"

namespace tudat
{

namespace observation_partials
{

// Forward declaration for PS scaling, used by PositionAngleScaling and SeparationScaling
class PositionAngleAndSeparationScaling;

//! Derived class for scaling three-dimensional position partial to position angle observable partial (size 1).
/*!
 *  Delegates to an internal PositionAngleAndSeparationScaling and extracts the first row (position angle).
 */
class PositionAngleScaling : public DirectPositionPartialScaling< 1 >
{
public:
    PositionAngleScaling( );

    ~PositionAngleScaling( );

    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation );

    Eigen::Matrix< double, 1, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        if( linkEndType == observation_models::transmitter )
            return referenceScalingFactorFirstTransmitter_;
        else if( linkEndType == observation_models::transmitter2 )
            return referenceScalingFactorSecondTransmitter_;
        else if( linkEndType == observation_models::receiver )
            return -( referenceScalingFactorFirstTransmitter_ + referenceScalingFactorSecondTransmitter_ );
        else
            throw std::runtime_error( "Error when getting position angle scaling factor, incorrect link end type." );
    }

    Eigen::Matrix< double, 1, 1 > getLightTimePartialScalingFactor( )
    {
        return referenceLightTimeCorrectionScaling_;
    }

    observation_models::LinkEndType getCurrentLinkEndType( )
    {
        return currentLinkEndType_;
    }

private:
    Eigen::Matrix< double, 1, 3 > referenceScalingFactorFirstTransmitter_;
    Eigen::Matrix< double, 1, 3 > referenceScalingFactorSecondTransmitter_;
    Eigen::Matrix< double, 1, 1 > referenceLightTimeCorrectionScaling_;
    observation_models::LinkEndType currentLinkEndType_;

    //! Internal PS scaling that performs the full computation
    std::shared_ptr< PositionAngleAndSeparationScaling > psScaling_;
};

//! Derived class for scaling three-dimensional position partial to angular separation observable partial (size 1).
/*!
 *  Delegates to an internal PositionAngleAndSeparationScaling and extracts the second row (angular separation).
 */
class SeparationScaling : public DirectPositionPartialScaling< 1 >
{
public:
    SeparationScaling( );

    ~SeparationScaling( );

    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation );

    Eigen::Matrix< double, 1, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        if( linkEndType == observation_models::transmitter )
            return referenceScalingFactorFirstTransmitter_;
        else if( linkEndType == observation_models::transmitter2 )
            return referenceScalingFactorSecondTransmitter_;
        else if( linkEndType == observation_models::receiver )
            return -( referenceScalingFactorFirstTransmitter_ + referenceScalingFactorSecondTransmitter_ );
        else
            throw std::runtime_error( "Error when getting separation scaling factor, incorrect link end type." );
    }

    Eigen::Matrix< double, 1, 1 > getLightTimePartialScalingFactor( )
    {
        return referenceLightTimeCorrectionScaling_;
    }

    observation_models::LinkEndType getCurrentLinkEndType( )
    {
        return currentLinkEndType_;
    }

private:
    Eigen::Matrix< double, 1, 3 > referenceScalingFactorFirstTransmitter_;
    Eigen::Matrix< double, 1, 3 > referenceScalingFactorSecondTransmitter_;
    Eigen::Matrix< double, 1, 1 > referenceLightTimeCorrectionScaling_;
    observation_models::LinkEndType currentLinkEndType_;

    //! Internal PS scaling that performs the full computation
    std::shared_ptr< PositionAngleAndSeparationScaling > psScaling_;
};

//! Derived class for scaling three-dimensional position partial to position angle and separation observable partial (size 2).
class PositionAngleAndSeparationScaling : public DirectPositionPartialScaling< 2 >
{
public:
    PositionAngleAndSeparationScaling( ): DirectPositionPartialScaling< 2 >( observation_models::position_angle_and_separation ) {}

    ~PositionAngleAndSeparationScaling( ) {}

    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation );

    Eigen::Matrix< double, 2, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType )
    {
        if( linkEndType == observation_models::transmitter )
            return referenceScalingFactorFirstTransmitter_;
        else if( linkEndType == observation_models::transmitter2 )
            return referenceScalingFactorSecondTransmitter_;
        else if( linkEndType == observation_models::receiver )
            return -( referenceScalingFactorFirstTransmitter_ + referenceScalingFactorSecondTransmitter_ );
        else
            throw std::runtime_error( "Error when getting position angle and separation scaling factor, incorrect link end type." );
    }

    Eigen::Vector2d getLightTimePartialScalingFactor( )
    {
        return referenceLightTimeCorrectionScaling_;
    }

    observation_models::LinkEndType getCurrentLinkEndType( )
    {
        return currentLinkEndType_;
    }

private:
    Eigen::Matrix< double, 2, 3 > referenceScalingFactorFirstTransmitter_;
    Eigen::Matrix< double, 2, 3 > referenceScalingFactorSecondTransmitter_;
    Eigen::Vector2d referenceLightTimeCorrectionScaling_;
    observation_models::LinkEndType currentLinkEndType_;
};

}  // namespace observation_partials

}  // namespace tudat

#endif  // TUDAT_POSITIONANGLEANDSEPARATIONPARTIAL_H
