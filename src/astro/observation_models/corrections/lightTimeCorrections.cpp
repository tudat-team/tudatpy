/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"

namespace tudat
{

namespace observation_models
{

bool requiresMultiLegIterations( const LightTimeCorrectionType& lightTimeCorrectionType )
{
    bool requiresMultiLegIterations = false;
    switch( lightTimeCorrectionType )
    {
        case first_order_relativistic:
        case tabulated_tropospheric:
        case saastamoinen_tropospheric:
            requiresMultiLegIterations = false;
            break;
        case tabulated_ionospheric:
        case jakowski_vtec_ionospheric:
        case inverse_power_series_solar_corona:
            requiresMultiLegIterations = true;
            break;
        default:
            throw std::runtime_error(
                    "Error when getting whether light time corrections require multi-leg iterations: could not find light"
                    "time correction type " +
                    std::to_string( lightTimeCorrectionType ) + "." );
    }

    return requiresMultiLegIterations;
}

std::string getLightTimeCorrectionName( const LightTimeCorrectionType& lightTimeCorrectionType )
{
    std::string name;

    switch( lightTimeCorrectionType )
    {
        case first_order_relativistic:
            name = "first order relativistic";
            break;
        case tabulated_tropospheric:
            name = "tabulated tropospheric";
            break;
        case tabulated_ionospheric:
            name = "tabulated ionospheric";
            break;
        case saastamoinen_tropospheric:
            name = "Saastamoinen tropospheric";
            break;
        case jakowski_vtec_ionospheric:
            name = "Jakowski VTEC ionospheric";
            break;
        case inverse_power_series_solar_corona:
            name = "inverse power series solar corona";
            break;
        default:
            throw std::runtime_error(
                    "Error when getting light time correction name: could not find light"
                    "time correction type " +
                    std::to_string( lightTimeCorrectionType ) + "." );
    }

    return name;
}


double LightTimeCorrection::calculateLightTimeCorrectionPartialDerivativeWrtLinkEndTime( const Eigen::Vector6d& transmitterState,
                                                                                                  const Eigen::Vector6d& receiverState,
                                                                                                  const double transmissionTime,
                                                                                                  const double receptionTime,
                                                                                                  const LinkEndType linkEndAtWhichPartialIsEvaluated )
{
    double upPerturbedCorrection = 0.0, downPerturbedCorrection = 0.0;
    if( ( linkEndAtWhichPartialIsEvaluated == transmitter ) )
    {
        double upPerturbedTransmissionTime = transmissionTime + timePerturbation_;
        upPerturbedCorrection = calculateLightTimeCorrection(
            transmitterState, receiverState, upPerturbedTransmissionTime, receptionTime );

        double downPerturbedTransmissionTime = transmissionTime - timePerturbation_;
        downPerturbedCorrection = calculateLightTimeCorrection(
            transmitterState, receiverState, downPerturbedTransmissionTime, receptionTime );

    }
    else if( ( linkEndAtWhichPartialIsEvaluated == receiver ) )
    {
        double upPerturbedReceptionTime = receptionTime + timePerturbation_;
        upPerturbedCorrection = calculateLightTimeCorrection(
            transmitterState, receiverState, transmissionTime, upPerturbedReceptionTime );

        double downPerturbedReceptionTime = receptionTime - timePerturbation_;
        downPerturbedCorrection = calculateLightTimeCorrection(
            transmitterState, receiverState, transmissionTime, downPerturbedReceptionTime );
    }

    return ( upPerturbedCorrection - downPerturbedCorrection ) / ( 2.0 * timePerturbation_ );
}

/*!
 * Function to compute the partial derivative of the light-time correction w.r.t. link end position. Partial is
 * currently not implemented, function returns 0.
 *
 * \param transmitterState State of transmitted at transmission time
 * \param receiverState State of receiver at reception time
 * \param transmissionTime Time of signal transmission
 * \param receptionTime Time of singal reception
 * \param linkEndAtWhichPartialIsEvaluated Link end at which the position partial is to be taken
 * \return Partial of ight-time correction w.r.t. link end position
 */
Eigen::Matrix< double, 3, 1 > LightTimeCorrection::calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
    const Eigen::Vector6d& transmitterState,
    const Eigen::Vector6d& receiverState,
    const double transmissionTime,
    const double receptionTime,
    const LinkEndType linkEndAtWhichPartialIsEvaluated )
{
    Eigen::Matrix< double, 3, 1 > positionPartial = Eigen::Matrix< double, 3, 1 >::Zero( );

    double positionPerturbation =
        ( linkEndAtWhichPartialIsEvaluated == receiver ) ?
        ( positionRelativePerturbation_ * receiverState.segment( 0, 3 ).norm( ) ) :
        ( positionRelativePerturbation_ * transmitterState.segment( 0, 3 ).norm( ) );

    Eigen::Vector6d perturbedState;
    for( int i = 0; i < 3; i++ )
    {
        double upPerturbedCorrection, downPerturbedCorrection;
        if ( linkEndAtWhichPartialIsEvaluated == receiver )
        {
            perturbedState = receiverState;
            perturbedState( i ) += positionPerturbation;
            upPerturbedCorrection = calculateLightTimeCorrection(
                transmitterState, perturbedState, transmissionTime, receptionTime );

            perturbedState = receiverState;
            perturbedState( i ) -= positionPerturbation;
            downPerturbedCorrection = calculateLightTimeCorrection(
                transmitterState, perturbedState, transmissionTime, receptionTime );

        }
        else
        {
            perturbedState = transmitterState;
            perturbedState( i ) += positionPerturbation;
            upPerturbedCorrection = calculateLightTimeCorrection(
                perturbedState, receiverState, transmissionTime, receptionTime );

            perturbedState = transmitterState;
            perturbedState( i ) -= positionPerturbation;
            downPerturbedCorrection = calculateLightTimeCorrection(
                perturbedState, receiverState, transmissionTime, receptionTime );
        }
        positionPartial( i ) = ( upPerturbedCorrection - downPerturbedCorrection ) / ( 2.0 * positionPerturbation );
    }
    return positionPartial;
}

}  // namespace observation_models

}  // namespace tudat
