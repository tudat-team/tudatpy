/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Relativity/relativisticLightTimeCorrection.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/firstOrderRelativisticLightTimeCorrection.h"


namespace tudat
{

namespace observation_models
{

//! Function to calculate first order relativistic light time correction due to set of gravitating point masses.
double FirstOrderLightTimeCorrectionCalculator::calculateLightTimeCorrection(
        const Eigen::Vector6d& transmitterState,
        const Eigen::Vector6d& receiverState,
        const double transmissionTime,
        const double receptionTime )
{
    // Retrieve ppn parameter gamma.
    double ppnParameterGamma = ppnParameterGammaFunction_( );

    // Initialize correction to zero.
    currentTotalLightTimeCorrection_ = 0.0;

    double evaluationTime = TUDAT_NAN;
    // Iterate over all gravitating bodies.
    for( unsigned int i = 0; i < perturbingBodyStateFunctions_.size( ); i++ )
    {
        evaluationTime = transmissionTime + lightTimeEvaluationContribution_.at( i ) * ( receptionTime - transmissionTime );
        // Calculate correction due to current body and add to total.
        currentLighTimeCorrectionComponents_[ i ] = relativity::calculateFirstOrderLightTimeCorrectionFromCentralBody(
                    perturbingBodyGravitationalParameterFunctions_[ i ]( ),
                    transmitterState.segment( 0, 3 ), receiverState.segment( 0, 3 ),
                    perturbingBodyStateFunctions_[ i ]( evaluationTime ).segment( 0, 3 ),
                    ppnParameterGamma );
        currentTotalLightTimeCorrection_ += currentLighTimeCorrectionComponents_[ i ];
    }

    return currentTotalLightTimeCorrection_;
}

//! Function to compute the partial derivative of the light-time correction w.r.t. link end position
Eigen::Matrix< double, 3, 1 > FirstOrderLightTimeCorrectionCalculator::
calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
        const Eigen::Vector6d& transmitterState,
        const Eigen::Vector6d& receiverState,
        const double transmissionTime,
        const double receptionTime,
        const LinkEndType linkEndAtWhichPartialIsEvaluated )
{
    // Retrieve ppn parameter gamma.
    double ppnParameterGamma = ppnParameterGammaFunction_( );

    // Initialize correction to zero.
    Eigen::Matrix< double, 3, 1 > currentTotalLightTimeCorrectionPartial_ = Eigen::Matrix< double, 3, 1 >::Zero( );

    double evaluationTime = TUDAT_NAN;

    Eigen::Vector6d perturbingBodyState;
    // Iterate over all gravitating bodies.
    for( unsigned int i = 0; i < perturbingBodyStateFunctions_.size( ); i++ )
    {
        evaluationTime = transmissionTime + lightTimeEvaluationContribution_.at( i ) * ( receptionTime - transmissionTime );

        perturbingBodyState = perturbingBodyStateFunctions_[ i ]( evaluationTime );

        // Calculate correction due to current body and add to total.
        currentTotalLightTimeCorrectionPartial_ += relativity::calculateFirstOrderCentralBodyLightTimeCorrectionGradient(
                    perturbingBodyGravitationalParameterFunctions_[ i ]( ),
                    transmitterState.segment( 0, 3 ), receiverState.segment( 0, 3 ),
                    perturbingBodyStateFunctions_[ i ]( evaluationTime ).segment( 0, 3 ),
                ( linkEndAtWhichPartialIsEvaluated == receiver ),
                ppnParameterGamma );
    }

    return currentTotalLightTimeCorrectionPartial_;
}

}

}

