#include <iostream>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Astrodynamics/Relativity/relativisticLightTimeCorrection.h"
#include "Astrodynamics/ObservationModels/ObservableCorrections/firstOrderRelativisticLightTimeCorrection.h"


namespace tudat
{

namespace observation_models
{

//! Function to calculate first order relativistic light time correction due to set of gravitating point masses.
double FirstOrderLightTimeCorrectionCalculator::calculateLightTimeCorrection(
        const basic_mathematics::Vector6d& transmitterState,
        const basic_mathematics::Vector6d& receiverState,
        const double transmissionTime,
        const double receptionTime )
{
    // Retrieve ppn parameter gamma.
    double ppnParameterGamma = ppnParameterGammaFunction_( );

    // Initialize correction to zero.
    currentTotalLightTimeCorrection_ = 0.0;

    // Iterate over all gravitating bodies.
    for( unsigned int i = 0; i < perturbingBodyStateFunctions_.size( ); i++ )
    {
        // Calculate correction due to current body and add to total.
        currentLighTimeCorrectionComponents_[ i ] = relativity::calculateFirstOrderLightTimeCorrectionFromCentralBody(
                    perturbingBodyGravitationalParameterFunctions_[ i ]( ),
                    transmitterState.segment( 0, 3 ), receiverState.segment( 0, 3 ),
                    perturbingBodyStateFunctions_[ i ]( ( transmissionTime + receptionTime ) / 2.0 ).segment( 0, 3 ),
                    ppnParameterGamma );
        currentTotalLightTimeCorrection_ += currentLighTimeCorrectionComponents_[ i ];
    }

    return currentTotalLightTimeCorrection_;
}

double FirstOrderLightTimeCorrectionDerivativeCalculator::calculateLightTimeDerivativeCorrection(
        const basic_mathematics::Vector6d& transmitterState,
        const basic_mathematics::Vector6d& receiverState,
        const double transmissionTime,
        const double receptionTime,
        const bool fixTransmissionTime )
{
    // Retrieve ppn parameter gamma.
    double ppnParameterGamma = ppnParameterGammaFunction_( );

    // Initialize correction to zero.
    currentTotalLightTimeDerivativeCorrection_ = 0.0;

    // Iterate over all gravitating bodies.
    for( unsigned int i = 0; i < perturbingBodyStateFunctions_.size( ); i++ )
    {
        // Calculate correction due to current body and add to total.
        currentLighTimeCorrectionDerivativeComponents_[ i ] = relativity::calculateFirstOrderLightTimeCorrectionDerivativeFromCentralBody(
                    perturbingBodyGravitationalParameterFunctions_[ i ]( ),
                    transmitterState, receiverState,
                    perturbingBodyStateFunctions_[ i ]( ( transmissionTime + receptionTime ) / 2.0 ),
                    ppnParameterGamma,
                    fixTransmissionTime );
        currentTotalLightTimeDerivativeCorrection_ += currentLighTimeCorrectionDerivativeComponents_[ i ];
    }

    return currentTotalLightTimeDerivativeCorrection_;
}

}

}

