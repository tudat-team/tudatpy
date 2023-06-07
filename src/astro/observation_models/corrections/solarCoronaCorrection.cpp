/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/observation_models/corrections/solarCoronaCorrection.h"

#include "tudat/math/basic/linearAlgebra.h"

namespace tudat
{

namespace observation_models
{

double SolarCoronaCorrection::computeMinimumDistanceOfLineOfSight(
        Eigen::Vector3d transmitterPositionWrtSun,
        Eigen::Vector3d receiverPositionWrtSun )
{

    // Moyer (2000), eq. 10-58
    Eigen::Vector3d lineOfSight = ( receiverPositionWrtSun - transmitterPositionWrtSun ).normalized( );

    // Moyer (2000), eq. 10-59
    Eigen::Vector3d minimumDistancePoint =
            transmitterPositionWrtSun - transmitterPositionWrtSun.dot( lineOfSight ) * lineOfSight;

    // Check if point is valid. If not, set minimum distance to NAN
    // Moyer (2000), eqs. 10-67, 10-68
    if ( transmitterPositionWrtSun.dot( lineOfSight ) < 0 && receiverPositionWrtSun.dot( lineOfSight ) > 0 )
    {
        return minimumDistancePoint.norm( );
    }
    else
    {
        return TUDAT_NAN;
    }

}

double SolarCoronaCorrection::getCurrentFrequency(
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings,
        const double currentTime )
{
    // Retrieve frequency bands
    std::vector< FrequencyBands > frequencyBands;
    if( ancillarySettings == nullptr )
    {
        throw std::runtime_error(
                "Error when computing solar corona corrections: no ancillary settings found. " );
    }
    try
    {
        frequencyBands = convertDoubleVectorToFrequencyBands( ancillarySettings->getAncilliaryDoubleVectorData( frequency_bands ) );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error(
                "Error when retrieving frequency bands for solar corona corrections: " +
                std::string( caughtException.what( ) ) );
    }

    return transmittedFrequencyFunction_( frequencyBands, currentTime );
}

double InversePowerSeriesSolarCoronaCorrection::calculateLightTimeCorrectionWithMultiLegLinkEndStates(
        const std::vector< Eigen::Vector6d >& linkEndsStates,
        const std::vector< double >& linkEndsTimes,
        const unsigned int currentMultiLegTransmitterIndex,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings )
{
    // Retrieve state and time of receiver and transmitter
    Eigen::Vector6d legTransmitterState, legReceiverState;
    double legTransmissionTime, legReceptionTime;
    getTransmissionReceptionTimesAndStates(
            linkEndsStates, linkEndsTimes, currentMultiLegTransmitterIndex, legTransmitterState, legReceiverState,
            legTransmissionTime, legReceptionTime );

    // Retrieve state wrt Sun
    Eigen::Vector3d transmitterPositionWrtSun = ( legTransmitterState - sunStateFunction_( legTransmissionTime ) ).segment( 0, 3 );
    Eigen::Vector3d receiverPositionWrtSun = ( legReceiverState - sunStateFunction_( legReceptionTime ) ).segment( 0, 3 );

    // Check if minimum distance along line of sight is valid
    if ( std::isnan( computeMinimumDistanceOfLineOfSight( transmitterPositionWrtSun, receiverPositionWrtSun ) ) )
    {
        throw std::runtime_error( "Error when computing solar corona correction: minimum distance along LOS is not valid." );
    }

    double sunReceiverTransmitterAngle = linear_algebra::computeAngleBetweenVectors(
            - receiverPositionWrtSun,
            transmitterPositionWrtSun - receiverPositionWrtSun );
    double receiverSunTransmitterAngle = linear_algebra::computeAngleBetweenVectors(
            receiverPositionWrtSun, transmitterPositionWrtSun );

    // Verma et al. (2013), eq. 3/4
    double electronDensityIntegral = 0.0;
    for ( unsigned int i = 0; i < coefficients_.size( ); ++i )
    {
        // If exponent is integer
        if ( std::fmod( positiveExponents_.at( i ), 1.0 ) == 0 )
        {
            electronDensityIntegral += coefficients_.at( i ) * computeSingleTermIntegralAnalytically(
                    receiverPositionWrtSun, sunReceiverTransmitterAngle, receiverSunTransmitterAngle,
                    static_cast< unsigned int >( positiveExponents_.at( i ) ) );
        }
        else
        {
            electronDensityIntegral += coefficients_.at( i ) * computeSingleTermIntegralViaTaylorExpansion(
                    receiverPositionWrtSun, sunReceiverTransmitterAngle, receiverSunTransmitterAngle,
                    positiveExponents_.at( i ) );
        }
    }

    // Verma et al. (2013), eq. 1
    return sign_ * criticalPlasmaDensityDelayCoefficient_ * std::pow( getCurrentFrequency( ancillarySettings, linkEndsTimes.front( ) ), 2.0 ) *
        electronDensityIntegral / physical_constants::getSpeedOfLight< double >( );
}

// Verma et al. (2006), eq. A.6
double InversePowerSeriesSolarCoronaCorrection::computeSingleTermIntegralViaTaylorExpansion(
        const Eigen::Vector3d& receiverPositionWrtSun,
        const double sunReceiverTransmitterAngle,
        const double receiverSunTransmitterAngle,
        const double positiveExponent )
{
    // Rename terms for simplicity
    double alpha = sunReceiverTransmitterAngle;
    double beta = receiverSunTransmitterAngle;

    double mulTerm = std::pow( sunRadius_, positiveExponent ) /
            std::pow( receiverPositionWrtSun.norm( ) * std::sin( alpha ), positiveExponent - 1.0 );

    double sumTerm1 = beta;

    double sumTerm2 = - ( positiveExponent - 2.0 ) / 6.0 * (
            std::pow( beta + alpha - mathematical_constants::PI / 2.0, 3.0 ) -
            std::pow( alpha - mathematical_constants::PI / 2.0, 3.0 ) );

    double sumTerm3 = ( 3.0 * std::pow( positiveExponent, 2.0 ) - 14.0 * positiveExponent + 16.0 ) / 12.0 *
            ( std::pow( beta + alpha - mathematical_constants::PI / 2.0, 5.0 ) -  std::pow( alpha - mathematical_constants::PI / 2.0, 5.0 ) );

    return mulTerm * ( sumTerm1 + sumTerm2 + sumTerm3 );
}

// Verma et al., eq. A.4
double InversePowerSeriesSolarCoronaCorrection::computeSingleTermIntegralAnalytically(
        const Eigen::Vector3d& receiverPositionWrtSun,
        const double sunReceiverTransmitterAngle,
        const double receiverSunTransmitterAngle,
        const unsigned int positiveExponent )
{
    // Rename terms for simplicity
    double alpha = sunReceiverTransmitterAngle;
    double beta = receiverSunTransmitterAngle;

    double upperBound = beta + alpha - mathematical_constants::PI / 2.0;
    double lowerBound = alpha - mathematical_constants::PI / 2.0;

    double mulTerm = std::pow( sunRadius_, positiveExponent ) /
            std::pow( receiverPositionWrtSun.norm( ) * std::sin( alpha ), positiveExponent - 1.0 );

    return mulTerm * computeCosinePowerIntegral( lowerBound, upperBound, positiveExponent );
}

// Wikipedia, https://en.wikipedia.org/wiki/List_of_integrals_of_trigonometric_functions
double InversePowerSeriesSolarCoronaCorrection::computeCosinePowerIntegral(
        const double lowerBound,
        const double upperBound,
        const unsigned int positiveExponent )
{
    double integral;

    // If integral hasn't been computed yet
    if ( cosinePowersIntegrals_.count( positiveExponent ) == 0 )
    {
        if ( positiveExponent == 0 )
        {
            integral = upperBound - lowerBound;
        }
        else if ( positiveExponent == 1 )
        {
            integral = std::sin( upperBound ) - std::sin( lowerBound );
        }
        else
        {
            std::function< double ( double ) > indefiniteIntegral = [=] ( const double x ) {
                return ( std::pow( std::cos( x ), positiveExponent - 1.0 ) * std::sin( x ) ) / positiveExponent; };

            double term1 = indefiniteIntegral( upperBound ) - indefiniteIntegral( lowerBound );

            double term2 = static_cast< double >( positiveExponent - 1 ) / positiveExponent * computeCosinePowerIntegral(
                    lowerBound, upperBound, positiveExponent - 2 );

            integral = term1 + term2;
        }

        // Set computed value of integral
        cosinePowersIntegrals_[ positiveExponent ] = integral;
    }
    // If integral has already been computed, use previous value
    else
    {
        integral = cosinePowersIntegrals_.at( positiveExponent );
    }

    return integral;
}


} // namespace observation_models

} // namespace tudat