/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/dualFrequencyDopplerCalibration.h"
#include "tudat/simulation/estimation_setup/createLightTimeCalculator.h"

namespace tudat
{
namespace observation_models
{

Eigen::Matrix< long double, Eigen::Dynamic, 6 > getDsnDopplerFrequencyCalibration(
        const LinkEnds& linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector< Time >& receptionUtcTimes,
        const std::vector< long double >& observations,
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings,
        const std::shared_ptr< ground_stations::StationFrequencyInterpolator >& transmittingFrequency,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& nonDispersiveCorrections,
        const std::shared_ptr< LightTimeCorrectionSettings >& ionosphericCorrection )
{
    if( linkEnds.size( ) != 3 )
    {
        throw std::runtime_error( "Error in DSN dual-frequency calibration, only two-way links are supported" );
    }

    // Retrieve turnaround ratios, count interval and reference frequency
    std::vector< FrequencyBands > frequencyBands =
            convertDoubleVectorToFrequencyBands( ancillarySettings->getAncillaryDoubleVectorData( frequency_bands ) );
    FrequencyBands referenceBand =
            convertDoubleToFrequencyBand( ancillarySettings->getAncillaryDoubleData( reception_reference_frequency_band ) );
    std::function< double( FrequencyBands, FrequencyBands ) > turnaroundRatio =
            bodies.at( linkEnds.at( reflector1 ).bodyName_ )->getVehicleSystems( )->getTransponderTurnaroundRatio( );
    long double currentTurnaroundRatio = turnaroundRatio( frequencyBands.at( 0 ), frequencyBands.at( 1 ) );
    long double referenceCarrier = static_cast< long double >( turnaroundRatio( referenceBand, frequencyBands.at( 1 ) ) ) *
            static_cast< long double >( ancillarySettings->getAncillaryDoubleData( doppler_reference_frequency ) );
    double integrationTime = ancillarySettings->getAncillaryDoubleData( doppler_integration_time );

    // Create light-time models
    std::shared_ptr< FullLinkLightTimeCalculator< long double, Time > > lightTimeCalculator =
            createFullLinkLightTimeCalculator< long double, Time >(
                    linkEnds, bodies, dsn_n_way_averaged_doppler, { nonDispersiveCorrections, nonDispersiveCorrections } );
    std::vector< std::shared_ptr< LightTimeCorrection > > ionosphericCorrections = {
        createLightTimeCorrections( ionosphericCorrection, bodies, linkEnds, transmitter, reflector1, dsn_n_way_averaged_doppler ),
        createLightTimeCorrections( ionosphericCorrection, bodies, linkEnds, reflector1, receiver, dsn_n_way_averaged_doppler )
    };

    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeConverter = earth_orientation::createDefaultTimeConverter( );
    auto getStationPosition = [ & ]( const LinkEndType linkEndType ) {
        return bodies.at( linkEnds.at( linkEndType ).bodyName_ )
                ->getGroundStation( linkEnds.at( linkEndType ).getReferencePointName( ) )
                ->getNominalStationState( )
                ->getNominalCartesianPosition( );
    };
    Eigen::Vector3d transmitterPosition = getStationPosition( transmitter );
    Eigen::Vector3d receiverPosition = getStationPosition( receiver );

    Eigen::Matrix< long double, Eigen::Dynamic, 6 > calibration( receptionUtcTimes.size( ), 6 );
    for( unsigned int i = 0; i < receptionUtcTimes.size( ); i++ )
    {
        // Compute emission time, transmitted frequency and ionospheric phase delay at start and end of count interval
        Time emissionUtcTimes[ 2 ];
        long double transmittedFrequencies[ 2 ], ionosphericPhaseDelays[ 2 ];
        for( unsigned int j = 0; j < 2; j++ )
        {
            Time receptionTdbTime = timeConverter->getCurrentTime( basic_astrodynamics::utc_scale,
                                                                   basic_astrodynamics::tdb_scale,
                                                                   receptionUtcTimes.at( i ) + ( j == 0 ? -0.5 : 0.5 ) * integrationTime,
                                                                   receiverPosition );
            std::vector< double > linkEndTimes;
            std::vector< Eigen::Vector6d > linkEndStates;
            long double lightTime = lightTimeCalculator->calculateLightTimeWithLinkEndsStates(
                    receptionTdbTime, receiver, linkEndTimes, linkEndStates, ancillarySettings );
            emissionUtcTimes[ j ] = timeConverter->getCurrentTime(
                    basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, receptionTdbTime - lightTime, transmitterPosition );

            transmittedFrequencies[ j ] = transmittingFrequency->getTemplatedCurrentFrequency< long double, Time >( emissionUtcTimes[ j ] );
            ancillarySettings->setIntermediateDoubleData( transmitter_frequency_intermediate, transmittedFrequencies[ j ] );
            ancillarySettings->setIntermediateDoubleData( received_frequency_intermediate,
                                                          currentTurnaroundRatio * transmittedFrequencies[ j ] );
            double ionosphericLightTime = 0.0;
            for( unsigned int leg = 0; leg < 2; leg++ )
            {
                ionosphericLightTime += ionosphericCorrections.at( leg )->calculateLightTimeCorrectionWithMultiLegLinkEndStates(
                        linkEndStates, linkEndTimes, 2 * leg, ancillarySettings );
            }
            ionosphericPhaseDelays[ j ] = transmittedFrequencies[ j ] * ionosphericLightTime;
        }

        long double emissionDuration = static_cast< long double >( emissionUtcTimes[ 1 ] - emissionUtcTimes[ 0 ] );
        long double meanTransmittedFrequency =
                transmittingFrequency->getTemplatedFrequencyIntegral< long double, Time >( emissionUtcTimes[ 0 ], emissionUtcTimes[ 1 ] ) /
                emissionDuration;
        long double carrierNormalization = currentTurnaroundRatio * meanTransmittedFrequency;

        // Received frequency (from DSN Doppler sign convention), normalized by carrier
        calibration( i, 0 ) = ( referenceCarrier - observations.at( i ) ) / carrierNormalization;
        calibration( i, 1 ) = carrierNormalization;
        // Change in ionospheric phase delay over count reduces the normalized received frequency
        calibration( i, 2 ) =
                -( ionosphericPhaseDelays[ 1 ] - ionosphericPhaseDelays[ 0 ] ) / ( meanTransmittedFrequency * integrationTime );
        calibration( i, 3 ) = std::pow( 1.0E9L / meanTransmittedFrequency, 2 );
        calibration( i, 4 ) = std::pow( 1.0E9L / carrierNormalization, 2 );
        calibration( i, 5 ) = emissionDuration;
    }
    return calibration;
}

}  // namespace observation_models
}  // namespace tudat
