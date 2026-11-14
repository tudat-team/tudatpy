/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATE_OBSERVATION_COLLECTION_H
#define TUDAT_CREATE_OBSERVATION_COLLECTION_H

#include <Eigen/Core>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedRotationalEphemeris.h"
#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/basics/timeType.h"
#include "tudat/basics/utilities.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/trackingSupplementaryData.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/trackingDataWeighting.h"
#include "tudat/simulation/estimation_setup/singleObservationSet.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/simulation/estimation_setup/observationCollectionConversion.h"

namespace tudat
{

namespace observation_models
{

observation_models::ObservableType getObservableTypeFromTrackingDataString( const std::string& observableTypeString );

observation_models::LinkEnds getLinkEndsFromTrackingData(
        const std::vector< std::pair< std::pair< std::string, std::string >, std::string > >& rawLinkEnds );

void checkTrackingDataLinkEnds( const observation_models::ObservableType observableType,
                                const observation_models::LinkEnds& linkEnds,
                                const observation_models::LinkEndType referenceLinkEnd );

bool shouldSkipObservationCollectionAncillarySetting( const std::string& ancillarySetting );

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > getAncillarySettingsFromTrackingData(
        const std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > trackingData )
{
    // Create empty ancillary simulation settings
    std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings =
            std::make_shared< observation_models::ObservationAncillarySimulationSettings >( );

    // Parse and add ancillary settings of type scalar doubles
    for( auto& it : trackingData->getAncillarySettingsDouble( ) )
    {
        if( trackingData->isObservationMetadata( it.first ) || shouldSkipObservationCollectionAncillarySetting( it.first ) )
        {
            continue;
        }
        ancillarySettings->setAncillaryDoubleData( ancillarySettings->getAncillaryVariableFromString( it.first ), it.second );
    }

    // Parse and add ancillary settings of type double vectors
    for( auto& it : trackingData->getAncillarySettingsDoubleVector( ) )
    {
        if( trackingData->isObservationMetadata( it.first ) || shouldSkipObservationCollectionAncillarySetting( it.first ) )
        {
            continue;
        }
        ancillarySettings->setAncillaryDoubleVectorData( ancillarySettings->getAncillaryVariableFromString( it.first ), it.second );
    }

    // Parse and add ancillary settings of type string vectors (frequency band(s))
    for( auto& it : trackingData->getAncillarySettingsStringVector( ) )
    {
        if( trackingData->isObservationMetadata( it.first ) || shouldSkipObservationCollectionAncillarySetting( it.first ) )
        {
            continue;
        }
        observation_models::ObservationAncillarySimulationVariable ancillaryVariable =
                ancillarySettings->getAncillaryVariableFromString( it.first );

        // Check consistency between ancillary variable and value type
        if( ancillaryVariable != observation_models::frequency_bands )
        {
            throw std::runtime_error(
                    "Error when setting ancillary simulation settings from tracking data, inconsistency between the ancillary variable (" +
                    it.first + ") and the type of ancillary value provided (string vector)." );
        }

        // Convert frequency bands (strings) to doubles
        std::vector< double > bandDoubles;
        for( const std::string& band : it.second )
        {
            bandDoubles.push_back(
                    basic_astrodynamics::convertFrequencyBandToDouble( basic_astrodynamics::getFrequencyBandFromString( band ) ) );
        }

        // Set frequency bands as double vector
        ancillarySettings->setAncillaryDoubleVectorData( ancillaryVariable, bandDoubles );
    }

    // Parse and add ancillary settings of type string (reception reference frequency band)
    for( auto& it : trackingData->getAncillarySettingsString( ) )
    {
        if( trackingData->isObservationMetadata( it.first ) || shouldSkipObservationCollectionAncillarySetting( it.first ) )
        {
            continue;
        }
        observation_models::ObservationAncillarySimulationVariable ancillaryVariable =
                ancillarySettings->getAncillaryVariableFromString( it.first );

        // Check consistency between ancillary variable and value type
        if( ancillaryVariable != observation_models::reception_reference_frequency_band )
        {
            throw std::runtime_error(
                    "Error when setting ancillary simulation settings from tracking data, inconsistency between the ancillary variable (" +
                    it.first + ") and the type of ancillary value provided (string)." );
        }

        // Convert reference frequency band (string) to double and set it
        ancillarySettings->setAncillaryDoubleData(
                ancillaryVariable,
                basic_astrodynamics::convertFrequencyBandToDouble( basic_astrodynamics::getFrequencyBandFromString( it.second ) ) );
    }

    // Metadata and reader-only settings may fill the input maps without producing
    // simulation settings. Angular observations require nullptr in that case.
    if( ancillarySettings->getDoubleData( ).empty( ) && ancillarySettings->getDoubleVectorData( ).empty( ) )
    {
        return nullptr;
    }

    return ancillarySettings;
}

//! Estimation-side evaluator; tracking data stores only plain correction settings.
template< typename ObservationScalarType, typename TimeType >
using TrackingDataCorrectionFunction = std::function< std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >(
        const std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > >&,
        const simulation_setup::SystemOfBodies&,
        const std::vector< TimeType >& ) >;

//! Add the observations in one tracking-data object as one dataset set.
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
int addTrackingDataToObservationDataset( const std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > trackingData,
                                         const simulation_setup::SystemOfBodies& bodies,
                                         ObservationDataset< ObservationScalarType, TimeType >& observationDataset,
                                         const bool applyCorrections = false,
                                         const TrackingDataCorrectionFunction< ObservationScalarType, TimeType >& correctionFunction = {} )
{
    // Identify observable type from tracking data object
    observation_models::ObservableType observableType = getObservableTypeFromTrackingDataString( trackingData->getObservableType( ) );

    // Get reference link end from tracking data object
    observation_models::LinkEndType referenceLinkEnd = observation_models::getLinkEndTypeFromString( trackingData->getReferenceLinkEnd( ) );

    // Identify link ends from tracking data object and check validity
    observation_models::LinkEnds rawLinkEnds = getLinkEndsFromTrackingData( trackingData->getLinkEnds( ) );
    checkTrackingDataLinkEnds( observableType, rawLinkEnds, referenceLinkEnd );
    LinkDefinition linkEnds = LinkDefinition( rawLinkEnds );

    // Get observations from tracking data
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations = trackingData->getObservations( );

    // Apply corrections if requested (and if they exist)
    if( applyCorrections && !trackingData->getCorrectionSettings( )->getObservationCorrections( ).empty( ) )
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > corrections =
                trackingData->getCorrectionSettings( )->getObservationCorrections( );

        // Check size consistency
        if( corrections.size( ) != observations.size( ) )
        {
            throw std::runtime_error( "Error when creating single observation set from tracking data, the size of the corrections (" +
                                      std::to_string( corrections.size( ) ) + ") is inconsistent with the number of observations (" +
                                      std::to_string( observations.size( ) ) + ")." );
        }

        // Apply corrections
        for( unsigned int i = 0; i < observations.size( ); i++ )
        {
            // Check size consistency of single correction
            if( corrections[ i ].size( ) != observations[ i ].size( ) )
            {
                throw std::runtime_error(
                        "Error when creating single observation set from tracking data, size of single observation "
                        "correction (" +
                        std::to_string( corrections[ i ].size( ) ) + ") does not match the single observation size (" +
                        std::to_string( observations[ i ].size( ) ) + ")." );
            }
            observations[ i ] += corrections[ i ];
        }
    }

    auto epochsInput = trackingData->getObservationEpochs( );
    auto epochsTdb = std::vector< TimeType >( epochsInput.size( ) );
    std::string referenceLinkEndName = trackingData->getReferencePointName( );

    if( trackingData->getTimeScale( ) != "TDB" )
    {
        auto timeScaleConverter = earth_orientation::createDefaultTimeConverter( );
        auto const& inputScale = basic_astrodynamics::timeScaleFromString( trackingData->getTimeScale( ) );
        Eigen::Vector3d timeScaleConversionPosition = Eigen::Vector3d::Zero( );
        const LinkEndId referenceLinkEndId = rawLinkEnds.at( referenceLinkEnd );
        if( referenceLinkEndId.bodyName_ == "Earth" && bodies.doesBodyExist( referenceLinkEndId.bodyName_ ) )
        {
            const std::shared_ptr< simulation_setup::Body > earthBody = bodies.getBody( referenceLinkEndId.bodyName_ );
            const std::map< std::string, std::shared_ptr< ground_stations::GroundStation > > groundStations =
                    earthBody->getGroundStationMap( );
            auto groundStationIterator = groundStations.find( referenceLinkEndName );
            if( groundStationIterator != groundStations.end( ) && groundStationIterator->second != nullptr )
            {
                timeScaleConversionPosition = groundStationIterator->second->getNominalStationState( )->getNominalCartesianPosition( );
            }
        }

        epochsTdb = timeScaleConverter->getCurrentTimesFromSinglePosition< TimeType >(
                inputScale, basic_astrodynamics::tdb_scale, epochsInput, timeScaleConversionPosition );
    }
    else
    {
        epochsTdb = epochsInput;
    }

    // Time biases and uncertainties are currently only supported for angular position observations, for which they are converted
    // to corrections and uncertainties of the observations using the rates of right ascension and declination
    const std::shared_ptr< data::TimingBasedCorrectionSettings > timeBiases =
            trackingData->getCorrectionSettings( )->getTimingBasedCorrectionSettings( );
    const std::shared_ptr< data::TimingBasedUncertaintySettings > timeUncertainties =
            trackingData->getWeightSettings( )->getTimingBasedUncertaintySettings( );
    if( ( timeBiases != nullptr || timeUncertainties != nullptr ) && observableType != observation_models::angular_position )
    {
        throw std::runtime_error(
                "Error when adding tracking data to an observation dataset, observation time biases and uncertainties are only "
                "supported for angular position observations." );
    }
    std::vector< Eigen::Vector2d > angularRates;
    if( ( applyCorrections && timeBiases != nullptr ) || timeUncertainties != nullptr )
    {
        angularRates = getRightAscensionAndDeclinationRates< TimeType >( bodies, rawLinkEnds, epochsTdb );
    }
    if( applyCorrections && timeBiases != nullptr )
    {
        for( unsigned int i = 0; i < timeBiases->timeBiases_.size( ); ++i )
        {
            observations[ i ] -= ( timeBiases->timeBiases_.at( i ) * angularRates.at( i ) ).template cast< ObservationScalarType >( );
        }
    }

    // Evaluate environment-dependent corrections at TDB reception times using the original observations.
    if( applyCorrections && trackingData->getCorrectionSettings( )->hasObservationCorrectionSettings( ) )
    {
        if( !correctionFunction )
        {
            throw std::runtime_error( "Tracking data requests environment-dependent angular corrections, but no evaluator was supplied." );
        }
        const auto corrections = correctionFunction( trackingData, bodies, epochsTdb );
        if( corrections.size( ) != observations.size( ) )
        {
            throw std::runtime_error( "Computed angular corrections have an inconsistent number of observations." );
        }
        for( unsigned int i = 0; i < observations.size( ); ++i )
        {
            if( corrections[ i ].size( ) != observations[ i ].size( ) )
            {
                throw std::runtime_error( "Computed angular corrections have an inconsistent observable size." );
            }
            observations[ i ] += corrections[ i ];
            const ObservationScalarType pi = std::acos( static_cast< ObservationScalarType >( -1 ) );
            observations[ i ]( 0 ) = std::fmod( observations[ i ]( 0 ) + pi, 2 * pi );
            if( observations[ i ]( 0 ) < 0 )
            {
                observations[ i ]( 0 ) += 2 * pi;
            }
            observations[ i ]( 0 ) -= pi;
        }
    }

    // Convert ancillary settings information from tracking data object to ObservationAncillarySimulationSettings
    std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings =
            getAncillarySettingsFromTrackingData< ObservationScalarType, TimeType >( trackingData );

    // Compute the weights from the weight settings of the tracking data
    const std::shared_ptr< ObservationWeightSettings > observationWeightSettings =
            computeObservationWeights( trackingData, bodies, rawLinkEnds, referenceLinkEnd, angularRates );

    int setId;
    if( observationWeightSettings != nullptr )
    {
        setId = observationDataset.addObservationSetWithWeights(
                observableType,
                linkEnds,
                observations,
                epochsTdb,
                referenceLinkEnd,
                observationWeightSettings,
                std::vector< Eigen::VectorXd >( ),
                nullptr,
                ancillarySettings,
                std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
                true );
    }
    else
    {
        setId = observationDataset.addObservationSet( observableType,
                                                      linkEnds,
                                                      observations,
                                                      epochsTdb,
                                                      referenceLinkEnd,
                                                      std::vector< Eigen::VectorXd >( ),
                                                      nullptr,
                                                      ancillarySettings,
                                                      std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                                                      std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
                                                      true );
    }

    // Match the stable time ordering used when the observation set was created.
    if( !trackingData->getNumericalObservationMetadata( ).empty( ) )
    {
        std::vector< std::size_t > permutation( epochsTdb.size( ) );
        for( std::size_t i = 0; i < permutation.size( ); ++i )
        {
            permutation[ i ] = i;
        }
        std::stable_sort( permutation.begin( ), permutation.end( ), [ &epochsTdb ]( std::size_t i, std::size_t j ) {
            return epochsTdb[ i ] < epochsTdb[ j ];
        } );
        const auto& observationIds = observationDataset.getObservationIdsForSet( setId );
        for( const auto& metadata : trackingData->getNumericalObservationMetadata( ) )
        {
            std::vector< double > values;
            values.reserve( permutation.size( ) );
            for( const auto index : permutation )
            {
                values.push_back( metadata.second.at( index ) );
            }
            observationDataset.addNumericalObservationMetadata( metadata.first, observationIds, values );
        }
    }
    return setId;
}

// Create single observation set object from tracking data object
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > createSingleObservationSetFromTrackingData(
        const std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > trackingData,
        const SystemOfBodies& bodies,
        const bool applyCorrections = false )
{
    auto dataset = std::make_shared< ObservationDataset< ObservationScalarType, TimeType > >( );
    const int setId = addTrackingDataToObservationDataset( trackingData, bodies, *dataset, applyCorrections );
    return createSingleObservationSet( dataset, setId );
}

//! Create a dataset from one or more generic tracking-data objects.
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createObservationDatasetFromTrackingData(
        const std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > > trackingDataList,
        simulation_setup::SystemOfBodies& bodies,
        const bool applyCorrections = false,
        const TrackingDataCorrectionFunction< ObservationScalarType, TimeType >& correctionFunction = {} )
{
    auto observationDataset = std::make_shared< ObservationDataset< ObservationScalarType, TimeType > >( );
    for( const auto& trackingData : trackingDataList )
    {
        addTrackingDataToObservationDataset( trackingData, bodies, *observationDataset, applyCorrections, correctionFunction );
    }
    return observationDataset;
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > createObservationCollection(
        const std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > > trackingDataList,
        SystemOfBodies& bodies,
        const bool applyCorrections = false )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            createObservationDatasetFromTrackingData< ObservationScalarType, TimeType >( trackingDataList, bodies, applyCorrections ) );
}

template< typename EphemerisScalarType, typename EphemerisTimeType >
inline void resetTabulatedEphemerisFromTrackingSupplementaryStateHistory(
        const std::map< double, Eigen::Vector6d >& stateHistory,
        const std::shared_ptr< ephemerides::TabulatedCartesianEphemeris< EphemerisScalarType, EphemerisTimeType > > tabulatedEphemeris )
{
    std::map< EphemerisTimeType, Eigen::Matrix< EphemerisScalarType, 6, 1 > > castStateHistory;
    utilities::castMatrixMap< double, double, EphemerisTimeType, EphemerisScalarType, 6, 1 >( stateHistory, castStateHistory );

    tabulatedEphemeris->resetInterpolator( interpolators::createOneDimensionalInterpolator(
            castStateHistory,
            interpolators::linearInterpolation( interpolators::huntingAlgorithm, interpolators::extrapolate_at_boundary ) ) );
}
void resetTabulatedEphemerisFromTrackingSupplementaryStateHistory( const std::map< double, Eigen::Vector6d >& stateHistory,
                                                                   const std::shared_ptr< ephemerides::Ephemeris > ephemeris,
                                                                   const std::string& bodyName );

template< typename EphemerisScalarType, typename EphemerisTimeType >
inline void resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory(
        const std::map< double, Eigen::Vector7d >& rotationalStateHistory,
        const std::shared_ptr< ephemerides::TabulatedRotationalEphemeris< EphemerisScalarType, EphemerisTimeType > >
                tabulatedRotationalEphemeris )
{
    std::map< EphemerisTimeType, Eigen::Matrix< EphemerisScalarType, 7, 1 > > castRotationalStateHistory;
    utilities::castMatrixMap< double, double, EphemerisTimeType, EphemerisScalarType, 7, 1 >( rotationalStateHistory,
                                                                                              castRotationalStateHistory );

    tabulatedRotationalEphemeris->reset(
            interpolators::createOneDimensionalInterpolator( castRotationalStateHistory, interpolators::linearInterpolation( ) ) );
}

void resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory(
        const std::map< double, Eigen::Vector7d >& rotationalStateHistory,
        const std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris,
        const std::string& bodyName );
std::map< double, Eigen::Vector6d > getTranslationalStateHistoryWithVelocity(
        const data::TranslationalStateSupplementaryData& translationalStateSupplementaryData );

void setTranslationalStateSupplementaryDataInBodies(
        simulation_setup::SystemOfBodies& bodies,
        const std::map< std::pair< std::string, std::string >, std::vector< data::TranslationalStateSupplementaryData > >&
                translationalStateSupplementaryData );

void setRotationalStateSupplementaryDataInBodies(
        simulation_setup::SystemOfBodies& bodies,
        const std::map< std::pair< std::string, std::string >, std::vector< data::RotationalStateSupplementaryData > >&
                rotationalStateSupplementaryData );

void setFrequencySupplementaryDataInBodies(
        simulation_setup::SystemOfBodies& bodies,
        const std::map< std::pair< std::string, std::string >, std::vector< std::shared_ptr< data::FrequencySupplementaryData > > >&
                frequencySupplementaryData );

void setInstrumentSupplementaryDataInBodies(
        simulation_setup::SystemOfBodies& bodies,
        const std::map< std::pair< std::string, std::string >, std::vector< std::shared_ptr< data::InstrumentSupplementaryData > > >&
                instrumentSupplementaryData );

//! Create a ground station for every reference point position, or check that an existing one is at the same position.
void setReferencePointSupplementaryDataInBodies(
        simulation_setup::SystemOfBodies& bodies,
        const std::map< std::pair< std::string, std::string >, std::vector< std::shared_ptr< data::ReferencePointSupplementaryData > > >&
                referencePointSupplementaryData );

void setTrackingSupplementaryDataInBodies( simulation_setup::SystemOfBodies& bodies,
                                           const std::vector< data::TrackingSupplementaryData >& supplementaryData );

void setTrackingSupplementaryDataInBodies( simulation_setup::SystemOfBodies& bodies,
                                           const std::vector< std::shared_ptr< data::TrackingSupplementaryData > >& supplementaryData );
}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_CREATE_OBSERVATION_COLLECTION_H
