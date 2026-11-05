/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <functional>
#include <memory>

#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/simulation/estimation_setup/singleObservationSet.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/defaultGroundStationSettings.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"

namespace tudat
{
namespace observation_models
{
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > compressDopplerData(
        const std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > originalDopplerData,
        const unsigned int compressionRatio,
        std::map< std::string, Eigen::Vector3d > approximateGroundStationPositions =
                simulation_setup::getCombinedApproximateGroundStationPositions( ) )
{
    if( compressionRatio == 0 )
    {
        throw std::runtime_error( "Error in Doppler data compression, compression ratio must be positive." );
    }

    const double cadenceTolerance = 0.01;
    ObservationScalarType floatingCompressionRatio =
            mathematical_constants::getFloatingInteger< ObservationScalarType >( compressionRatio );

    double currentCompressionTime = originalDopplerData->getAncillarySettings( )->getAncillaryDoubleData(
            observation_models::ObservationAncillarySimulationVariable::doppler_integration_time );

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > originalObservations =
            originalDopplerData->getObservationsReference( );
    std::vector< TimeType > originalObservationTimesTdb = originalDopplerData->getObservationTimesReference( );

    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter =
            earth_orientation::createDefaultTimeConverter( );
    std::string stationName = originalDopplerData->getLinkEnds( ).at( observation_models::LinkEndType::receiver ).getReferencePointName( );
    auto stationPositionIterator = approximateGroundStationPositions.find( stationName );
    if( stationPositionIterator == approximateGroundStationPositions.end( ) )
    {
        throw std::runtime_error(
                "Error in Doppler data compression, could not retrieve approximate station "
                "position for " +
                stationName );
    }
    Eigen::Vector3d stationPosition = stationPositionIterator->second;

    std::vector< TimeType > originalObservationTimesUtc =
            timeScaleConverter->getCurrentTimesFromSinglePosition< TimeType >( basic_astrodynamics::TimeScales::tdb_scale,
                                                                               basic_astrodynamics::TimeScales::utc_scale,
                                                                               originalObservationTimesTdb,
                                                                               stationPosition );

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > compressedObservations;
    std::vector< TimeType > compressedObservationTimesUtc;

    for( unsigned int runStart = 0; runStart < originalObservations.size( ); )
    {
        unsigned int runEnd = runStart + 1;
        while( runEnd < originalObservations.size( ) )
        {
            double currentTimeStep =
                    static_cast< double >( originalObservationTimesUtc.at( runEnd ) - originalObservationTimesUtc.at( runEnd - 1 ) );
            if( std::fabs( currentTimeStep - currentCompressionTime ) <= cadenceTolerance )
            {
                runEnd++;
            }
            else
            {
                break;
            }
        }

        for( unsigned int blockStart = runStart; blockStart + compressionRatio <= runEnd; blockStart += compressionRatio )
        {
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > newObservable = originalObservations.at( blockStart );
            TimeType newTime = originalObservationTimesUtc.at( blockStart );

            for( unsigned int j = 1; j < compressionRatio; j++ )
            {
                newObservable += originalObservations.at( blockStart + j );
                newTime += originalObservationTimesUtc.at( blockStart + j );
            }

            newObservable /= floatingCompressionRatio;
            newTime = newTime / floatingCompressionRatio;

            compressedObservations.push_back( newObservable );
            compressedObservationTimesUtc.push_back( newTime );
        }

        runStart = runEnd;
    }

    std::vector< Eigen::Vector3d > compressedEarthFixedPositions;
    for( unsigned int i = 0; i < compressedObservationTimesUtc.size( ); ++i )
    {
        compressedEarthFixedPositions.push_back( stationPosition );
    }
    std::vector< TimeType > compressedObservationTimesTdb =
            timeScaleConverter->getCurrentTimes< TimeType >( basic_astrodynamics::TimeScales::utc_scale,
                                                             basic_astrodynamics::TimeScales::tdb_scale,
                                                             compressedObservationTimesUtc,
                                                             compressedEarthFixedPositions );

    std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySimulationSettings =
            std::make_shared< ObservationAncillarySimulationSettings >( *( originalDopplerData->getAncillarySettings( ) ) );
    double originalIntegrationTime = ancillarySimulationSettings->getAncillaryDoubleData(
            observation_models::ObservationAncillarySimulationVariable::doppler_integration_time );
    ancillarySimulationSettings->setAncillaryDoubleData(
            observation_models::ObservationAncillarySimulationVariable::doppler_integration_time,
            originalIntegrationTime * static_cast< double >( floatingCompressionRatio ) );
    return std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
            originalDopplerData->getObservableType( ),
            originalDopplerData->getLinkEnds( ),
            compressedObservations,
            compressedObservationTimesTdb,
            originalDopplerData->getReferenceLinkEnd( ),
            std::vector< Eigen::VectorXd >( ),
            originalDopplerData->getDependentVariableBookkeeping( ),
            ancillarySimulationSettings );
}

//! Function to compute the number of Doppler observations to average to obtain a given count time
template< typename ObservationScalarType = double, typename TimeType = double >
unsigned int getDopplerCompressionRatio(
        const std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > originalDopplerData,
        const double targetCountTime )
{
    const double countTimeTolerance = 1.0e-4;
    double currentCountTime = originalDopplerData->getAncillarySettings( )->getAncillaryDoubleData(
            observation_models::ObservationAncillarySimulationVariable::doppler_integration_time );

    unsigned int compressionRatio = static_cast< unsigned int >( std::round( targetCountTime / currentCountTime ) );
    if( compressionRatio == 0 ||
        std::fabs( static_cast< double >( compressionRatio ) * currentCountTime - targetCountTime ) > countTimeTolerance )
    {
        throw std::runtime_error( "Error in Doppler data compression, target count time " + std::to_string( targetCountTime ) +
                                  " s is not an integer multiple of the count time " + std::to_string( currentCountTime ) + " s." );
    }
    return compressionRatio;
}

//! Function to compress Doppler data, with the compression ratio of each (split) observation set provided by a function
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > >
createCompressedDopplerCollectionFromRatioFunction(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > originalDopplerData,
        const std::function< unsigned int( const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > ) >
                compressionRatioFunction,
        const unsigned int minNumberObservations,
        const double maxArcGap,
        std::map< std::string, Eigen::Vector3d > approximateGroundStationPositions )
{
    // Split Doppler observation sets into arcs
    std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > compressedData =
            splitObservationSets( originalDopplerData,
                                  observationSetSplitter( time_interval_splitter, maxArcGap, minNumberObservations ),
                                  observationParser( observation_models::ObservableType::dsn_n_way_averaged_doppler ) );

    std::map< LinkEnds, std::vector< std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > > >
            uncompressedObservationSets =
                    compressedData->getObservationsSets( ).at( observation_models::ObservableType::dsn_n_way_averaged_doppler );

    std::map< LinkEnds, std::vector< unsigned int > > indicesSetsToRemove;
    for( auto const& [ linkEnds, observationSets ] : uncompressedObservationSets )
    {
        for( unsigned int index = 0; index < observationSets.size( ); index++ )
        {
            std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > compressedDataSet =
                    compressDopplerData< ObservationScalarType, TimeType >( observationSets.at( index ),
                                                                            compressionRatioFunction( observationSets.at( index ) ),
                                                                            approximateGroundStationPositions );

            if( compressedDataSet->getObservationTimes( ).size( ) )
            {
                compressedData->replaceSingleObservationSet( compressedDataSet, index );
            }
            else
            {
                indicesSetsToRemove[ linkEnds ].push_back( index );
            }
        }
    }

    // Remove empty compressed sets if any
    compressedData->removeSingleObservationSets(
            { { observation_models::ObservableType::dsn_n_way_averaged_doppler, indicesSetsToRemove } } );

    return compressedData;
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > createCompressedDopplerCollection(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > originalDopplerData,
        const unsigned int compressionRatio,
        const unsigned int minNumberObservations = 10,
        const double maxArcGap = 300.0,
        std::map< std::string, Eigen::Vector3d > approximateGroundStationPositions =
                simulation_setup::getCombinedApproximateGroundStationPositions( ) )
{
    return createCompressedDopplerCollectionFromRatioFunction< ObservationScalarType, TimeType >(
            originalDopplerData,
            [ = ]( const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > ) { return compressionRatio; },
            minNumberObservations,
            maxArcGap,
            approximateGroundStationPositions );
}

//! Function to retrieve the nominal Earth-fixed positions of all ground stations on Earth
inline std::map< std::string, Eigen::Vector3d > getEarthGroundStationPositions( const simulation_setup::SystemOfBodies& bodies )
{
    std::map< std::string, Eigen::Vector3d > groundStationPositions;
    for( const auto& groundStation : bodies.at( "Earth" )->getGroundStationMap( ) )
    {
        groundStationPositions[ groundStation.first ] = groundStation.second->getNominalStationState( )->getNominalCartesianPosition( );
    }
    return groundStationPositions;
}

//! Function to compress Doppler data to a target count time, using the ground station positions in the bodies
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > createCountTimeCompressedDopplerCollection(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > originalDopplerData,
        const simulation_setup::SystemOfBodies& bodies,
        const double targetCountTime,
        const unsigned int minNumberObservations = 10,
        const double maxArcGap = 300.0 )
{
    return createCompressedDopplerCollectionFromRatioFunction< ObservationScalarType, TimeType >(
            originalDopplerData,
            [ = ]( const std::shared_ptr< SingleObservationSet< ObservationScalarType, TimeType > > observationSet ) {
                return getDopplerCompressionRatio< ObservationScalarType, TimeType >( observationSet, targetCountTime );
            },
            minNumberObservations,
            maxArcGap,
            getEarthGroundStationPositions( bodies ) );
}
}  // namespace observation_models
}  // namespace tudat
