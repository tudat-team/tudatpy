/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TRACKING_DATA_WEIGHTING_H
#define TUDAT_TRACKING_DATA_WEIGHTING_H

#include <Eigen/Core>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/trackingDataWeights.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace observation_models
{

Eigen::Vector3d getGroundStationPositionForTrackingDataLinkEnd( const simulation_setup::SystemOfBodies& bodies,
                                                                const LinkEndId& linkEndId );

//! Function to retrieve, for each observation, the metadata used by the VFCC17 weights (observatory, note2 and catalog)
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::vector< std::map< std::string, std::string > > getVfcc17Metadata(
        const std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > trackingData,
        const std::string& stationCode )
{
    std::vector< std::map< std::string, std::string > > stringMetadata( trackingData->getNumberOfObservations( ) );
    for( unsigned int i = 0; i < trackingData->getNumberOfObservations( ); ++i )
    {
        stringMetadata.at( i )[ "observatory" ] = stationCode;
    }

    const std::map< std::string, std::vector< std::string > >& ancillarySettingsStringVector =
            trackingData->getAncillarySettingsStringVector( );
    for( const std::string metadataKey : { "note2", "catalog" } )
    {
        auto metadataIterator = ancillarySettingsStringVector.find( metadataKey );
        if( metadataIterator == ancillarySettingsStringVector.end( ) )
        {
            continue;
        }
        if( metadataIterator->second.size( ) != trackingData->getNumberOfObservations( ) )
        {
            throw std::runtime_error( "Error when setting VFCC17 weights from TrackingData, metadata key '" + metadataKey +
                                      "' has an inconsistent number of entries." );
        }
        for( unsigned int i = 0; i < trackingData->getNumberOfObservations( ); ++i )
        {
            stringMetadata.at( i )[ metadataKey ] = metadataIterator->second.at( i );
        }
    }
    return stringMetadata;
}

//! Function to compute the weights of the observations in a TrackingData object from its weight settings, when creating an
//! observation dataset
/*!
 * Function to compute the weights of the observations in a TrackingData object from its weight settings: the explicit weights, or
 * (if these are not defined) the named weights, and the reduction of the weights per block of observations.
 * \param trackingData Tracking data with the weight settings
 * \param bodies System of bodies (used for the station position)
 * \param rawLinkEnds Link ends of the observations
 * \param referenceLinkEnd Reference link end of the observations
 * \return Weights of the observations (nullptr if no weights are defined, in which case unit weights are used)
 */
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< ObservationWeightSettings > computeObservationWeights(
        const std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > trackingData,
        const simulation_setup::SystemOfBodies& bodies,
        const observation_models::LinkEnds& rawLinkEnds,
        const observation_models::LinkEndType referenceLinkEnd )
{
    const std::shared_ptr< data::TrackingDataWeightSettings > weightSettings = trackingData->getWeightSettings( );
    const std::shared_ptr< ObservationWeightSettings > explicitWeights = weightSettings->getObservationWeightSettings( );
    if( explicitWeights == nullptr && weightSettings->getNamedWeightingMethod( ) == data::no_named_weighting )
    {
        if( weightSettings->getBlockDeweighting( ) != nullptr )
        {
            throw std::runtime_error(
                    "Error when computing TrackingData weights, a reduction of the weights per block is defined, but no explicit or "
                    "named weights." );
        }
        return nullptr;
    }

    // Explicit weights, without further contributions, are used directly
    if( explicitWeights != nullptr && weightSettings->getBlockDeweighting( ) == nullptr )
    {
        return explicitWeights;
    }

    const unsigned int numberOfObservations = trackingData->getNumberOfObservations( );
    const unsigned int singleObservationSize = trackingData->getSingleObservationSize( );
    const std::vector< TimeType >& epochs = trackingData->getObservationEpochs( );

    // Number of observations in the same local night (used by the VFCC17 weights and the reduction of the weights per local night)
    const std::vector< int > observationsPerLocalNight = data::getNumberOfObservationsPerLocalNight(
            epochs,
            std::vector< Eigen::Vector3d >(
                    numberOfObservations, getGroundStationPositionForTrackingDataLinkEnd( bodies, rawLinkEnds.at( referenceLinkEnd ) ) ) );
    if( weightSettings->getBlockDeweighting( ) != nullptr &&
        ( weightSettings->getBlockDeweighting( )->blockDefinition_ != data::local_night ||
          weightSettings->getBlockDeweighting( )->rule_ != data::divide_by_number_of_observations ) )
    {
        throw std::runtime_error( "Error when computing TrackingData weights, block definition or deweighting rule not recognized." );
    }
    const double deweightingExponent = weightSettings->getBlockDeweighting( ) != nullptr ? 1.0 : 0.0;

    // Set-level explicit weights: covariance of each observation multiplied by N (W -> D W D, with D containing 1/sqrt(N) for each
    // observation)
    if( std::shared_ptr< SetBlockWeightSettings > setLevelWeights = std::dynamic_pointer_cast< SetBlockWeightSettings >( explicitWeights ) )
    {
        Eigen::VectorXd scaling( numberOfObservations * singleObservationSize );
        for( unsigned int i = 0; i < numberOfObservations; ++i )
        {
            scaling.segment( i * singleObservationSize, singleObservationSize )
                    .setConstant( 1.0 / std::sqrt( static_cast< double >( observationsPerLocalNight.at( i ) ) ) );
        }
        return setBlockWeights( scaling.asDiagonal( ) * setLevelWeights->getWeightMatrix( ) * scaling.asDiagonal( ) );
    }

    // Weights per observation: explicit weights, or named weights
    std::vector< std::map< std::string, std::string > > vfcc17Metadata;
    if( explicitWeights == nullptr )
    {
        if( weightSettings->getNamedWeightingMethod( ) != data::vfcc17 )
        {
            throw std::runtime_error( "Error when computing TrackingData weights, named weighting method not recognized." );
        }
        vfcc17Metadata = getVfcc17Metadata( trackingData, rawLinkEnds.at( referenceLinkEnd ).getReferencePointName( ) );
    }
    std::vector< Eigen::MatrixXd > weightBlocks( numberOfObservations );
    for( unsigned int i = 0; i < numberOfObservations; ++i )
    {
        weightBlocks.at( i ) = explicitWeights != nullptr
                ? explicitWeights->getWeightBlock( i, singleObservationSize )
                : Eigen::MatrixXd( data::getVfcc17Weight( epochs.at( i ), vfcc17Metadata.at( i ), observationsPerLocalNight.at( i ) ) *
                                   Eigen::MatrixXd::Identity( singleObservationSize, singleObservationSize ) );
        weightBlocks.at( i ) /= std::pow( static_cast< double >( observationsPerLocalNight.at( i ) ), deweightingExponent );
    }

    // Keep diagonal weights (named weights, or explicit weights other than weight blocks) diagonal
    if( std::dynamic_pointer_cast< ConstantBlockWeightSettings >( explicitWeights ) == nullptr &&
        std::dynamic_pointer_cast< BlockPerObservationWeightSettings >( explicitWeights ) == nullptr )
    {
        std::vector< Eigen::VectorXd > diagonalWeights;
        for( const Eigen::MatrixXd& weightBlock : weightBlocks )
        {
            diagonalWeights.push_back( weightBlock.diagonal( ) );
        }
        return diagonalPerObservationWeights( diagonalWeights );
    }
    return blockPerObservationWeights( weightBlocks );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_TRACKING_DATA_WEIGHTING_H
