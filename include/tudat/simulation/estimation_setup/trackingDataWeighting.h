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
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/orbit_determination/observation_partials/angularPositionPartial.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/trackingDataWeights.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"

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

//! Function to compute the rates of right ascension and declination (w.r.t. J2000) of the observed body (transmitter), as seen by
//! the observer (receiver), at a list of times, from the geometric states of both link ends in the environment
template< typename TimeType = double >
std::vector< Eigen::Vector2d > getRightAscensionAndDeclinationRates( const simulation_setup::SystemOfBodies& bodies,
                                                                     const observation_models::LinkEnds& rawLinkEnds,
                                                                     const std::vector< TimeType >& observationTimes )
{
    const std::function< Eigen::Vector6d( const TimeType ) > observerStateFunction =
            simulation_setup::getLinkEndCompleteEphemerisFunction< TimeType, double >( rawLinkEnds.at( receiver ), bodies );
    const std::function< Eigen::Vector6d( const TimeType ) > observedBodyStateFunction =
            simulation_setup::getLinkEndCompleteEphemerisFunction< TimeType, double >( rawLinkEnds.at( transmitter ), bodies );

    // Rotation from the global frame to J2000, in which the right ascension and declination are defined
    Eigen::Matrix3d rotationToJ2000;
    try
    {
        rotationToJ2000 = reference_frames::getJ2000ToFrameOrientationTransformationMatrix( bodies.getFrameOrientation( ) ).transpose( );
    }
    catch( const std::runtime_error& caughtException )
    {
        throw std::runtime_error(
                "Error when computing right ascension and declination rates for observation time biases or "
                "uncertainties: " +
                std::string( caughtException.what( ) ) );
    }

    std::vector< Eigen::Vector2d > angularRates;
    angularRates.reserve( observationTimes.size( ) );
    for( const TimeType& observationTime : observationTimes )
    {
        // The rates are the partials w.r.t. the observer position, multiplied by the velocity of the observer w.r.t. the observed body
        const Eigen::Vector6d relativeState = observerStateFunction( observationTime ) - observedBodyStateFunction( observationTime );
        angularRates.push_back( observation_partials::calculatePartialOfAngularPositionWrtLinkEndPosition(
                                        rotationToJ2000 * relativeState.segment< 3 >( 0 ), true ) *
                                rotationToJ2000 * relativeState.segment< 3 >( 3 ) );
    }
    return angularRates;
}

//! Function to add the contribution of observation epoch uncertainties to the weights of angular position observations
/*!
 * Function to add the contribution of observation epoch uncertainties to the weights of angular position observations, following
 * Farnocchia et al. (2022). For each observation, the covariance C_P (inverse of the positional weights) is augmented to
 * C = C_P + sigma_T^2 r r^T, where sigma_T is the epoch uncertainty and r the vector of right ascension and declination rates.
 * The returned weights are the inverse of C, for each observation. For set-level weights (correlating different observations), the
 * covariance of the full set is augmented with these 2x2 blocks along its diagonal, and its inverse is returned as set-level weights.
 * \param positionalWeightSettings Weights of the angular position observations, without contribution of the epoch uncertainties
 * \param timeUncertainties Uncertainties of the observation epochs (in s)
 * \param angularRates Rates of right ascension and declination (in rad/s), at each observation epoch
 * \return Weights of the angular position observations, including the contribution of the epoch uncertainties
 */
inline ObservationWeightSettings getAngularPositionWeightsWithTimeUncertainty( const ObservationWeightSettings& positionalWeightSettings,
                                                                               const std::vector< double >& timeUncertainties,
                                                                               const std::vector< Eigen::Vector2d >& angularRates )
{
    if( positionalWeightSettings.type_ == ObservationWeightSettings::WeightsBlockType::set_block )
    {
        Eigen::MatrixXd covariance = positionalWeightSettings.weightBlock_.inverse( );
        for( unsigned int i = 0; i < angularRates.size( ); ++i )
        {
            covariance.block( 2 * i, 2 * i, 2, 2 ) +=
                    timeUncertainties.at( i ) * timeUncertainties.at( i ) * angularRates.at( i ) * angularRates.at( i ).transpose( );
        }
        return ObservationWeightSettings::setBlock( covariance.inverse( ) );
    }

    std::vector< Eigen::MatrixXd > weightBlocks;
    for( unsigned int i = 0; i < angularRates.size( ); ++i )
    {
        const Eigen::Matrix2d covariance = positionalWeightSettings.getObservationWeightBlock( i, 2 ).inverse( ) +
                timeUncertainties.at( i ) * timeUncertainties.at( i ) * angularRates.at( i ) * angularRates.at( i ).transpose( );
        weightBlocks.push_back( covariance.inverse( ) );
    }
    return ObservationWeightSettings::blockPerObservation( weightBlocks );
}

//! Function to compute the weights of the observations in a TrackingData object from its weight settings, when creating an
//! observation dataset
/*!
 * Function to compute the weights of the observations in a TrackingData object from its weight settings: the positional weights,
 * the named weights (for the observations without positional weights), the reduction of the weights per block of observations,
 * and the contribution of the uncertainties of the observation epochs (in this order).
 * \param trackingData Tracking data with the weight settings
 * \param bodies System of bodies (used for the station position and the rates of the observations)
 * \param rawLinkEnds Link ends of the observations
 * \param referenceLinkEnd Reference link end of the observations
 * \param angularRates Rates of right ascension and declination at each observation epoch (only used for astrometric time
 * uncertainties)
 * \return Weights of the observations (empty if no weights are defined, in which case unit weights are used)
 */
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::optional< ObservationWeightSettings > computeObservationWeights(
        const std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > trackingData,
        const simulation_setup::SystemOfBodies& bodies,
        const observation_models::LinkEnds& rawLinkEnds,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::vector< Eigen::Vector2d >& angularRates )
{
    const std::shared_ptr< data::TrackingDataWeightSettings > weightSettings = trackingData->getWeightSettings( );
    const std::shared_ptr< data::PositionalWeightSettings > positionalWeights = weightSettings->positionalWeights_;
    if( positionalWeights == nullptr && weightSettings->namedWeights_ == nullptr )
    {
        if( weightSettings->blockDeweighting_ != nullptr || weightSettings->timingUncertainties_ != nullptr )
        {
            throw std::runtime_error(
                    "Error when computing TrackingData weights, a reduction of the weights per block or timing uncertainties are "
                    "defined, but no positional or named weights." );
        }
        return std::nullopt;
    }

    // Positional weights for all observations, without further contributions, are used directly
    if( positionalWeights != nullptr && positionalWeights->observationIndices_.empty( ) && weightSettings->blockDeweighting_ == nullptr &&
        weightSettings->timingUncertainties_ == nullptr )
    {
        return positionalWeights->weights_;
    }

    const unsigned int numberOfObservations = trackingData->getNumberOfObservations( );
    const unsigned int singleObservationSize = trackingData->getSingleObservationSize( );
    const std::vector< TimeType >& epochs = trackingData->getObservationEpochs( );

    // Number of observations in the same local night (used by the VFCC17 weights and the reduction of the weights per local night)
    std::vector< int > observationsPerLocalNight;
    if( weightSettings->namedWeights_ != nullptr || weightSettings->blockDeweighting_ != nullptr )
    {
        observationsPerLocalNight = data::getNumberOfObservationsPerLocalNight(
                epochs,
                std::vector< Eigen::Vector3d >(
                        numberOfObservations,
                        getGroundStationPositionForTrackingDataLinkEnd( bodies, rawLinkEnds.at( referenceLinkEnd ) ) ) );
    }
    if( weightSettings->blockDeweighting_ != nullptr &&
        ( weightSettings->blockDeweighting_->blockDefinition_ != data::ObservationBlockDefinition::local_night ||
          weightSettings->blockDeweighting_->rule_ != data::BlockDeweightingRule::divide_by_number_of_observations ) )
    {
        throw std::runtime_error( "Error when computing TrackingData weights, block definition or deweighting rule not recognized." );
    }

    ObservationWeightSettings weights;
    if( positionalWeights != nullptr && positionalWeights->weights_.type_ == ObservationWeightSettings::WeightsBlockType::set_block )
    {
        Eigen::MatrixXd weightMatrix = positionalWeights->weights_.weightBlock_;
        if( weightSettings->blockDeweighting_ != nullptr )
        {
            // Covariance of each observation multiplied by N (W -> D W D, with D containing 1/sqrt(N) for each observation)
            Eigen::VectorXd scaling( numberOfObservations * singleObservationSize );
            for( unsigned int i = 0; i < numberOfObservations; ++i )
            {
                scaling.segment( i * singleObservationSize, singleObservationSize )
                        .setConstant( 1.0 / std::sqrt( static_cast< double >( observationsPerLocalNight.at( i ) ) ) );
            }
            weightMatrix = scaling.asDiagonal( ) * weightMatrix * scaling.asDiagonal( );
        }
        weights = ObservationWeightSettings::setBlock( weightMatrix );
    }
    else
    {
        std::vector< Eigen::MatrixXd > weightBlocks( numberOfObservations );
        std::vector< bool > hasPositionalWeights( numberOfObservations, false );
        if( positionalWeights != nullptr )
        {
            const std::vector< unsigned int >& observationIndices = positionalWeights->observationIndices_;
            const unsigned int numberOfWeights = observationIndices.empty( ) ? numberOfObservations : observationIndices.size( );
            for( unsigned int i = 0; i < numberOfWeights; ++i )
            {
                const unsigned int observationIndex = observationIndices.empty( ) ? i : observationIndices.at( i );
                weightBlocks.at( observationIndex ) = positionalWeights->weights_.getObservationWeightBlock( i, singleObservationSize );
                hasPositionalWeights.at( observationIndex ) = true;
            }
        }

        // Named weights for the observations without positional weights
        std::vector< std::map< std::string, std::string > > vfcc17Metadata;
        if( weightSettings->namedWeights_ != nullptr )
        {
            if( weightSettings->namedWeights_->method_ != data::NamedWeightingMethod::vfcc17 )
            {
                throw std::runtime_error( "Error when computing TrackingData weights, named weighting method not recognized." );
            }
            vfcc17Metadata = getVfcc17Metadata( trackingData, rawLinkEnds.at( referenceLinkEnd ).getReferencePointName( ) );
        }
        for( unsigned int i = 0; i < numberOfObservations; ++i )
        {
            if( !hasPositionalWeights.at( i ) )
            {
                if( weightSettings->namedWeights_ == nullptr )
                {
                    throw std::runtime_error( "Error when computing TrackingData weights, observation " + std::to_string( i ) +
                                              " has no positional weights, and no named weights are defined." );
                }
                weightBlocks.at( i ) = data::getVfcc17Weight( epochs.at( i ), vfcc17Metadata.at( i ), observationsPerLocalNight.at( i ) ) *
                        Eigen::MatrixXd::Identity( singleObservationSize, singleObservationSize );
            }
            if( weightSettings->blockDeweighting_ != nullptr )
            {
                weightBlocks.at( i ) /= static_cast< double >( observationsPerLocalNight.at( i ) );
            }
        }
        weights = ObservationWeightSettings::blockPerObservation( weightBlocks );
    }

    // Contribution of the uncertainties of the observation epochs
    if( weightSettings->timingUncertainties_ != nullptr )
    {
        if( std::dynamic_pointer_cast< data::AstrometricTimeUncertaintySettings >( weightSettings->timingUncertainties_ ) == nullptr )
        {
            throw std::runtime_error( "Error when computing TrackingData weights, type of timing uncertainties not recognized." );
        }
        weights = getAngularPositionWeightsWithTimeUncertainty(
                weights, weightSettings->timingUncertainties_->timeUncertainties_, angularRates );
    }
    return weights;
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_TRACKING_DATA_WEIGHTING_H
