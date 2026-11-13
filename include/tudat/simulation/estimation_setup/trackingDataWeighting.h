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

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void setObservationWeightsFromTrackingDataScheme(
        const std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > trackingData,
        const simulation_setup::SystemOfBodies& bodies,
        const observation_models::LinkEnds& rawLinkEnds,
        const observation_models::LinkEndType referenceLinkEnd )
{
    if( trackingData->getWeighingScheme( ).empty( ) )
    {
        return;
    }

    // A named scheme is only evaluated when the tracking data contains no explicit weights yet.
    if( trackingData->hasObservationWeightSettings( ) )
    {
        return;
    }

    const std::string weighingScheme = trackingData->getWeighingScheme( );
    if( weighingScheme != "VFCC17" && weighingScheme != "ADES" && weighingScheme != "ADES_NIGHT" )
    {
        throw std::runtime_error( "Error when setting observation weights from TrackingData, weighing scheme '" +
                                  trackingData->getWeighingScheme( ) + "' is not recognized." );
    }

    const Eigen::Vector3d stationPosition = getGroundStationPositionForTrackingDataLinkEnd( bodies, rawLinkEnds.at( referenceLinkEnd ) );
    const std::string stationCode = rawLinkEnds.at( referenceLinkEnd ).getReferencePointName( );
    std::vector< Eigen::Vector3d > stationPositions( trackingData->getNumberOfObservations( ), stationPosition );

    std::vector< std::map< std::string, std::string > > stringMetadata( trackingData->getNumberOfObservations( ) );
    for( unsigned int i = 0; i < trackingData->getNumberOfObservations( ); ++i )
    {
        stringMetadata.at( i )[ "observatory" ] = stationCode;
    }

    const std::map< std::string, std::vector< std::string > > ancillaryStringVectorSettings =
            trackingData->getAncillarySettingsStringVector( );

    for( const std::string metadataKey : { "note2", "catalog", "rmsRA", "rmsDec", "rmsCorr" } )
    {
        auto metadataIterator = ancillaryStringVectorSettings.find( metadataKey );
        if( metadataIterator == ancillaryStringVectorSettings.end( ) )
        {
            continue;
        }
        if( metadataIterator->second.size( ) != trackingData->getNumberOfObservations( ) )
        {
            throw std::runtime_error( "Error when setting " + weighingScheme + " weights from TrackingData, metadata key '" + metadataKey +
                                      "' has an inconsistent number of entries." );
        }
        for( unsigned int i = 0; i < trackingData->getNumberOfObservations( ); ++i )
        {
            stringMetadata.at( i )[ metadataKey ] = metadataIterator->second.at( i );
        }
    }

    if( weighingScheme == "VFCC17" )
    {
        data::setVFCC17Weights< ObservationScalarType, TimeType >( trackingData, stationPositions, stringMetadata );
    }
    else
    {
        data::setAdesWeights< ObservationScalarType, TimeType >(
                trackingData, stationPositions, stringMetadata, weighingScheme == "ADES_NIGHT" );
    }
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

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_TRACKING_DATA_WEIGHTING_H
