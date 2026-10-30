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

    for( const std::string metadataKey : { "note2", "catalog", "rmsRA", "rmsDec" } )
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
        data::setAdesWeights< ObservationScalarType, TimeType >( trackingData, stationPositions, stringMetadata, weighingScheme == "ADES_NIGHT" );
    }
}

//! Function to compute the rates of right ascension and declination of an observed body
/*!
 * Function to compute the rates of right ascension and declination of an observed body, from its Cartesian state relative to the
 * observer (expressed in the frame in which the right ascension and declination are defined).
 * \param relativeState Cartesian state of the observed body w.r.t. the observer
 * \return Rates of right ascension and declination (in rad/s)
 */
inline Eigen::Vector2d computeRightAscensionAndDeclinationRates( const Eigen::Vector6d& relativeState )
{
    const double x = relativeState( 0 ), y = relativeState( 1 ), z = relativeState( 2 );
    const double vx = relativeState( 3 ), vy = relativeState( 4 ), vz = relativeState( 5 );
    const double squaredEquatorialDistance = x * x + y * y;
    const double squaredDistance = squaredEquatorialDistance + z * z;

    return Eigen::Vector2d(
            ( x * vy - y * vx ) / squaredEquatorialDistance,
            ( vz * squaredEquatorialDistance - z * ( x * vx + y * vy ) ) / ( squaredDistance * std::sqrt( squaredEquatorialDistance ) ) );
}

//! Function to compute the rates of right ascension and declination of the observed body (transmitter), as seen by the observer
//! (receiver), at a list of times, from the geometric states of both link ends in the environment
template< typename TimeType = double >
std::vector< Eigen::Vector2d > getRightAscensionAndDeclinationRates( const simulation_setup::SystemOfBodies& bodies,
                                                                     const observation_models::LinkEnds& rawLinkEnds,
                                                                     const std::vector< TimeType >& observationTimes )
{
    const std::function< Eigen::Vector6d( const TimeType ) > observerStateFunction =
            simulation_setup::getLinkEndCompleteEphemerisFunction< TimeType, double >( rawLinkEnds.at( receiver ), bodies );
    const std::function< Eigen::Vector6d( const TimeType ) > observedBodyStateFunction =
            simulation_setup::getLinkEndCompleteEphemerisFunction< TimeType, double >( rawLinkEnds.at( transmitter ), bodies );

    std::vector< Eigen::Vector2d > angularRates;
    for( const TimeType& observationTime : observationTimes )
    {
        angularRates.push_back( computeRightAscensionAndDeclinationRates( observedBodyStateFunction( observationTime ) -
                                                                          observerStateFunction( observationTime ) ) );
    }
    return angularRates;
}

//! Function to add the contribution of observation epoch uncertainties to the weights of angular position observations
/*!
 * Function to add the contribution of observation epoch uncertainties to the weights of angular position observations, following
 * Farnocchia et al. (2022). For each observation, the covariance C_P (inverse of the positional weights) is augmented to
 * C = C_P + sigma_T^2 r r^T, where sigma_T is the epoch uncertainty and r the vector of right ascension and declination rates.
 * The returned weights are the inverse of C, for each observation.
 * \param positionalWeightSettings Weights of the angular position observations, without contribution of the epoch uncertainties
 * \param timeUncertainties Uncertainties of the observation epochs (in s)
 * \param angularRates Rates of right ascension and declination (in rad/s), at each observation epoch
 * \return Weights of the angular position observations, including the contribution of the epoch uncertainties
 */
inline ObservationWeightSettings getAngularPositionWeightsWithTimeUncertainty( const ObservationWeightSettings& positionalWeightSettings,
                                                                                const std::vector< double >& timeUncertainties,
                                                                                const std::vector< Eigen::Vector2d >& angularRates )
{
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
