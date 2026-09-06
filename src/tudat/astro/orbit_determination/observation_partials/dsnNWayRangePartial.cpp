/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          T. Moyer (2003), Formulation for Observed and Computed Values of Deep Space Network Data
 * Types for Navigation, DEEP SPACE COMMUNICATIONS AND NAVIGATION SERIES, JPL/NASA
 */

#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/observation_models/dsnNWayRangeObservationModel.h"
#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/astro/observation_models/observationFrequencies.h"
#include "tudat/astro/orbit_determination/observation_partials/dsnNWayRangePartial.h"

namespace tudat
{

namespace observation_partials
{

double getDsnNWayRangeScalingFactor(
        const std::function< double( std::vector< observation_models::FrequencyBands >, double ) > transmittedFrequencyFunction,
        const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter,
        const std::vector< double >& linkEndTimes,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings )
{
    if( ancillarySettings == nullptr )
    {
        throw std::runtime_error(
                "Error when computing DSN N-way range partials scaling factor; no ancillary "
                "settings found. " );
    }

    std::vector< observation_models::FrequencyBands > frequencyBands;
    try
    {
        frequencyBands = observation_models::convertDoubleVectorToFrequencyBands(
                ancillarySettings->getAncillaryDoubleVectorData( observation_models::frequency_bands ) );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error when retrieving ancillary settings for DSN N-way range partial: " +
                                  std::string( caughtException.what( ) ) );
    }

    if( frequencyBands.empty( ) )
    {
        throw std::runtime_error(
                "Error when computing DSN N-way range partials scaling factor; no frequency bands "
                "found in ancillary settings." );
    }
    double conversionFactor = observation_models::getDsnRangeUnitConversionFactor( frequencyBands.at( 0 ) );

    double utcTransmissionTime = timeScaleConverter->getCurrentTime< double >(
            basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, linkEndTimes.at( 0 ), Eigen::Vector3d::Zero( ) );
    double transmittedFrequency = transmittedFrequencyFunction( frequencyBands, utcTransmissionTime );

    return conversionFactor * transmittedFrequency / physical_constants::getSpeedOfLight< double >( );
}

DsnNWayRangePartial::DsnNWayRangePartial(
        const std::shared_ptr< NWayRangePartial > nWayRangePartial,
        const std::function< double( const std::vector< double >&,
                                     const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ) >
                scalingFactorFunction ):
    ObservationPartial< 1 >( nWayRangePartial->getParameterIdentifier( ) ), nWayRangePartial_( nWayRangePartial ),
    scalingFactorFunction_( scalingFactorFunction )
{}

//! Function to calculate the observation partial(s) at required time and state
/*!
 *  Function to calculate the observation partial(s) at required time and state. State and time
 *  are typically obtained from evaluation of observation model.
 *  \param states Link end states. Index maps to link end for a given ObsevableType through getLinkEndIndex function.
 *  \param times Link end time.
 *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
 *  \param ancillarySettings Observation ancillary simulation settings.
 *  \param currentObservation Value of the observation for which the partial is to be computed.
 *  \return Vector of pairs containing partial values and associated times.
 */
std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > DsnNWayRangePartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings,
        const Eigen::Vector1d& currentObservation )
{
    // The observation model (and hence the scaling factor computation) is only valid for observations
    // referenced to the reception time
    if( linkEndOfFixedTime != observation_models::receiver )
    {
        throw std::runtime_error( "Error when computing DSN N-way range partials: the selected reference link end (" +
                                  observation_models::getLinkEndTypeString( linkEndOfFixedTime ) + ") is not valid." );
    }

    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > scaledPartials =
            nWayRangePartial_->calculatePartial( states, times, linkEndOfFixedTime, ancillarySettings, currentObservation );

    double scalingFactor = scalingFactorFunction_( times, ancillarySettings );
    for( auto& scaledPartial : scaledPartials )
    {
        scaledPartial.first *= scalingFactor;
    }
    return scaledPartials;
}

}  // namespace observation_partials

}  // namespace tudat
