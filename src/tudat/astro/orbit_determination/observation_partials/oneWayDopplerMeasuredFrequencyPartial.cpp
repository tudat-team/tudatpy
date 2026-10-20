/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/observation_partials/oneWayDopplerMeasuredFrequencyPartial.h"

namespace tudat
{

namespace observation_partials
{

double getOneWayDopplerMeasuredFrequencyScalingFactor(
        const std::shared_ptr< ground_stations::StationFrequencyInterpolator > frequencyInterpolator,
        const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter,
        const std::vector< Eigen::Vector6d >& linkEndStates,
        const std::vector< double >& linkEndTimes,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings )
{
    // Frequency interpolator is queried in UTC to match the observation model
    const double transmitterUtcTime = timeScaleConverter->getCurrentTime< double >(
            basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, linkEndTimes.at( 0 ), linkEndStates.at( 0 ).segment< 3 >( 0 ) );
    return frequencyInterpolator->getTemplatedCurrentFrequency< double, double >( transmitterUtcTime );
}

OneWayDopplerMeasuredFrequencyPartial::OneWayDopplerMeasuredFrequencyPartial(
        const std::shared_ptr< ObservationPartial< 1 > > oneWayDopplerPartial,
        const std::function< double( const std::vector< Eigen::Vector6d >&,
                                     const std::vector< double >&,
                                     const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ) >
                scalingFactorFunction ):
    ObservationPartial< 1 >( oneWayDopplerPartial->getParameterIdentifier( ) ), oneWayDopplerPartial_( oneWayDopplerPartial ),
    scalingFactorFunction_( scalingFactorFunction )
{}

std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > OneWayDopplerMeasuredFrequencyPartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings,
        const Eigen::Vector1d& currentObservation )
{
    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > scaledPartials =
            oneWayDopplerPartial_->calculatePartial( states, times, linkEndOfFixedTime, ancillarySettings, currentObservation );

    double scalingFactor = scalingFactorFunction_( states, times, ancillarySettings );
    for( auto& scaledPartial : scaledPartials )
    {
        scaledPartial.first *= scalingFactor;
    }
    return scaledPartials;
}

}  // namespace observation_partials

}  // namespace tudat
