/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TRANSMISSIONFREQUENCYINTERFACE_h
#define TUDAT_TRANSMISSIONFREQUENCYINTERFACE_h

#include <map>

#include <functional>

#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/physicalConstants.h"

#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType = double, typename TimeType = Time >
void setTransmissionFrequency(
    const std::shared_ptr< LightTimeCalculator< ObservationScalarType, TimeType > > lightTimeCalculator,
    const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter,
    const std::shared_ptr< ground_stations::StationFrequencyInterpolator > transmittingFrequencyCalculator,
    const TimeType observationTime,
    const LinkEndType linkEndAssociatedWithTime,
    const std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings )
{
    TimeType approximateTransmissionTime, approximateReceptionTime;
    if( linkEndAssociatedWithTime == receiver )
    {
        approximateReceptionTime = observationTime;
        approximateTransmissionTime = observationTime -
          lightTimeCalculator->calculateFirstIterationLightTime( observationTime, true );
    }
    else if( linkEndAssociatedWithTime == transmitter )
    {
        approximateTransmissionTime = observationTime;
        approximateReceptionTime = observationTime +
                                         lightTimeCalculator->calculateFirstIterationLightTime( observationTime, false );
    }
    else
    {
        throw std::runtime_error( "Error when getting transmission frequency for one-way, reference link end is incompatible. " );
    }

    TimeType approximateUtcTransmissionTime = timeScaleConverter->getCurrentTime< TimeType >(
        basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, approximateTdbTransmissionTime );

    double approximateTransmissionFrequency =
        transmittingFrequencyCalculator->getTemplatedCurrentFrequency< double, TimeType >( approximateUtcTransmissionTime );
    ancillarySettings->setIntermediateDoubleData( transmitter_frequency_intermediate, approximateTransmissionFrequency );
    ancillarySettings->setIntermediateDoubleData( received_frequency_intermediate, approximateTransmissionFrequency );

}


template< typename ObservationScalarType = double, typename TimeType = Time >
void setTransmissionReceptionFrequencies(
    const std::shared_ptr< MultiLegLightTimeCalculator< ObservationScalarType, TimeType > > multiLegLightTimeCalculator,
    const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter,
    const std::shared_ptr< ground_stations::StationFrequencyInterpolator > transmittingFrequencyCalculator,
    const TimeType receptionTdbTime,
    const LinkEndType referenceLinkEndType,
    const std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings,
    const double turnAroundratio )
{
    if( referenceLinkEndType != receiver )
    {
        throw std::runtime_error( "Error when getting n-way transmission frequency, reference link end is incompatible. " );
    }
    TimeType approximateTdbTransmissionTime = receptionTdbTime -
                                              multiLegLightTimeCalculator->calculateFirstIterationLightTime( receptionTdbTime, receiver );

    TimeType approximateUtcTransmissionTime = timeScaleConverter->getCurrentTime< TimeType >(
        basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, approximateTdbTransmissionTime );

    double approximateTransmissionFrequency =
        transmittingFrequencyCalculator->getTemplatedCurrentFrequency< double, TimeType >( approximateUtcTransmissionTime );
    ancillarySettings->setIntermediateDoubleData( received_frequency_intermediate, approximateTransmissionFrequency * turnAroundratio );
    ancillarySettings->setIntermediateDoubleData( transmitter_frequency_intermediate, approximateTransmissionFrequency );
}
}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_TRANSMISSIONFREQUENCYINTERFACE_h
