/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DIFFERENCEDFREQUENCYOFARRIVALOBSERVATIONMODEL_H
#define TUDAT_DIFFERENCEDFREQUENCYOFARRIVALOBSERVATIONMODEL_H

#include <map>
#include <iostream>

#include <functional>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/observation_models/oneWayDopplerMeasuredFrequencyObservationModel.h"
#include "tudat/astro/basic_astro/timeConversions.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType = double, typename TimeType = double >
class OneWayDifferencedFrequencyOfArrivalObservationModel : public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;

    using ObservationModel< 1, ObservationScalarType, TimeType >::timeScaleConverter_;
    using ObservationModel< 1, ObservationScalarType, TimeType >::frequencyInterpolator_;

    OneWayDifferencedFrequencyOfArrivalObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< OneWayDopplerMeasuredFrequencyObservationModel< ObservationScalarType, TimeType > > firstDopplerModel,
            const std::shared_ptr< OneWayDopplerMeasuredFrequencyObservationModel< ObservationScalarType, TimeType > > secondDopplerModel,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr,
            const std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > >& stationStates =
                    std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > >( ),
            const basic_astrodynamics::TimeScales observableTimeScale = basic_astrodynamics::tdb_scale ):
        ObservationModel< 1, ObservationScalarType, TimeType >( differenced_frequency_of_arrival, linkEnds, observationBiasCalculator ),
        firstDopplerModel_( firstDopplerModel ), secondDopplerModel_( secondDopplerModel ), stationStates_( stationStates ),
        observableTimeScale_( observableTimeScale )
    {
        if( linkEnds.size() != 3 )
        {
            throw std::runtime_error(
                    "Error when creating differenced frequency of arrival observation model, exactly 3 link ends required" );
        }
        if( observableTimeScale_ != basic_astrodynamics::tdb_scale )
        {
            throw std::runtime_error(
                    "Error when creating differenced frequency of arrival observation model, only TDB time scale is currently supported" );
        }
        if( firstDopplerModel_ == nullptr || secondDopplerModel_ == nullptr )
        {
            throw std::runtime_error(
                    "Error when creating differenced frequency of arrival observation model, null doppler model pointer" );
        }
        if( firstDopplerModel_->frequencyInterpolator_ == nullptr || secondDopplerModel_->frequencyInterpolator_ == nullptr )
        {
            if( frequencyInterpolator_ == nullptr )
            {
                throw std::runtime_error(
                        "Error when creating differenced frequency of arrival observation model, no frequency interpolator found" );
            }
            firstDopplerModel_->setFrequencyInterpolator( frequencyInterpolator_ );
            secondDopplerModel_->setFrequencyInterpolator( frequencyInterpolator_ );
        }
    }

    //! Destructor
    ~OneWayDifferencedFrequencyOfArrivalObservationModel( ) {}

    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettingsInput = nullptr )
    {
        if( linkEndAssociatedWithTime == transmitter || linkEndAssociatedWithTime == receiver2 )
        {
            throw std::runtime_error(
                    "Error when computing Differenced Frequency of Arrival observables: the selected "
                    "reference link end (" +
                    getLinkEndTypeString( linkEndAssociatedWithTime ) + ") is not valid. Must be the receiver." );
        }

        // Separate link ends for both doppler models
        LinkEnds firstDopplerLinkEnds;
        firstDopplerLinkEnds[ transmitter ] = this->linkEnds_.at( transmitter );
        firstDopplerLinkEnds[ receiver ] = this->linkEnds_.at( receiver );
        LinkEnds secondDopplerLinkEnds;
        secondDopplerLinkEnds[ transmitter ] = this->linkEnds_.at( transmitter );
        secondDopplerLinkEnds[ receiver ] = this->linkEnds_.at( receiver2 );

        // Compute both doppler observables
        std::vector< double > firstLinkEndTimes, secondLinkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > firstLinkEndStates, secondLinkEndStates;
        const auto ancillarySettingsInput = ancilliarySetingsInput;
        Eigen::Matrix< ObservationScalarType, 1, 1 > firstDopplerObservation = firstDopplerModel_->computeIdealObservationsWithLinkEndData(
                time, linkEndAssociatedWithTime, firstLinkEndTimes, firstLinkEndStates, ancillarySettingsInput );
        Eigen::Matrix< ObservationScalarType, 1, 1 > secondDopplerObservation =
                secondDopplerModel_->computeIdealObservationsWithLinkEndData(
                        time, linkEndAssociatedWithTime, secondLinkEndTimes, secondLinkEndStates, ancillarySettingsInput );
        // Combine link end times and states
        // Store all 4 link end times/states: first transmitter, first receiver, second transmitter, second receiver
        linkEndTimes.resize( 4 );
        linkEndTimes[ 0 ] = firstLinkEndTimes[ 0 ];
        linkEndTimes[ 1 ] = firstLinkEndTimes[ 1 ];
        linkEndTimes[ 2 ] = secondLinkEndTimes[ 0 ];
        linkEndTimes[ 3 ] = secondLinkEndTimes[ 1 ];
        linkEndStates.resize( 4 );
        linkEndStates[ 0 ] = firstLinkEndStates[ 0 ];
        linkEndStates[ 1 ] = firstLinkEndStates[ 1 ];
        linkEndStates[ 2 ] = secondLinkEndStates[ 0 ];
        linkEndStates[ 3 ] = secondLinkEndStates[ 1 ];

        // Compute differenced frequency of arrival observation
        Eigen::Matrix< ObservationScalarType, 1, 1 > observation = firstDopplerObservation - secondDopplerObservation;

        return observation;
    }

    [[nodiscard]] auto getFirstDopplerMeasuredFrequencyModel( )
    {
        return firstDopplerModel_;
    }

    [[nodiscard]] auto getSecondDopplerMeasuredFrequencyModel( )
    {
        return secondDopplerModel_;
    }

private:
    const std::shared_ptr< OneWayDopplerMeasuredFrequencyObservationModel< ObservationScalarType, TimeType > > firstDopplerModel_;
    const std::shared_ptr< OneWayDopplerMeasuredFrequencyObservationModel< ObservationScalarType, TimeType > > secondDopplerModel_;

    std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > > stationStates_;

    basic_astrodynamics::TimeScales observableTimeScale_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_DIFFERENCEDFREQUENCYOFARRIVALOBSERVATIONMODEL_H
