/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DIFFERENCEDTIMEOFARRIVALOBSERVATIONMODEL_H
#define TUDAT_DIFFERENCEDTIMEOFARRIVALOBSERVATIONMODEL_H

#include <map>

#include <functional>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType = double, typename TimeType = double >
class OneWayDifferencedTimeOfArrivalObservationModel : public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;

    static std::vector< std::shared_ptr< FullLinkLightTimeCalculator< ObservationScalarType, TimeType > > >
    createFullLinkLightTimeCalculators(
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                    firstReceiverLightTimeCalculator,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                    secondReceiverLightTimeCalculator )
    {
        return std::vector< std::shared_ptr< FullLinkLightTimeCalculator< ObservationScalarType, TimeType > > >{
                std::make_shared< FullLinkLightTimeCalculator< ObservationScalarType, TimeType > >(
                        std::vector< std::shared_ptr<
                                observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > >{
                                firstReceiverLightTimeCalculator },
                        std::make_shared< LightTimeConvergenceCriteria >( ),
                        false ),
                std::make_shared< FullLinkLightTimeCalculator< ObservationScalarType, TimeType > >(
                        std::vector< std::shared_ptr<
                                observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > >{
                                secondReceiverLightTimeCalculator },
                        std::make_shared< LightTimeConvergenceCriteria >( ),
                        false ) };
    }

    OneWayDifferencedTimeOfArrivalObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                    firstReceiverLightTimeCalculator,
            const std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
                    secondReceiverLightTimeCalculator,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr,
            const std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > >& stationStates =
                    std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > >( ),
            const basic_astrodynamics::TimeScales observableTimeScale = basic_astrodynamics::tdb_scale ):
        ObservationModel< 1, ObservationScalarType, TimeType >(
                differenced_time_of_arrival,
                linkEnds,
                observationBiasCalculator,
                createFullLinkLightTimeCalculators( firstReceiverLightTimeCalculator, secondReceiverLightTimeCalculator ) ),
        stationStates_( stationStates ),
        observableTimeScale_( observableTimeScale )
    {
        if( observableTimeScale_ == basic_astrodynamics::utc_scale || observableTimeScale_ == basic_astrodynamics::ut1_scale )
        {
            if( stationStates.count( receiver ) == 0 )
            {
                throw std::runtime_error(
                        "Error when making differenced time of arrival observation model, no state model found for receiver " +
                        linkEnds.at( receiver ).bodyName_ + ", " + linkEnds.at( receiver ).stationName_ );
            }

            if( stationStates.count( receiver2 ) == 0 )
            {
                throw std::runtime_error(
                        "Error when making differenced time of arrival observation model, no state model found for receiver2 " +
                        linkEnds.at( receiver2 ).bodyName_ + ", " + linkEnds.at( receiver2 ).stationName_ );
            }
        }
    }

    //! Destructor
    ~OneWayDifferencedTimeOfArrivalObservationModel( ) {}

    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySetingsInput = nullptr ) override
    {
        linkEndTimes.resize( 3 );
        linkEndStates.resize( 3 );

        TimeType fullPrecisionTimeAtReceiver2;
        if( linkEndAssociatedWithTime == receiver )
        {
            // Calculate reception time at ground station at the start and end of the count interval at reception time.
            std::vector< double > firstLinkEndTimes;
            std::vector< double > secondLinkEndTimes;
            std::vector< Eigen::Matrix< double, 6, 1 > > firstLinkEndStates;
            std::vector< Eigen::Matrix< double, 6, 1 > > secondLinkEndStates;

            std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySetings;
            std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > firstReceiverLightTimeCalculator =
                    getFirstReceiverLightTimeCalculator( );
            std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > secondReceiverLightTimeCalculator =
                    getSecondReceiverLightTimeCalculator( );
            this->setFrequencyProperties( time, receiver, firstReceiverLightTimeCalculator, ancillarySetingsInput, ancillarySetings );
            ObservationScalarType lightTimeForFirstReceiver =
                    this->getFullLinkLightTimeCalculatorFromBase( 0 )->calculateLightTimeWithLinkEndsStates(
                            time, receiver, firstLinkEndTimes, firstLinkEndStates, ancillarySetings );
            TimeType fullPrecisionTransmissionTime = time - static_cast< TimeType >( lightTimeForFirstReceiver );

            linkEndTimes[ 0 ] = firstLinkEndTimes.at( 0 );
            linkEndTimes[ 1 ] = firstLinkEndTimes.at( 1 );

            this->setFrequencyProperties( time, receiver, secondReceiverLightTimeCalculator, ancillarySetingsInput, ancillarySetings );
            ObservationScalarType lightTimeForSecondReceiver =
                    this->getFullLinkLightTimeCalculatorFromBase( 1 )->calculateLightTimeWithLinkEndsStates(
                            fullPrecisionTransmissionTime, transmitter, secondLinkEndTimes, secondLinkEndStates, ancillarySetings );
            fullPrecisionTimeAtReceiver2 =
                    fullPrecisionTransmissionTime + static_cast< TimeType >( lightTimeForSecondReceiver );
            linkEndTimes[ 2 ] = secondLinkEndTimes.at( 1 );

            linkEndStates[ 0 ] = firstLinkEndStates.at( 0 );
            linkEndStates[ 1 ] = firstLinkEndStates.at( 1 );
            linkEndStates[ 2 ] = secondLinkEndStates.at( 1 );
        }
        else
        {
            throw std::runtime_error( "Error in differenced time of arrival, reference link end not recognized" );
        }

        if( observableTimeScale_ == basic_astrodynamics::tdb_scale )
        {
            return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( )
                     << static_cast< ObservationScalarType >( time - fullPrecisionTimeAtReceiver2 ) )
                    .finished( );
        }
        else if( observableTimeScale_ == basic_astrodynamics::utc_scale || observableTimeScale_ == basic_astrodynamics::ut1_scale )
        {
            Eigen::Vector3d nominalReceivingStationState = ( stationStates_.count( receiver ) == 0 )
                    ? Eigen::Vector3d::Zero( )
                    : stationStates_.at( receiver )->getNominalCartesianPosition( );
            Eigen::Vector3d nominalReceivingStationState2 = ( stationStates_.count( receiver2 ) == 0 )
                    ? Eigen::Vector3d::Zero( )
                    : stationStates_.at( receiver2 )->getNominalCartesianPosition( );

            TimeType utTimeAtReceiver = this->timeScaleConverter_->template getCurrentTime< TimeType >(
                    basic_astrodynamics::tdb_scale, observableTimeScale_, time, nominalReceivingStationState );
            TimeType utTimeAtReceiver2 = this->timeScaleConverter_->template getCurrentTime< TimeType >(
                    basic_astrodynamics::tdb_scale, observableTimeScale_, fullPrecisionTimeAtReceiver2, nominalReceivingStationState2 );

            return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( )
                     << static_cast< ObservationScalarType >( utTimeAtReceiver - utTimeAtReceiver2 ) )
                    .finished( );
        }
        else
        {
            TimeType convertedTimeAtReceiver = this->timeScaleConverter_->template getCurrentTime< TimeType >(
                    basic_astrodynamics::tdb_scale, observableTimeScale_, time );
            TimeType convertedTimeAtReceiver2 = this->timeScaleConverter_->template getCurrentTime< TimeType >(
                    basic_astrodynamics::tdb_scale, observableTimeScale_, fullPrecisionTimeAtReceiver2 );

            return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( )
                     << static_cast< ObservationScalarType >( convertedTimeAtReceiver - convertedTimeAtReceiver2 ) )
                    .finished( );
        }
    }

    //! Light time calculator to compute light time at the beginning of the integration time
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getFirstReceiverLightTimeCalculator( )
    {
        return this->getSingleLegLightTimeCalculator( 0, 0 );
    }

    //! Light time calculator to compute light time at the end of the integration time
    std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > > getSecondReceiverLightTimeCalculator( )
    {
        return this->getSingleLegLightTimeCalculator( 1, 0 );
    }

    std::map< std::pair< LinkEndType, LinkEndType >, std::vector< std::shared_ptr< LightTimeCalculatorBase > > > getLegLightTimeCalculators( ) const override
    {
        return { { std::make_pair( transmitter, receiver ), { this->getSingleLegLightTimeCalculator( 0, 0 ) } },
                 { std::make_pair( transmitter, receiver2 ), { this->getSingleLegLightTimeCalculator( 1, 0 ) } } };
    }

private:
    std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > > stationStates_;

    basic_astrodynamics::TimeScales observableTimeScale_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_DIFFERENCEDTIMEOFARRIVALOBSERVATIONMODEL_H
