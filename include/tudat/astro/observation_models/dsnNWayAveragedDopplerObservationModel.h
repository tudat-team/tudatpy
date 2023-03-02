/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DSNNWAYAVERAGEDDOPPLEROBSERVATIONMODEL_H
#define TUDAT_DSNNWAYAVERAGEDDOPPLEROBSERVATIONMODEL_H

#include <stdexcept>
#include <string>

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationFrequencies.h"
#include "tudat/astro/observation_models/nWayRangeObservationModel.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType = double, typename TimeType = double >
class DsnNWayAveragedDopplerObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;

    DsnNWayAveragedDopplerObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcStartObservationModel,
            const std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcEndObservationModel,
            const std::shared_ptr< simulation_setup::Body > bodyWithGroundStations,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr ):
        ObservationModel< 1, ObservationScalarType, TimeType >( dsn_n_way_averaged_doppler , linkEnds, observationBiasCalculator),
        arcStartObservationModel_( arcStartObservationModel ),
        arcEndObservationModel_( arcEndObservationModel ),
        bodyWithGroundStations_( bodyWithGroundStations ),
        numberOfLinkEnds_( linkEnds.size( ) )
    {
        if ( numberOfLinkEnds_ < 1 || numberOfLinkEnds_ > 3 )
        {
            throw std::runtime_error(
                    "Error when defining DSN N-way averaged Doppler observation model: the selected number of link ends (" +
                    std::to_string( numberOfLinkEnds_ ) + ") is not valid. Allowed values: 1, 2, and 3.");
        }
    }

    //! Destructor
    ~DsnNWayAveragedDopplerObservationModel( ) { }

    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncilliarySimulationSettings< TimeType > > ancillarySettings = nullptr )
    {
        std::vector< double > arcStartLinkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > arcStartLinkEndStates;
        std::vector< double > arcEndLinkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > arcEndLinkEndStates;

        TimeType integrationTime;
        ObservationScalarType referenceFrequency;
        FrequencyBands uplinkBand, downlinkBand;
        try
        {
            integrationTime = ancillarySettings->getAncilliaryDoubleData( doppler_integration_time );
            referenceFrequency = ancillarySettings->getAncilliaryDoubleData( doppler_reference_frequency );
            uplinkBand = static_cast< FrequencyBands >( ancillarySettings->getAncilliaryDoubleData( uplink_band ) );
            downlinkBand = static_cast< FrequencyBands >( ancillarySettings->getAncilliaryDoubleData( downlink_band ) );
        }
        catch( std::runtime_error& caughtException )
        {
            throw std::runtime_error(
                    "Error when retrieving integration ancillary settings for DSN N-way averaged Doppler observable: " +
                    std::string( caughtException.what( ) ) );
        }

        TimeType receptionStartTime = time - integrationTime / 2.0;
        TimeType receptionEndTime = time + integrationTime / 2.0;

        TimeType startLightTime = arcStartObservationModel_->computeIdealObservationsWithLinkEndData(
                receptionStartTime, linkEndAssociatedWithTime, arcStartLinkEndTimes, arcStartLinkEndStates,
                ancillarySettings )( 0, 0 ) / physical_constants::getSpeedOfLight< ObservationScalarType >( );
        TimeType endLightTime = arcEndObservationModel_->computeIdealObservationsWithLinkEndData(
                receptionEndTime, linkEndAssociatedWithTime, arcStartLinkEndTimes, arcStartLinkEndStates,
                ancillarySettings )( 0, 0 ) / physical_constants::getSpeedOfLight< ObservationScalarType >( );

        TimeType transmissionStartTime = receptionStartTime - startLightTime;
        TimeType transmissionEndTime = receptionEndTime - endLightTime;

        Eigen::Matrix< ObservationScalarType, 1, 1 > observation =
                ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << referenceFrequency -
                getDsnDefaultTurnaroundRatios( uplinkBand, downlinkBand ) /
                static_cast< ObservationScalarType >( integrationTime ) *
                bodyWithGroundStations_->getGroundStation(
                        this->getLinkEnds( ).at( observation_models::transmitter ).stationName_
                        )->getTransmittingFrequencyCalculator( )->getFrequencyIntegral(
                                transmissionStartTime, transmissionEndTime ) ).finished( );

        std::cerr << time << " " << transmissionStartTime << " " << transmissionEndTime << std::endl;

        linkEndTimes.clear( );
        linkEndStates.clear( );
        linkEndTimes.resize( 4 * ( numberOfLinkEnds_ - 1 ) );
        linkEndStates.resize( 4 * ( numberOfLinkEnds_ - 1 ) );

        for( unsigned int i = 0; i < 2 * ( numberOfLinkEnds_ - 1 ) ; i++ )
        {
            linkEndTimes[ i ] = arcStartLinkEndTimes[ i ];
            linkEndTimes[ i + 2 * ( numberOfLinkEnds_ - 1 ) ] = arcEndLinkEndTimes[ i ];

            linkEndStates[ i ] = arcStartLinkEndStates[ i ];
            linkEndStates[ i + 2 * ( numberOfLinkEnds_ - 1 ) ] = arcEndLinkEndStates[ i ];
        }

        return observation;
    }

    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > getArcEndObservationModel( )
    {
        return arcEndObservationModel_;
    }

    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > getArcStartObservationModel( )
    {
        return arcStartObservationModel_;
    }

private:

    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcStartObservationModel_;

    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcEndObservationModel_;

    std::shared_ptr< simulation_setup::Body > bodyWithGroundStations_;

    unsigned int numberOfLinkEnds_;
};



} // namespace observation_models

} // namespace tudat


#endif //TUDAT_DSNNWAYAVERAGEDDOPPLEROBSERVATIONMODEL_H
