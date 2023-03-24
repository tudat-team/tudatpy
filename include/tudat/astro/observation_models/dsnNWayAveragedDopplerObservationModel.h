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

inline double getDsnNWayAveragedDopplerScalingFactor(
        const std::shared_ptr< simulation_setup::SystemOfBodies > bodies,
        const LinkEnds& linkEnds,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::vector< Eigen::Vector6d >& linkEndStates,
        const std::vector< double >& linkEndTimes,
        const std::shared_ptr< ObservationAncilliarySimulationSettings< double > > ancillarySettings,
        const bool isFirstPartial )
{
    double integrationTime;
    double turnaroundRatio;
    try
    {
        integrationTime = ancillarySettings->getAncilliaryDoubleData( doppler_integration_time );
        turnaroundRatio = ancillarySettings->getAncilliaryDoubleData( turnaround_ratio );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error(
                "Error when retrieving integration ancillary settings for DSN N-way averaged Doppler observable: " +
                std::string( caughtException.what( ) ) );
    }

    double transmissionTime;
    if ( isFirstPartial )
    {
        transmissionTime = linkEndTimes.at( 3 );
    }
    else
    {
        transmissionTime = linkEndTimes.at( 0 );
    }

    double frequency = bodies->getBody( linkEnds.at( observation_models::transmitter ).bodyName_ )->getGroundStation(
                linkEnds.at( observation_models::transmitter ).stationName_ )->getTransmittingFrequencyCalculator( )->
                        template getTemplatedCurrentFrequency< double >( transmissionTime );

    return turnaroundRatio * frequency / integrationTime;
}

template< typename ObservationScalarType = double, typename TimeType = Time >
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
        if( !std::is_same< Time, TimeType >::value )
        {
            throw std::runtime_error(
                    "Error when defining DSN N-way averaged Doppler observation model: the selected time type "
                    "is not valid, using it would lead to large numerical errors.");
        }

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
        double turnaroundRatio;
        try
        {
            integrationTime = ancillarySettings->getAncilliaryDoubleData( doppler_integration_time );
            referenceFrequency = ancillarySettings->getAncilliaryDoubleData( doppler_reference_frequency );
            turnaroundRatio = ancillarySettings->getAncilliaryDoubleData( turnaround_ratio );
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
                receptionEndTime, linkEndAssociatedWithTime, arcEndLinkEndTimes, arcEndLinkEndStates,
                ancillarySettings )( 0, 0 ) / physical_constants::getSpeedOfLight< ObservationScalarType >( );

        TimeType transmissionStartTime = receptionStartTime - startLightTime;
        TimeType transmissionEndTime = receptionEndTime - endLightTime;

        double transmitterFrequencyIntegral = bodyWithGroundStations_->getGroundStation(
                this->getLinkEnds( ).at( observation_models::transmitter ).stationName_
                )->getTransmittingFrequencyCalculator( )->template getTemplatedFrequencyIntegral< ObservationScalarType, TimeType >(
                        transmissionStartTime, transmissionEndTime );

        Eigen::Matrix< ObservationScalarType, 1, 1 > observation = ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) <<
                turnaroundRatio * ( referenceFrequency - 1.0 / static_cast< ObservationScalarType >( integrationTime ) *
                transmitterFrequencyIntegral ) ).finished( );

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
