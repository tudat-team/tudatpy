/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <map>

#include <functional>



#include "tudat/simulation/estimation_setup/createObservationModel.h"

namespace tudat
{

namespace observation_models
{


std::vector< LinkDefinition > getObservationModelListLinkEnds(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationModelList )
{
    std::vector< LinkDefinition > linkEndsList;
    for( unsigned int i = 0; i < observationModelList.size( ); i++ )
    {
        linkEndsList.push_back( observationModelList.at( i )->linkEnds_ );
    }
    return linkEndsList;
}

std::map< ObservableType, std::vector< std::shared_ptr< ObservationModelSettings > > > sortObservationModelSettingsByType(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationModelSettings )
{
    std::map< ObservableType, std::vector< std::shared_ptr< ObservationModelSettings > > > sortedObservationModelSettings;
    for( unsigned int i = 0; i < observationModelSettings.size( ); i++ )
    {
        sortedObservationModelSettings[ observationModelSettings.at( i )->observableType_ ].push_back(
                observationModelSettings.at( i ) );
    }
    return sortedObservationModelSettings;
}


std::function< double ( observation_models::FrequencyBands, observation_models::FrequencyBands ) > getTurnaroundFunction(
    const simulation_setup::SystemOfBodies &bodies,
    const LinkEnds& linkEnds )
{
    std::function< double ( observation_models::FrequencyBands, observation_models::FrequencyBands ) > turnaroundRatioFunction;
    // Check if retransmitter is a body
    if ( linkEnds.at( observation_models::retransmitter ).stationName_ == "" || !simulation_setup::isReferencePointGroundStation(
        bodies, linkEnds.at( observation_models::retransmitter ).bodyName_, linkEnds.at( observation_models::retransmitter ).stationName_ ) )
    {
        if ( bodies.getBody( linkEnds.at( observation_models::retransmitter ).bodyName_ )->getVehicleSystems( ) == nullptr )
        {
            throw std::runtime_error(
                "Error when retrieving turnaround ratio: vehicle systems are not "
                "defined for retransmitter link end body " + linkEnds.at( observation_models::retransmitter ).bodyName_ + "." );
        }
        turnaroundRatioFunction = bodies.getBody( linkEnds.at( observation_models::retransmitter ).bodyName_ )->getVehicleSystems(
        )->getTransponderTurnaroundRatio( );
    }
        // If retransmitter is a ground station of the body
    else
    {
        if ( bodies.getBody( linkEnds.at( observation_models::retransmitter ).bodyName_ )->getGroundStation(
            linkEnds.at( observation_models::retransmitter ).stationName_ )->getVehicleSystems( ) == nullptr )
        {
            throw std::runtime_error(
                "Error when retrieving turnaround ratio: vehicle systems are not "
                "defined for retransmitter link end ID " + linkEnds.at( observation_models::retransmitter ).stationName_ + "." );
        }
        turnaroundRatioFunction = bodies.getBody( linkEnds.at( observation_models::retransmitter ).bodyName_ )->getGroundStation(
            linkEnds.at( observation_models::retransmitter ).stationName_ )->getVehicleSystems( )->getTransponderTurnaroundRatio( );
    }
    return turnaroundRatioFunction;
}

std::shared_ptr< ground_stations::StationFrequencyInterpolator > getTransmittingFrequencyInterpolator(
    const simulation_setup::SystemOfBodies &bodies,
    const LinkEnds& linkEnds )
{
    if ( bodies.getBody( linkEnds.at( observation_models::transmitter ).bodyName_ )->getGroundStation(
        linkEnds.at( observation_models::transmitter ).stationName_ )->getTransmittingFrequencyCalculator( ) ==
         nullptr )
    {
        throw std::runtime_error(
            "Error when creating DSN N-way averaged Doppler observation model: transmitted frequency not  "
            "defined for link end station " + linkEnds.at( observation_models::transmitter ).bodyName_ + ", " +
            linkEnds.at( observation_models::transmitter ).stationName_ );
    }
    return bodies.getBody( linkEnds.at( observation_models::transmitter ).bodyName_ )->getGroundStation(
        linkEnds.at( observation_models::transmitter ).stationName_ )->getTransmittingFrequencyCalculator( );
}


} // namespace observation_models

} // namespace tudat

