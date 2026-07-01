/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SIMULATEPSEUDOOBSERVATIONS_H
#define TUDAT_SIMULATEPSEUDOOBSERVATIONS_H

#include <memory>
#include <vector>

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tudat
{

namespace simulation_setup
{

template< typename TimeType = double, typename StateScalarType = double >
std::pair< std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >,
           std::shared_ptr< observation_models::ObservationDataset< StateScalarType, TimeType > > >
simulatePseudoObservationDataset( const SystemOfBodies& bodies,
                                  const std::vector< std::string >& bodiesToPropagate,
                                  const std::vector< std::string >& centralBodies,
                                  const std::vector< TimeType >& observationTimes )
{
    using namespace observation_models;

    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;
    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > measurementSimulationInput;

    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ observed_body ] = bodiesToPropagate.at( i );
        linkEnds[ observer ] = centralBodies.at( i );

        observationModelSettingsList.push_back( relativePositionObservableSettings( linkEnds ) );

        measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                relative_position_observable, linkEnds, observationTimes, observed_body ) );
    }

    std::shared_ptr< ObservationSimulatorBase< StateScalarType, TimeType > > observationSimulator =
            createObservationSimulators< StateScalarType, TimeType >( observationModelSettingsList, bodies ).at( 0 );

    std::shared_ptr< observation_models::ObservationDataset< StateScalarType, TimeType > > observationDataset =
            simulateObservationDataset< StateScalarType, TimeType >( measurementSimulationInput, { observationSimulator }, bodies );

    return std::make_pair( observationModelSettingsList, observationDataset );
}

template< typename TimeType = double, typename StateScalarType = double >
std::pair< std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >,
           std::shared_ptr< observation_models::ObservationDataset< StateScalarType, TimeType > > >
simulatePseudoObservationDataset( const SystemOfBodies& bodies,
                                  const std::vector< std::string >& bodiesToPropagate,
                                  const std::vector< std::string >& centralBodies,
                                  const TimeType initialTime,
                                  const TimeType finalTime,
                                  const TimeType dataPointInterval )
{
    std::vector< TimeType > observationTimes;
    TimeType currentTime = initialTime + 3600.0;
    while( currentTime < finalTime - 3600.0 )
    {
        observationTimes.push_back( currentTime );
        currentTime += static_cast< double >( dataPointInterval );
    }

    return simulatePseudoObservationDataset< TimeType, StateScalarType >( bodies, bodiesToPropagate, centralBodies, observationTimes );
}

template< typename TimeType = double, typename StateScalarType = double >
std::pair< std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >,
           std::shared_ptr< observation_models::ObservationCollection< StateScalarType, TimeType > > >
simulatePseudoObservations( const SystemOfBodies& bodies,
                            const std::vector< std::string >& bodiesToPropagate,
                            const std::vector< std::string >& centralBodies,
                            const TimeType initialTime,
                            const TimeType finalTime,
                            const TimeType dataPointInterval )
{
    std::pair< std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >,
               std::shared_ptr< observation_models::ObservationDataset< StateScalarType, TimeType > > >
            datasetResult = simulatePseudoObservationDataset< TimeType, StateScalarType >(
                    bodies, bodiesToPropagate, centralBodies, initialTime, finalTime, dataPointInterval );

    return std::make_pair( datasetResult.first,
                           observation_models::createObservationCollection< StateScalarType, TimeType >( datasetResult.second ) );
}

template< typename TimeType = double, typename StateScalarType = double >
std::pair< std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >,
           std::shared_ptr< observation_models::ObservationCollection< StateScalarType, TimeType > > >
simulatePseudoObservations( const SystemOfBodies& bodies,
                            const std::vector< std::string >& bodiesToPropagate,
                            const std::vector< std::string >& centralBodies,
                            const std::vector< TimeType >& observationTimes )
{
    std::pair< std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >,
               std::shared_ptr< observation_models::ObservationDataset< StateScalarType, TimeType > > >
            datasetResult = simulatePseudoObservationDataset< TimeType, StateScalarType >(
                    bodies, bodiesToPropagate, centralBodies, observationTimes );

    return std::make_pair( datasetResult.first,
                           observation_models::createObservationCollection< StateScalarType, TimeType >( datasetResult.second ) );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_SIMULATEPSEUDOOBSERVATIONS_H
