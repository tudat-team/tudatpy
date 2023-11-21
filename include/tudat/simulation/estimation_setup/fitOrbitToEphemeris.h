/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FITORBITTOEPHEMERIS_H
#define TUDAT_FITORBITTOEPHEMERIS_H

#include <vector>

#include <memory>
#include <functional>

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/timeType.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/basics/utilities.h"

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tudat
{

namespace simulation_setup
{

template< typename TimeType = double, typename StateScalarType = double >
std::pair< std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >,
std::shared_ptr< observation_models::ObservationCollection< StateScalarType, TimeType > > > simulatePseudoObservations(
    const SystemOfBodies& bodies,
    const std::vector< std::string >& bodiesToPropagate,
    const std::vector< std::string >& centralBodies,
    const TimeType initialTime,
    const TimeType finalTime,
    const TimeType dataPointInterval  )
{
    using namespace observation_models;

    std::vector< TimeType > observationTimes;
    double currentTime = initialTime;
    while( currentTime < finalTime )
    {
        observationTimes.push_back( currentTime );
        currentTime += dataPointInterval;
    }


    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;
    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > measurementSimulationInput;

    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        // Create link ends
        LinkEnds linkEnds;
        linkEnds[ observed_body ] = bodiesToPropagate.at( i );
        linkEnds[ observer ] = centralBodies.at( i );

        // Create observation model settings
        observationModelSettingsList.push_back( relativePositionObservableSettings( linkEnds ) );

        measurementSimulationInput.push_back(
            std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                relative_position_observable, linkEnds, observationTimes, observed_body ) );
    }

    std::shared_ptr< ObservationSimulatorBase<StateScalarType, TimeType> > observationSimulator =
        createObservationSimulators< StateScalarType, TimeType >( observationModelSettingsList, bodies ).at( 0 );

    std::shared_ptr< observation_models::ObservationCollection< StateScalarType, TimeType > > observationCollection =
        simulateObservations< StateScalarType, TimeType >( measurementSimulationInput, { observationSimulator }, bodies );

    return std::make_pair( observationModelSettingsList, observationCollection );
}

template< typename TimeType = double, typename StateScalarType = double >
std::shared_ptr< EstimationOutput< > > createBestFitToCurrentEphemeris(
    const SystemOfBodies& bodies,
    const basic_astrodynamics::AccelerationMap& accelerationModelMap,
    const std::vector< std::string >& bodiesToPropagate,
    const std::vector< std::string >& centralBodies,
    const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > >& integratorSettings,
    const TimeType initialTime,
    const TimeType finalTime,
    const TimeType dataPointInterval,
    const std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > additionalParameterNames =
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >( ),
    const int numberOfIterations = 3 )
{
    using namespace observation_models;
    using namespace estimatable_parameters;
    using namespace propagators;

    double initialPropagationTime = initialTime;

   Eigen::VectorXd initialState = getInitialStatesOfBodies( bodiesToPropagate, centralBodies, bodies, initialPropagationTime );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
        std::make_shared< TranslationalStatePropagatorSettings< double > >
        ( centralBodies, accelerationModelMap, bodiesToPropagate, initialState, initialPropagationTime, integratorSettings,
            std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
        getInitialStateParameterSettings< double >( propagatorSettings, bodies );
    parameterNames.insert(parameterNames.end(), additionalParameterNames.begin(), additionalParameterNames.end());

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
        createParametersToEstimate< StateScalarType, TimeType >( parameterNames, bodies, propagatorSettings );


    std::pair< std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >,
        std::shared_ptr< observation_models::ObservationCollection< StateScalarType, TimeType > > >
        observationCollectionAndModelSettings = simulatePseudoObservations(
            bodies, bodiesToPropagate, centralBodies, initialTime, finalTime, dataPointInterval  );
    std::shared_ptr< observation_models::ObservationCollection< > > observationCollection = observationCollectionAndModelSettings.second;

    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList =
        observationCollectionAndModelSettings.first;

    OrbitDeterminationManager< StateScalarType, TimeType > orbitDeterminationManager = OrbitDeterminationManager< StateScalarType, TimeType >(
        bodies, parametersToEstimate, observationModelSettingsList, propagatorSettings );

    std::shared_ptr< EstimationInput< > > estimationInput = std::make_shared< EstimationInput< > >(
        observationCollection );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( numberOfIterations ) );
    estimationInput->defineEstimationSettings(
        0, 1, 0, 1, 1, 1 );
    return  orbitDeterminationManager.estimateParameters( estimationInput );

}

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_FITORBITTOEPHEMERIS_H
