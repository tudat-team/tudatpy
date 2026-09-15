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
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"

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
#include "tudat/simulation/estimation_setup/simulatePseudoObservations.h"

namespace tudat
{

namespace simulation_setup
{

template< typename TimeType = double, typename StateScalarType = double >
std::shared_ptr< EstimationOutput< StateScalarType, TimeType > > createBestFitToCurrentEphemeris(
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
        const int numberOfIterations = 3,
        const bool reintegrateVariationalEquations = true,
        const double resultsPrintFrequency = 0.0 )
{
    using namespace observation_models;
    using namespace estimatable_parameters;
    using namespace propagators;

    TimeType initialPropagationTime = initialTime;

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialState =
            getInitialStatesOfBodies< TimeType, StateScalarType >( bodiesToPropagate, centralBodies, bodies, initialPropagationTime );

    std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >(
                    centralBodies,
                    accelerationModelMap,
                    bodiesToPropagate,
                    initialState,
                    initialPropagationTime,
                    integratorSettings,
                    std::make_shared< PropagationTimeTerminationSettings >( finalTime ) );
    if( resultsPrintFrequency > 0.0 )
    {
        propagatorSettings->getPrintSettings( )->setResultsPrintFrequencyInSteps( resultsPrintFrequency );
    }

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialStateParameterSettings< StateScalarType, TimeType >( propagatorSettings, bodies );
    parameterNames.insert( parameterNames.end( ), additionalParameterNames.begin( ), additionalParameterNames.end( ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType, TimeType >( parameterNames, bodies, propagatorSettings );
    printEstimatableParameterEntries( parametersToEstimate );

    std::pair< std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >,
               std::shared_ptr< observation_models::ObservationDataset< StateScalarType, TimeType > > >
            observationDatasetAndModelSettings = simulatePseudoObservationDataset< TimeType, StateScalarType >(
                    bodies, bodiesToPropagate, centralBodies, initialTime, finalTime, dataPointInterval );
    std::shared_ptr< observation_models::ObservationDataset< StateScalarType, TimeType > > observationDataset =
            observationDatasetAndModelSettings.second;

    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList =
            observationDatasetAndModelSettings.first;

    OrbitDeterminationManager< StateScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< StateScalarType, TimeType >(
                    bodies, parametersToEstimate, observationModelSettingsList, propagatorSettings );

    std::shared_ptr< EstimationInput< StateScalarType, TimeType > > estimationInput =
            std::make_shared< EstimationInput< StateScalarType, TimeType > >( observationDataset );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( numberOfIterations ) );
    estimationInput->defineEstimationSettings( 0, 1, 0, 1, 1, 1 );
    return orbitDeterminationManager.estimateParameters( estimationInput );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_FITORBITTOEPHEMERIS_H
