/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudatpy/kernelExternTemplatesDynamicsSimulators.h"
#include "tudatpy/kernelExternTemplatesEnvironmentUpdater.h"
#include "tudatpy/kernelExternTemplatesObservationManagers.h"
#include "tudatpy/kernelExternTemplatesObservationSimulators.h"
#include "tudatpy/kernelExternTemplatesOrbitDetermination.h"
#include "tudatpy/kernelExternTemplatesParameters.h"
#include "tudatpy/kernelExternTemplatesPropagationOutput.h"
#include "tudatpy/kernelExternTemplatesSimulateObservations.h"
#include "tudatpy/kernelExternTemplatesStateDerivativeModels.h"
#include "tudatpy/kernelExternTemplatesVariationalEquationsSolvers.h"

namespace tudat
{

namespace observation_models
{

template std::map< ObservableType, std::shared_ptr< ObservationManagerBase< STATE_SCALAR_TYPE, TIME_TYPE > > >
createObservationManagersBase< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationSettingsList,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< STATE_SCALAR_TYPE > > fullParameters,
        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const std::shared_ptr< propagators::DependentVariablesInterface< TIME_TYPE > > dependentVariablesInterface );

}  // namespace observation_models

namespace simulation_setup
{

template std::shared_ptr< propagators::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >
createDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >(
        const SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< STATE_SCALAR_TYPE > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated );

template std::shared_ptr< propagators::VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE > >
createVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >(
        const SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< STATE_SCALAR_TYPE > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< STATE_SCALAR_TYPE > > parametersToEstimate,
        const bool integrateEquationsOnCreation );

template std::shared_ptr< EstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE > >
createBestFitToCurrentEphemeris< TIME_TYPE, STATE_SCALAR_TYPE >(
        const SystemOfBodies& bodies,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TIME_TYPE > >& integratorSettings,
        const TIME_TYPE initialTime,
        const TIME_TYPE finalTime,
        const TIME_TYPE dataPointInterval,
        const std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > additionalParameterNames,
        const int numberOfIterations,
        const bool reintegrateVariationalEquations,
        const double resultsPrintFrequency );

template std::shared_ptr< observation_models::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >
simulateObservations< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationsToSimulate,
        const std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >&
                observationSimulators,
        const SystemOfBodies bodies );

}  // namespace simulation_setup

}  // namespace tudat
