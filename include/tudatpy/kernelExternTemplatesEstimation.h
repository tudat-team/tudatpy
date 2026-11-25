/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_ESTIMATION_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_ESTIMATION_H

// Explicit instantiation declarations of the heavy estimation templates for the kernel flavour
// < STATE_SCALAR_TYPE, TIME_TYPE >. Each is instantiated once, in src/tudatpy/kernelInstantiation*.cpp.

#include "tudatpy/kernelExternTemplatesParameters.h"
#include "tudatpy/kernelExternTemplatesPropagation.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/estimation_setup/createObservationManager.h"
#include "tudat/simulation/estimation_setup/createStateDerivativePartials.h"
#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tudat
{

namespace propagators
{

extern template class VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class MultiArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class HybridArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >;

}  // namespace propagators

namespace observation_models
{

extern template std::shared_ptr< ObservationManagerBase< STATE_SCALAR_TYPE, TIME_TYPE > >
createObservationManager< 1, STATE_SCALAR_TYPE, TIME_TYPE >(
        const ObservableType observableType,
        const std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettingsList,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< STATE_SCALAR_TYPE > > parametersToEstimate,
        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const std::shared_ptr< propagators::DependentVariablesInterface< TIME_TYPE > > dependentVariablesInterface );

extern template std::shared_ptr< ObservationManagerBase< STATE_SCALAR_TYPE, TIME_TYPE > >
createObservationManager< 2, STATE_SCALAR_TYPE, TIME_TYPE >(
        const ObservableType observableType,
        const std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettingsList,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< STATE_SCALAR_TYPE > > parametersToEstimate,
        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const std::shared_ptr< propagators::DependentVariablesInterface< TIME_TYPE > > dependentVariablesInterface );

extern template std::shared_ptr< ObservationManagerBase< STATE_SCALAR_TYPE, TIME_TYPE > >
createObservationManager< 3, STATE_SCALAR_TYPE, TIME_TYPE >(
        const ObservableType observableType,
        const std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettingsList,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< STATE_SCALAR_TYPE > > parametersToEstimate,
        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const std::shared_ptr< propagators::DependentVariablesInterface< TIME_TYPE > > dependentVariablesInterface );

extern template std::map< ObservableType, std::shared_ptr< ObservationManagerBase< STATE_SCALAR_TYPE, TIME_TYPE > > >
createObservationManagersBase< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationSettingsList,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< STATE_SCALAR_TYPE > > fullParameters,
        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const std::shared_ptr< propagators::DependentVariablesInterface< TIME_TYPE > > dependentVariablesInterface );

extern template std::vector< std::shared_ptr< ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >
createObservationSimulators< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationSettingsList,
        const simulation_setup::SystemOfBodies& bodies );

}  // namespace observation_models

namespace simulation_setup
{

extern template class OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE >;

extern template std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >
createStateDerivativePartials< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::unordered_map< propagators::IntegratedStateType,
                                  std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< STATE_SCALAR_TYPE, TIME_TYPE > > > >
                stateDerivativeModels,
        const SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< STATE_SCALAR_TYPE > > parametersToEstimate );

extern template std::shared_ptr< propagators::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >
createDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >(
        const SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< STATE_SCALAR_TYPE > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated );

extern template std::shared_ptr< propagators::VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE > >
createVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >(
        const SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< STATE_SCALAR_TYPE > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< STATE_SCALAR_TYPE > > parametersToEstimate,
        const bool integrateEquationsOnCreation );

extern template std::shared_ptr< EstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE > >
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

extern template std::shared_ptr< observation_models::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >
simulateObservations< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationsToSimulate,
        const std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >&
                observationSimulators,
        const SystemOfBodies bodies );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_ESTIMATION_H
