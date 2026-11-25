/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_PROPAGATION_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_PROPAGATION_H

// Explicit instantiation declarations of the heavy propagation templates for the kernel flavour
// < STATE_SCALAR_TYPE, TIME_TYPE >. Each is instantiated once, in src/tudatpy/kernelInstantiation*.cpp.

#include "tudatpy/scalarTypes.h"
#include "tudat/astro/propagators/getZeroProperModeRotationalInitialState.h"
#include "tudatpy/kernelExternTemplatesRelativisticTimeConverter.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulatorBase.h"
#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"
#include "tudat/simulation/propagation_setup/multiArcDynamicsSimulator.h"
#include "tudat/simulation/propagation_setup/hybridArcDynamicsSimulator.h"

namespace tudat
{

namespace propagators
{

extern template class DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class MultiArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class HybridArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >;

extern template std::vector< std::shared_ptr< SingleStateTypeDerivative< STATE_SCALAR_TYPE, TIME_TYPE > > >
createStateDerivativeModels< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::shared_ptr< SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > > propagatorSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const TIME_TYPE propagationStartTime );

extern template std::shared_ptr< EnvironmentUpdater< STATE_SCALAR_TYPE, TIME_TYPE > >
createEnvironmentUpdaterForDynamicalEquations< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::shared_ptr< SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > > propagatorSettings,
        const simulation_setup::SystemOfBodies& bodies );

extern template std::pair< std::function< Eigen::VectorXd( ) >, std::map< std::pair< int, int >, std::string > >
createDependentVariableListFunction< TIME_TYPE, STATE_SCALAR_TYPE >(
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables,
        const simulation_setup::SystemOfBodies& bodies,
        std::map< std::pair< int, int >, std::shared_ptr< SingleDependentVariableSaveSettings > >& orderedDependentVariables,
        const std::unordered_map< IntegratedStateType,
                                  std::vector< std::shared_ptr< SingleStateTypeDerivative< STATE_SCALAR_TYPE, TIME_TYPE > > > >&
                stateDerivativeModels,
        const std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >& stateDerivativePartials,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< STATE_SCALAR_TYPE > > parametersToEstimate );

extern template std::shared_ptr< PropagationTerminationCondition > createPropagationTerminationConditions< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const double initialTimeStep,
        const std::unordered_map< IntegratedStateType,
                                  std::vector< std::shared_ptr< SingleStateTypeDerivative< STATE_SCALAR_TYPE, TIME_TYPE > > > >&
                stateDerivativeModels,
        const std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >& stateDerivativePartials );

extern template std::shared_ptr< DampedInitialRotationalStateResults< STATE_SCALAR_TYPE > >
getZeroProperModeRotationalStateWithStruct< TIME_TYPE, STATE_SCALAR_TYPE >(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > > propagatorSettings,
        const double bodyMeanRotationRate,
        const std::vector< double > dissipationTimes,
        const bool propagateNominal );

}  // namespace propagators

}  // namespace tudat

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_PROPAGATION_H
