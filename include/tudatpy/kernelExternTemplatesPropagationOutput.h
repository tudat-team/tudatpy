/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_PROPAGATION_OUTPUT_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_PROPAGATION_OUTPUT_H

// Explicit instantiation declarations of the dependent variable and termination condition creation for the kernel
// flavour < STATE_SCALAR_TYPE, TIME_TYPE >, instantiated once in src/tudatpy/kernelInstantiationPropagationFunctions.cpp.

#include "tudatpy/scalarTypes.h"
#include "tudat/simulation/propagation_setup/propagationOutput.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"

namespace tudat
{

namespace propagators
{

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

}  // namespace propagators

}  // namespace tudat

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_PROPAGATION_OUTPUT_H
