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
#include "tudatpy/kernelExternTemplatesPropagationOutput.h"
#include "tudatpy/kernelExternTemplatesRelativisticTimeConverter.h"
#include "tudat/simulation/environment_setup/createRelativisticTimeConverter.h"

namespace tudat
{

namespace propagators
{

template std::shared_ptr< EnvironmentUpdater< STATE_SCALAR_TYPE, TIME_TYPE > >
createEnvironmentUpdaterForDynamicalEquations< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::shared_ptr< SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > > propagatorSettings,
        const simulation_setup::SystemOfBodies& bodies );

template std::pair< std::function< Eigen::VectorXd( ) >, std::map< std::pair< int, int >, std::string > >
createDependentVariableListFunction< TIME_TYPE, STATE_SCALAR_TYPE >(
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables,
        const simulation_setup::SystemOfBodies& bodies,
        std::map< std::pair< int, int >, std::shared_ptr< SingleDependentVariableSaveSettings > >& orderedDependentVariables,
        const std::unordered_map< IntegratedStateType,
                                  std::vector< std::shared_ptr< SingleStateTypeDerivative< STATE_SCALAR_TYPE, TIME_TYPE > > > >&
                stateDerivativeModels,
        const std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >& stateDerivativePartials,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< STATE_SCALAR_TYPE > > parametersToEstimate );

template std::shared_ptr< PropagationTerminationCondition > createPropagationTerminationConditions< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const double initialTimeStep,
        const std::unordered_map< IntegratedStateType,
                                  std::vector< std::shared_ptr< SingleStateTypeDerivative< STATE_SCALAR_TYPE, TIME_TYPE > > > >&
                stateDerivativeModels,
        const std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >& stateDerivativePartials );

template std::shared_ptr< DampedInitialRotationalStateResults< STATE_SCALAR_TYPE > >
getZeroProperModeRotationalStateWithStruct< TIME_TYPE, STATE_SCALAR_TYPE >(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > > propagatorSettings,
        const double bodyMeanRotationRate,
        const std::vector< double > dissipationTimes,
        const bool propagateNominal );

}  // namespace propagators

namespace simulation_setup
{

template void setRelativisticTimeConverters< STATE_SCALAR_TYPE, TIME_TYPE >(
        const SystemOfBodies& bodies,
        const std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings< STATE_SCALAR_TYPE, TIME_TYPE > > >&
                converterSettings );

}  // namespace simulation_setup

}  // namespace tudat
