/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_OBSERVATION_MANAGERS_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_OBSERVATION_MANAGERS_H

// Explicit instantiation declarations of the observation manager creation for the kernel flavour
// < STATE_SCALAR_TYPE, TIME_TYPE >, instantiated once in src/tudatpy/kernelInstantiationObservationManager*.cpp and
// kernelInstantiationEstimationFunctions.cpp.

#include "tudatpy/scalarTypes.h"
#include "tudat/simulation/estimation_setup/createObservationManager.h"

namespace tudat
{

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

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_OBSERVATION_MANAGERS_H
