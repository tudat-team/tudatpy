/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudatpy/kernelExternTemplatesObservationManagers.h"

namespace tudat
{

namespace observation_models
{

template std::shared_ptr< ObservationManagerBase< STATE_SCALAR_TYPE, TIME_TYPE > >
createObservationManager< 3, STATE_SCALAR_TYPE, TIME_TYPE >(
        const ObservableType observableType,
        const std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettingsList,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< STATE_SCALAR_TYPE > > parametersToEstimate,
        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const std::shared_ptr< propagators::DependentVariablesInterface< TIME_TYPE > > dependentVariablesInterface );

}  // namespace observation_models

}  // namespace tudat
