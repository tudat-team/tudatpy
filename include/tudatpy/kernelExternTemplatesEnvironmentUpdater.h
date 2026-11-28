/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_ENVIRONMENT_UPDATER_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_ENVIRONMENT_UPDATER_H

// Explicit instantiation declaration of the environment updater creation for the kernel flavour
// < STATE_SCALAR_TYPE, TIME_TYPE >, instantiated once in src/tudatpy/kernelInstantiationPropagationFunctions.cpp.

#include "tudatpy/scalarTypes.h"
#include "tudat/simulation/propagation_setup/createEnvironmentUpdater.h"

namespace tudat
{

namespace propagators
{

extern template std::shared_ptr< EnvironmentUpdater< STATE_SCALAR_TYPE, TIME_TYPE > >
createEnvironmentUpdaterForDynamicalEquations< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::shared_ptr< SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > > propagatorSettings,
        const simulation_setup::SystemOfBodies& bodies );

}  // namespace propagators

}  // namespace tudat

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_ENVIRONMENT_UPDATER_H
