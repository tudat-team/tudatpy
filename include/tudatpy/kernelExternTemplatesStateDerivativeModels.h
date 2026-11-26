/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_STATE_DERIVATIVE_MODELS_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_STATE_DERIVATIVE_MODELS_H

// Explicit instantiation declaration of the state derivative model creation for the kernel flavour
// < STATE_SCALAR_TYPE, TIME_TYPE >, instantiated once in src/tudatpy/kernelInstantiationStateDerivativeModels.cpp.

#include "tudatpy/scalarTypes.h"
#include "tudat/simulation/propagation_setup/createStateDerivativeModel.h"

namespace tudat
{

namespace propagators
{

extern template std::vector< std::shared_ptr< SingleStateTypeDerivative< STATE_SCALAR_TYPE, TIME_TYPE > > >
createStateDerivativeModels< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::shared_ptr< SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > > propagatorSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const TIME_TYPE propagationStartTime );

}  // namespace propagators

}  // namespace tudat

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_STATE_DERIVATIVE_MODELS_H
