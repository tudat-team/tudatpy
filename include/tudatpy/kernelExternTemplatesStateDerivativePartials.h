/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_STATE_DERIVATIVE_PARTIALS_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_STATE_DERIVATIVE_PARTIALS_H

// Explicit instantiation declaration of the state derivative partial creation for the kernel flavour
// < STATE_SCALAR_TYPE, TIME_TYPE >, instantiated once in src/tudatpy/kernelInstantiationStateDerivativePartials.cpp.

#include "tudatpy/scalarTypes.h"
#include "tudat/simulation/estimation_setup/createStateDerivativePartials.h"

namespace tudat
{

namespace simulation_setup
{

extern template std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >
createStateDerivativePartials< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::unordered_map< propagators::IntegratedStateType,
                                  std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< STATE_SCALAR_TYPE, TIME_TYPE > > > >
                stateDerivativeModels,
        const SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< STATE_SCALAR_TYPE > > parametersToEstimate );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_STATE_DERIVATIVE_PARTIALS_H
