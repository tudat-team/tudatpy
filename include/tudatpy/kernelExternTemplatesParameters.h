/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_PARAMETERS_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_PARAMETERS_H

// Explicit instantiation declaration of the estimated-parameter creation for the kernel flavour
// < STATE_SCALAR_TYPE, TIME_TYPE >, instantiated once in src/tudatpy/kernelInstantiationParametersToEstimate.cpp.

#include "tudatpy/scalarTypes.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"

namespace tudat
{

namespace simulation_setup
{

extern template std::shared_ptr< estimatable_parameters::EstimatableParameterSet< STATE_SCALAR_TYPE > >
createParametersToEstimate< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >& parameterNames,
        const SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< STATE_SCALAR_TYPE > > propagatorSettings,
        const std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >& considerParameterNames,
        bool printParameterOrderWarning );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_PARAMETERS_H
