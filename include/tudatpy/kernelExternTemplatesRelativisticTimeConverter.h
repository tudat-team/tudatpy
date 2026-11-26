/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_RELATIVISTIC_TIME_CONVERTER_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_RELATIVISTIC_TIME_CONVERTER_H

// Explicit instantiation declaration of the relativistic time converter setup for the kernel flavour
// < STATE_SCALAR_TYPE, TIME_TYPE >, instantiated once in src/tudatpy/kernelInstantiationPropagationFunctions.cpp.

#include "tudatpy/scalarTypes.h"
#include "tudat/simulation/environment_setup/setRelativisticTimeConverters.h"

namespace tudat
{

namespace simulation_setup
{

extern template void setRelativisticTimeConverters< STATE_SCALAR_TYPE, TIME_TYPE >(
        const SystemOfBodies& bodies,
        const std::map< std::string, std::shared_ptr< DirectRelativisticTimeConverterSettings< STATE_SCALAR_TYPE, TIME_TYPE > > >&
                converterSettings );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_RELATIVISTIC_TIME_CONVERTER_H
