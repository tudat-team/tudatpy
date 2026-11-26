/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_OBSERVATION_SIMULATORS_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_OBSERVATION_SIMULATORS_H

// Explicit instantiation declaration of the observation simulator creation for the kernel flavour
// < STATE_SCALAR_TYPE, TIME_TYPE >, instantiated once in src/tudatpy/kernelInstantiationObservationSimulators.cpp.

#include "tudatpy/scalarTypes.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"

namespace tudat
{

namespace observation_models
{

extern template std::vector< std::shared_ptr< ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >
createObservationSimulators< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationSettingsList,
        const simulation_setup::SystemOfBodies& bodies );

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_OBSERVATION_SIMULATORS_H
