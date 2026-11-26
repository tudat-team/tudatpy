/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_SIMULATE_OBSERVATIONS_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_SIMULATE_OBSERVATIONS_H

// Explicit instantiation declaration of the observation simulation for the kernel flavour < STATE_SCALAR_TYPE, TIME_TYPE >,
// instantiated once in src/tudatpy/kernelInstantiationEstimationFunctions.cpp.

#include "tudatpy/scalarTypes.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tudat
{

namespace simulation_setup
{

extern template std::shared_ptr< observation_models::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >
simulateObservations< STATE_SCALAR_TYPE, TIME_TYPE >(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationsToSimulate,
        const std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >&
                observationSimulators,
        const SystemOfBodies bodies );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_SIMULATE_OBSERVATIONS_H
