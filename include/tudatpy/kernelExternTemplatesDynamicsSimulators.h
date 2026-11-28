/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_DYNAMICS_SIMULATORS_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_DYNAMICS_SIMULATORS_H

// Explicit instantiation declarations of the dynamics simulators for the kernel flavour < STATE_SCALAR_TYPE, TIME_TYPE >,
// instantiated once in src/tudatpy/kernelInstantiation*DynamicsSimulator.cpp and kernelInstantiationPropagationFunctions.cpp.

#include "tudatpy/scalarTypes.h"
#include "tudat/astro/propagators/getZeroProperModeRotationalInitialState.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulatorBase.h"
#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"
#include "tudat/simulation/propagation_setup/multiArcDynamicsSimulator.h"
#include "tudat/simulation/propagation_setup/hybridArcDynamicsSimulator.h"

namespace tudat
{

namespace propagators
{

extern template class DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class MultiArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class HybridArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >;

extern template std::shared_ptr< DampedInitialRotationalStateResults< STATE_SCALAR_TYPE > >
getZeroProperModeRotationalStateWithStruct< TIME_TYPE, STATE_SCALAR_TYPE >(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > > propagatorSettings,
        const double bodyMeanRotationRate,
        const std::vector< double > dissipationTimes,
        const bool propagateNominal );

}  // namespace propagators

}  // namespace tudat

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_DYNAMICS_SIMULATORS_H
