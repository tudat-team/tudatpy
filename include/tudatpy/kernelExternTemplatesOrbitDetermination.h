/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_ORBIT_DETERMINATION_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_ORBIT_DETERMINATION_H

// Explicit instantiation declarations of the orbit determination manager and the ephemeris fit for the kernel flavour
// < STATE_SCALAR_TYPE, TIME_TYPE >, instantiated once in src/tudatpy/kernelInstantiationOrbitDeterminationManager.cpp
// and kernelInstantiationEstimationFunctions.cpp.

#include "tudatpy/scalarTypes.h"
#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"

namespace tudat
{

namespace simulation_setup
{

extern template class OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE >;

extern template std::shared_ptr< EstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE > >
createBestFitToCurrentEphemeris< TIME_TYPE, STATE_SCALAR_TYPE >(
        const SystemOfBodies& bodies,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TIME_TYPE > >& integratorSettings,
        const TIME_TYPE initialTime,
        const TIME_TYPE finalTime,
        const TIME_TYPE dataPointInterval,
        const std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > additionalParameterNames,
        const int numberOfIterations,
        const bool reintegrateVariationalEquations,
        const double resultsPrintFrequency );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_ORBIT_DETERMINATION_H
