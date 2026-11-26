/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_VARIATIONAL_EQUATIONS_SOLVERS_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_VARIATIONAL_EQUATIONS_SOLVERS_H

// Explicit instantiation declarations of the variational equations solvers and of the numerical simulator creation for
// the kernel flavour < STATE_SCALAR_TYPE, TIME_TYPE >, instantiated once in
// src/tudatpy/kernelInstantiation*VariationalEquationsSolver.cpp and kernelInstantiationEstimationFunctions.cpp.

#include "tudatpy/scalarTypes.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"

namespace tudat
{

namespace propagators
{

extern template class VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class MultiArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class HybridArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >;

}  // namespace propagators

namespace simulation_setup
{

extern template std::shared_ptr< propagators::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >
createDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >(
        const SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< STATE_SCALAR_TYPE > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated );

extern template std::shared_ptr< propagators::VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE > >
createVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >(
        const SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< STATE_SCALAR_TYPE > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< STATE_SCALAR_TYPE > > parametersToEstimate,
        const bool integrateEquationsOnCreation );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_VARIATIONAL_EQUATIONS_SOLVERS_H
