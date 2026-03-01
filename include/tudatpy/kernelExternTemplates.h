/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_H

#include "tudatpy/scalarTypes.h"

#if TUDAT_BUILD_KERNEL_EXTERN_TEMPLATES
#include "tudat/simulation/propagation_setup/dynamicsSimulatorForwardDeclarations.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolverForwardDeclarations.h"
#include "tudat/simulation/estimation_setup/estimationInterfacesForwardDeclarations.h"
#endif

#if TUDAT_BUILD_KERNEL_EXTERN_TEMPLATES

namespace tudat
{

namespace propagators
{

extern template class DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE, 0 >;
extern template class SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class MultiArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class HybridArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >;

extern template class VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE, 0 >;
extern template class SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class MultiArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >;
extern template class HybridArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >;

}  // namespace propagators

}  // namespace tudat

#endif

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_H
