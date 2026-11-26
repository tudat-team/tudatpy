/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_DYNAMICSSIMULATOREXTERNTEMPLATESTESTCASE_H
#define TUDAT_DYNAMICSSIMULATOREXTERNTEMPLATESTESTCASE_H

#include "tudat/basics/timeType.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"

// Explicit instantiations of the dynamics simulators for the state scalar and time types that the unit tests use next to
// < double, double >, compiled in tudat_estimation_test_case_support. Separate from estimationExternTemplatesTestCase.h,
// so that propagation tests do not parse the estimation headers.
#define TUDAT_TEST_DYNAMICS_SIMULATOR_INSTANTIATIONS( PREFIX, StateScalarType, TimeType )              \
    PREFIX template class tudat::propagators::DynamicsSimulator< StateScalarType, TimeType >;          \
    PREFIX template class tudat::propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType >; \
    PREFIX template class tudat::propagators::MultiArcDynamicsSimulator< StateScalarType, TimeType >;  \
    PREFIX template class tudat::propagators::HybridArcDynamicsSimulator< StateScalarType, TimeType >;

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS && TUDAT_BUILD_WITH_ESTIMATION_TOOLS
TUDAT_TEST_DYNAMICS_SIMULATOR_INSTANTIATIONS( extern, long double, tudat::Time )
#if TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS
TUDAT_TEST_DYNAMICS_SIMULATOR_INSTANTIATIONS( extern, long double, double )
#endif
TUDAT_TEST_DYNAMICS_SIMULATOR_INSTANTIATIONS( extern, double, tudat::Time )
#endif

#endif  // TUDAT_DYNAMICSSIMULATOREXTERNTEMPLATESTESTCASE_H
