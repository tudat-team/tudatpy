/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_VARIATIONALEQUATIONSSOLVERFORWARDDECLARATIONS_H
#define TUDAT_VARIATIONALEQUATIONSSOLVERFORWARDDECLARATIONS_H

#include <type_traits>

#include "tudat/basics/tudatTypeTraits.h"

namespace tudat
{

namespace propagators
{

template< typename StateScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< StateScalarType, TimeType >::value, int >::type Dummy >
class VariationalEquationsSolver;

template< typename StateScalarType, typename TimeType >
class SingleArcVariationalEquationsSolver;

template< typename StateScalarType, typename TimeType >
class MultiArcVariationalEquationsSolver;

template< typename StateScalarType, typename TimeType >
class HybridArcVariationalEquationsSolver;

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_VARIATIONALEQUATIONSSOLVERFORWARDDECLARATIONS_H
