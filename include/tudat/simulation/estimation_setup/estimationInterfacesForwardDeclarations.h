/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ESTIMATIONINTERFACESFORWARDDECLARATIONS_H
#define TUDAT_ESTIMATIONINTERFACESFORWARDDECLARATIONS_H

#include <type_traits>

#include "tudat/basics/tudatTypeTraits.h"

namespace tudat
{

namespace simulation_setup
{

class EstimationConvergenceChecker;

template< typename ObservationScalarType, typename TimeType >
class CovarianceAnalysisInput;

template< typename ObservationScalarType, typename TimeType >
class EstimationInput;

template< typename ObservationScalarType, typename TimeType >
struct CovarianceAnalysisOutput;

template< typename ObservationScalarType, typename TimeType >
struct EstimationOutput;

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
class OrbitDeterminationManager;

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_ESTIMATIONINTERFACESFORWARDDECLARATIONS_H
