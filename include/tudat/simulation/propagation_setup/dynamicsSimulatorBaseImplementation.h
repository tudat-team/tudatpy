/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DYNAMICSSIMULATORBASEIMPLEMENTATION_H
#define TUDAT_DYNAMICSSIMULATORBASEIMPLEMENTATION_H

#include "tudat/simulation/propagation_setup/dynamicsSimulatorBase.h"

namespace tudat
{

namespace propagators
{

template< typename StateScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< StateScalarType, TimeType >::value, int >::type Dummy >
bool DynamicsSimulator< StateScalarType, TimeType, Dummy >::getSetIntegratedResult( )
{
    return propagatorSettingsBase_->getOutputSettingsBase( )->getSetIntegratedResult( );
}

template< typename StateScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< StateScalarType, TimeType >::value, int >::type Dummy >
std::shared_ptr< DependentVariablesInterface< TimeType > >
DynamicsSimulator< StateScalarType, TimeType, Dummy >::getDependentVariablesInterface( )
{
    return getPropagationResults( )->getDependentVariablesInterface( );
}

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_DYNAMICSSIMULATORBASEIMPLEMENTATION_H
