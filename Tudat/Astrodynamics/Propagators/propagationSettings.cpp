/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"

namespace tudat
{

namespace propagators
{

//! Get size of state for single propagated state of given type.
int getSingleIntegrationSize( const IntegratedStateType stateType )
{
    int singleStateSize = 0;
    switch( stateType )
    {
    case transational_state:
        singleStateSize = 6;
        break;
    default:
        std::cerr<<"Did not recognize state type when getting size"<<std::endl;
    }
    return singleStateSize;
}

//! Get order of differential equation for governing equations of dynamics of given type.
int getSingleIntegrationDifferentialEquationOrder( const IntegratedStateType stateType )
{
    int singleStateSize = 0;
    switch( stateType )
    {
    case transational_state:
        singleStateSize = 2;
        break;
    default:
        std::cerr<<"Did not recognize state type when getting order"<<std::endl;
    }
    return singleStateSize;
}

} // namespace propagators

} // namespace tudat
