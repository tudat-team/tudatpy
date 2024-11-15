/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include <iostream>
#include "tudat/astro/propagators/singleStateTypeDerivative.h"

namespace tudat
{

namespace propagators
{
std::string getIntegratedStateTypString( const IntegratedStateType stateType )
{
    std::string stateTypeString;
    switch( stateType )
    {
    case translational_state:
        stateTypeString = "Translational state";
        break;
    case rotational_state:
        stateTypeString = "Rotational state";
        break;
    case body_mass_state:
        stateTypeString = "Mass state";
        break;
    case custom_state:
        stateTypeString = "Custom state";
        break;
    case hybrid:
        stateTypeString = "Multi-type state";
        break;
    default:
        std::string errorMessage =
                "Did not recognize state type " + std::to_string( stateType ) + "when getting string";
       throw std::runtime_error( errorMessage );
    }
    return stateTypeString;
}


//! Get size of state for single propagated state of given type.
int getSingleIntegrationSize( const IntegratedStateType stateType )
{
    int singleStateSize = 0;
    switch( stateType )
    {
    case translational_state:
        singleStateSize = 6;
        break;
    case rotational_state:
        singleStateSize = 7;
        break;
    case body_mass_state:
        singleStateSize = 1;
        break;
    case custom_state:
        singleStateSize = 0;
        std::cerr<<"Warning when requesting state size of custom state, size is unknown. Returning value of 0"<<std::endl;
        break;
    default:
        std::string errorMessage =
                "Did not recognize state type " + std::to_string( stateType ) + "when getting size";
       throw std::runtime_error( errorMessage );
    }
    return singleStateSize;
}


//!Get order of differential equation for governing equations of dynamics of given type.
int getSingleIntegrationDifferentialEquationOrder( const IntegratedStateType stateType )
{
    int singleStateSize = 0;
    switch( stateType )
    {
    case translational_state:
        singleStateSize = 2;
        break;
    case body_mass_state:
        singleStateSize = 1;
        break;
    case rotational_state:
        singleStateSize = 1;
        break;
    default:
        std::string errorMessage =
                "Did not recognize state type " + std::to_string( stateType ) + "when getting order";
       throw std::runtime_error( errorMessage );
    }
    return singleStateSize;
}

//! Function to get the size of the generalized acceleration for a given state type
int getGeneralizedAccelerationSize( const IntegratedStateType stateType )
{
    int accelerationSize = 0;
    switch( stateType )
    {
    case translational_state:
        accelerationSize = 3;
        break;
    case body_mass_state:
        accelerationSize = 1;
        break;
    case rotational_state:
        accelerationSize = 3;
        break;
    default:
        std::string errorMessage =
                "Did not recognize state type " + std::to_string( stateType ) + "when getting acceleration sizw";
       throw std::runtime_error( errorMessage );
    }
    return accelerationSize;
}


//template class SingleStateTypeDerivative< double, double >;

}

}
