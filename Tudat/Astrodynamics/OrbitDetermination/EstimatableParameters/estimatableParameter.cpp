/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Function to determine whether the given parameter represents an initial dynamical state, or a static parameter.
bool isParameterDynamicalPropertyInitialState( const EstimatebleParametersEnum parameterType )
{
    bool flag;
    switch( parameterType )
    {
    case initial_body_state:
        flag = true;
        break;
    default:
        flag = false;
        break;
    }
    return flag;
}

//! Function to determine whether the given (non-dynamical) parameter is a double or vector parameter.
bool isDoubleParameter( const EstimatebleParametersEnum parameterType )
{
    bool isDoubleParameter;
    switch( parameterType )
    {
    case gravitational_parameter:
        isDoubleParameter = true;
        break;
    case constant_drag_coefficient:
        isDoubleParameter = true;
        break;
    case radiation_pressure_coefficient:
        isDoubleParameter = true;
        break;
    default:
        throw std::runtime_error( "Error, parameter type " + boost::lexical_cast< std::string >( parameterType ) +
                                  " not found when getting parameter type" );
    }
    return isDoubleParameter;
}




}

}

