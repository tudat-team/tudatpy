/*    Copyright (c) 2010-2017, Delft University of Technology
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
    case spherical_harmonics_cosine_coefficient_block:
        isDoubleParameter = false;
        break;
    case spherical_harmonics_sine_coefficient_block:
        isDoubleParameter = false;
        break;
    case constant_rotation_rate:
        isDoubleParameter = true;
        break;
    case rotation_pole_position:
        isDoubleParameter = false;
        break;
    case constant_additive_observation_bias:
        isDoubleParameter = false;
        break;
    case constant_relative_observation_bias:
        isDoubleParameter = false;
        break;
    default:
        throw std::runtime_error( "Error, parameter type " + boost::lexical_cast< std::string >( parameterType ) +
                                  " not found when getting parameter type" );
    }
    return isDoubleParameter;
}

//! Function to determine whether the given (non-dynamical) parameter influences a body's orientation.
bool isParameterRotationMatrixProperty( const EstimatebleParametersEnum parameterType )
{
    bool flag;
    switch( parameterType )
    {
    case constant_rotation_rate:
        flag = true;
        break;
    case rotation_pole_position:
        flag = true;
        break;
    default:
        flag = false;
        break;
    }
    return flag;
}



}

}

