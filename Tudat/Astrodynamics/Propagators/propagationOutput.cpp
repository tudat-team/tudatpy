/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/Astrodynamics/Propagators/propagationOutput.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluate a set of double and vector-returning functions and concatenate the results.
Eigen::VectorXd evaluateListOfFunctions(
        const std::vector< boost::function< double( ) > >& doubleFunctionList,
        const std::vector< std::pair< boost::function< Eigen::VectorXd( ) >, int > > vectorFunctionList,
        const int totalSize)
{
    Eigen::VectorXd variableList = Eigen::VectorXd::Zero( totalSize );
    int currentIndex = 0;

    for( unsigned int i = 0; i < doubleFunctionList.size( ); i++ )
    {
        variableList( i ) = doubleFunctionList.at( i )( );
        currentIndex++;
    }

    for( unsigned int i = 0; i < vectorFunctionList.size( ); i++ )
    {
        variableList.segment( currentIndex, vectorFunctionList.at( i ).second ) =
                vectorFunctionList.at( i ).first( );
        currentIndex += vectorFunctionList.at( i ).second;
    }

    // Check consistency with input
    if( currentIndex != totalSize )
    {
        std::string errorMessage = "Error when evaluating lists of functions, sizes are inconsistent: " +
                boost::lexical_cast< std::string >( currentIndex ) + " and " +
                boost::lexical_cast< std::string >( totalSize );
        throw std::runtime_error( errorMessage );
    }

    return variableList;
}

//! Funtion to get the size of a dependent variable
int getDependentVariableSize(
        const PropagationDependentVariables dependentVariableSettings )
{
    int variableSize = -1;
    switch( dependentVariableSettings )
    {
    case mach_number_dependent_variable:
        variableSize = 1;
        break;
    case altitude_dependent_variable:
        variableSize = 1;
        break;
    case airspeed_dependent_variable:
        variableSize = 1;
        break;
    case local_density_dependent_variable:
        variableSize = 1;
        break;
    case relative_speed_dependent_variable:
        variableSize = 1;
        break;
    case relative_position_dependent_variable:
        variableSize = 3;
        break;
    case relative_distance_dependent_variable:
        variableSize = 1;
        break;
    case relative_velocity_dependent_variable:
        variableSize = 3;
        break;
    case radiation_pressure_dependent_variable:
        variableSize = 1;
        break;
    case total_acceleration_norm_dependent_variable:
        variableSize = 1;
        break;
    case single_acceleration_norm_dependent_variable:
        variableSize = 1;
        break;
    case total_acceleration_dependent_variable:
        variableSize = 3;
        break;
    case single_acceleration_dependent_variable:
        variableSize = 3;
        break;
    default:
        std::string errorMessage = "Error, did not recognize dependent variable size of type: " +
                boost::lexical_cast< std::string >( dependentVariableSettings );
        throw std::runtime_error( errorMessage );
    }
    return variableSize;
}


} // namespace propagators

} // namespace tudat
