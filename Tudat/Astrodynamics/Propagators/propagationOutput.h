/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATIONOUTPUT_H
#define TUDAT_PROPAGATIONOUTPUT_H

#include <boost/function.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/Astrodynamics/Propagators/dynamicsStateDerivativeModel.h"
#include "Tudat/SimulationSetup/body.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluate a function with two input variables (by reference) from function pointers
/*!
 *  Function to evaluate a function with two input variables (by reference) from function pointers that return these
 *  two input variables.
 *  \param functionToEvaluate Function that is to be evaluated with input from function pointers.
 *  \param firstInput Function returning first input to functionToEvaluate.
 *  \param secondInput Function returning second input to functionToEvaluate.
 *  \return Output from functionToEvaluate, using functions firstInput and secondInput as input.
 */
template< typename OutputType, typename InputType >
OutputType evaluateReferenceFunction(
        const boost::function< OutputType( const InputType&, const InputType& ) > functionToEvaluate,
        const boost::function< InputType( ) > firstInput,
        const boost::function< InputType( ) > secondInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ) );
}

//! Function to evaluate a function with two input variables from function pointers
/*!
 *  Function to evaluate a function with two input variables from function pointers that return these
 *  two input variables.
 *  \param functionToEvaluate Function that is to be evaluated with input from function pointers.
 *  \param firstInput Function returning first input to functionToEvaluate.
 *  \param secondInput Function returning second input to functionToEvaluate.
 *  \return Output from functionToEvaluate, using functions firstInput and secondInput as input.
 */
template< typename OutputType, typename InputType >
OutputType evaluateFunction(
        const boost::function< OutputType( const InputType, const InputType ) > functionToEvaluate,
        const boost::function< InputType( ) > firstInput,
        const boost::function< InputType( ) > secondInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ) );
}

//! Funtion to get the size of a dependent variable
/*!
 * Funtion to get the size (i.e. number of values in variable: one for altitude, three for position, etc.)
 * of a dependent variable
 * \param dependentVariableSettings Dependent variable type for which size is to be determined.
 * \return Size of requested dependent variable.
 */
int getDependentVariableSize(
        const PropagationDependentVariables dependentVariableSettings );

//! Function to create a function returning a requested dependent variable value (of type double).
/*!
 *  Function to create a function returning a requested dependent variable value (of type double), retrieved from
 *  environment and/or state derivative models.
 *  \param dependentVariableSettings Settings for dependent variable that is to be returned by function created here.
 *  \param bodyMap List of bodies to use in simulations (containing full environment).
 *  \param stateDerivativeModels List of state derivative models used in simulations (sorted by dynamics type as key)
 *  \return Function returning requested dependent variable. NOTE: The environment and state derivative models need to
 *  be updated to current state and independent variable before computation is performed.
 */
template< typename TimeType = double, typename StateScalarType = double >
boost::function< double( ) > getDoubleDependentVariableFunction(
        const boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::unordered_map< IntegratedStateType,
        std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > > stateDerivativeModels =
        std::unordered_map< IntegratedStateType,
        std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >( ) )
{
    boost::function< double( ) > variableFunction;

    // Retrieve base information on dependent variable
    PropagationDependentVariables dependentVariable = dependentVariableSettings->variableType_;
    const std::string& bodyWithProperty = dependentVariableSettings->associatedBody_;
    const std::string& secondaryBody = dependentVariableSettings->secondaryBody_;

    // Check dependent variable type and create function accordingly.
    switch( dependentVariable )
    {
    case mach_number_dependent_variable:
    {
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {

        }

        boost::function< double( const double, const double ) > functionToEvaluate =
                boost::bind( &aerodynamics::computeMachNumber, _1, _2 );

        // Retrieve functions for airspeed and speed of sound.
        boost::function< double( ) > firstInput =
                boost::bind( &aerodynamics::FlightConditions::getCurrentAirspeed,
                             bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
        boost::function< double( ) > secondInput =
                boost::bind( &aerodynamics::FlightConditions::getCurrentSpeedOfSound,
                             bodyMap.at( bodyWithProperty )->getFlightConditions( ) );


        variableFunction = boost::bind( &evaluateFunction< double, double >,
                                        functionToEvaluate, firstInput, secondInput );
        break;
    }
    case altitude_dependent_variable:
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {

        }
        variableFunction = boost::bind( &aerodynamics::FlightConditions::getCurrentAltitude,
                                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
        break;
    case airspeed_dependent_variable:
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {

        }
        variableFunction = boost::bind( &aerodynamics::FlightConditions::getCurrentAirspeed,
                                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
        break;
    case local_density_dependent_variable:
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {

        }
        variableFunction = boost::bind( &aerodynamics::FlightConditions::getCurrentDensity,
                                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
        break;
    case radiation_pressure_dependent_variable:
        if( bodyMap.at( bodyWithProperty )->getRadiationPressureInterfaces( ).count( secondaryBody ) == 0 )
        {

        }
        variableFunction = boost::bind( &electro_magnetism::RadiationPressureInterface::getCurrentRadiationPressure,
                                        bodyMap.at( bodyWithProperty )->getRadiationPressureInterfaces( ).at( secondaryBody ) );
        break;
    case relative_distance_dependent_variable:
    {
        // Retrieve functions for positions of two bodies.
        boost::function< double( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                boost::bind( &linear_algebra::computeNormOfVectorDifference, _1, _2 );
        boost::function< Eigen::Vector3d( ) > firstInput =
                boost::bind( &simulation_setup::Body::getPosition, bodyMap.at( bodyWithProperty ) );
        boost::function< Eigen::Vector3d( ) > secondInput =
                boost::bind( &simulation_setup::Body::getPosition, bodyMap.at( secondaryBody ) );

        variableFunction = boost::bind(
                    &evaluateReferenceFunction< double, Eigen::Vector3d >, functionToEvaluate, firstInput, secondInput );
        break;
    }
    case relative_speed_dependent_variable:
    {
        // Retrieve functions for velicoty of two bodies.
        boost::function< double( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                boost::bind( &linear_algebra::computeNormOfVectorDifference, _1, _2 );
        boost::function< Eigen::Vector3d( ) > firstInput =
                boost::bind( &simulation_setup::Body::getVelocity, bodyMap.at( bodyWithProperty ) );
        boost::function< Eigen::Vector3d( ) > secondInput =
                boost::bind( &simulation_setup::Body::getVelocity, bodyMap.at( secondaryBody ) );

        variableFunction = boost::bind(
                    &evaluateReferenceFunction< double, Eigen::Vector3d >, functionToEvaluate, firstInput, secondInput );

        break;
    }
    case single_acceleration_norm_dependent_variable:
    {
        // Check input consistency
        boost::shared_ptr< SingleAccelerationDependentVariableSaveSettings > accelerationDependentVariableSettings =
                boost::dynamic_pointer_cast< SingleAccelerationDependentVariableSaveSettings >(
                    dependentVariableSettings );
        if( accelerationDependentVariableSettings == NULL )
        {
            std::string errorMessage = "Error, inconsistent inout when creating dependent variable function of type single_acceleration_norm_dependent_variable";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            // Retrieve list of suitable acceleration models (size should be one to avoid ambiguities)
            std::vector< boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                    listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                        accelerationDependentVariableSettings->associatedBody_,
                        accelerationDependentVariableSettings->secondaryBody_,
                        stateDerivativeModels, accelerationDependentVariableSettings->accelerationModeType_ );
            if( listOfSuitableAccelerationModels.size( ) != 1 )
            {
                std::string errorMessage = "Error when getting acceleration between bodies " +
                        accelerationDependentVariableSettings->associatedBody_ + " and " +
                        accelerationDependentVariableSettings->secondaryBody_ + " of type " +
                        boost::lexical_cast< std::string >(
                                            accelerationDependentVariableSettings->accelerationModeType_ ) +
                        ", no such acceleration found";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                boost::function< Eigen::Vector3d( ) > vectorFunction =
                        boost::bind( &basic_astrodynamics::AccelerationModel3d::getAcceleration,
                                     listOfSuitableAccelerationModels.at( 0 ) );
                variableFunction = boost::bind( &linear_algebra::getVectorNormFromFunction, vectorFunction );
            }
        }
        break;
    }
    case total_acceleration_norm_dependent_variable:
    {
        // Retrieve model responsible for computing accelerations of requested bodies.
        boost::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > nBodyModel =
                getTranslationalStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
        boost::function< Eigen::Vector3d( ) > vectorFunction =
                boost::bind( &NBodyStateDerivative< StateScalarType, TimeType >::getTotalAccelerationForBody,
                             nBodyModel, bodyWithProperty );
        variableFunction = boost::bind( &linear_algebra::getVectorNormFromFunction, vectorFunction );

        break;
    }
    default:
        std::string errorMessage =
                "Error, did not recognize double dependent variable type when making variable function: " +
                boost::lexical_cast< std::string >( dependentVariableSettings->variableType_ );
        throw std::runtime_error( errorMessage );
    }
    return variableFunction;
}

//! Function to create a function returning a requested dependent variable value (of type VectorXd).
/*!
 *  Function to create a function returning a requested dependent variable value (of type VectorXd), retrieved from
 *  environment and/or state derivative models.
 *  \param dependentVariableSettings Settings for dependent variable that is to be returned by function created here.
 *  \param bodyMap List of bodies to use in simulations (containing full environment).
 *  \param stateDerivativeModels List of state derivative models used in simulations (sorted by dynamics type as key)
 *  \return Function returning requested dependent variable. NOTE: The environment and state derivative models need to
 *  be updated to current state and independent variable before computation is performed.
 */
template< typename TimeType = double, typename StateScalarType = double >
std::pair< boost::function< Eigen::VectorXd( ) >, int > getVectorDependentVariableFunction(
        const boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::unordered_map< IntegratedStateType,
        std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > > stateDerivativeModels =
        std::unordered_map< IntegratedStateType,
        std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >( ) )
{
    boost::function< Eigen::VectorXd( ) > variableFunction;
    int parameterSize;

    // Retrieve base information on dependent variable
    PropagationDependentVariables dependentVariable = dependentVariableSettings->variableType_;
    const std::string& bodyWithProperty = dependentVariableSettings->associatedBody_;
    const std::string& secondaryBody = dependentVariableSettings->secondaryBody_;

    // Check dependent variable type and create function accordingly.
    switch( dependentVariable )
    {
    case relative_position_dependent_variable:
    {
        // Retrieve functions for positions of two bodies.
        boost::function< Eigen::Vector3d( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                boost::bind( &linear_algebra::computeVectorDifference, _1, _2 );
        boost::function< Eigen::Vector3d( ) > firstInput =
                boost::bind( &simulation_setup::Body::getPosition, bodyMap.at( bodyWithProperty ) );
        boost::function< Eigen::Vector3d( ) > secondInput =
                boost::bind( &simulation_setup::Body::getPosition, bodyMap.at( secondaryBody ) );

        variableFunction = boost::bind(
                    &evaluateReferenceFunction< Eigen::Vector3d, Eigen::Vector3d >,
                    functionToEvaluate, firstInput, secondInput );
        parameterSize = 3;
        break;
    }
    case relative_velocity_dependent_variable:
    {
        // Retrieve functions for velocities of two bodies.
        boost::function< Eigen::Vector3d( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                boost::bind( &linear_algebra::computeVectorDifference, _1, _2 );
        boost::function< Eigen::Vector3d( ) > firstInput =
                boost::bind( &simulation_setup::Body::getVelocity, bodyMap.at( bodyWithProperty ) );
        boost::function< Eigen::Vector3d( ) > secondInput =
                boost::bind( &simulation_setup::Body::getVelocity, bodyMap.at( secondaryBody ) );

        variableFunction = boost::bind(
                    &evaluateReferenceFunction< Eigen::Vector3d, Eigen::Vector3d >,
                    functionToEvaluate, firstInput, secondInput );
        parameterSize = 3;


        break;
    }
    case total_acceleration_dependent_variable:
    {
        // Retrieve model responsible for computing accelerations of requested bodies.
        boost::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > nBodyModel =
                getTranslationalStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
        variableFunction =
                boost::bind( &NBodyStateDerivative< StateScalarType, TimeType >::getTotalAccelerationForBody, nBodyModel,
                             bodyWithProperty );
        parameterSize = 3;


        break;
    }
    case single_acceleration_dependent_variable:
    {
        // Check input consistency.
        boost::shared_ptr< SingleAccelerationDependentVariableSaveSettings > accelerationDependentVariableSettings =
                boost::dynamic_pointer_cast< SingleAccelerationDependentVariableSaveSettings >( dependentVariableSettings );
        if( accelerationDependentVariableSettings == NULL )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type single_acceleration_dependent_variable";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            // Retrieve list of suitable acceleration models (size should be one to avoid ambiguities)
            std::vector< boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                    listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                        accelerationDependentVariableSettings->associatedBody_,
                        accelerationDependentVariableSettings->secondaryBody_,
                        stateDerivativeModels, accelerationDependentVariableSettings->accelerationModeType_ );
            if( listOfSuitableAccelerationModels.size( ) != 1 )
            {
                std::string errorMessage = "Error when getting acceleration between bodies " +
                        accelerationDependentVariableSettings->associatedBody_ + " and " +
                        accelerationDependentVariableSettings->secondaryBody_ + " of type " +
                        boost::lexical_cast< std::string >(
                                            accelerationDependentVariableSettings->accelerationModeType_ ) +
                        ", no such acceleration found";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                //boost::function< Eigen::Vector3d( ) > vectorFunction =
                variableFunction = boost::bind( &basic_astrodynamics::AccelerationModel3d::getAcceleration,
                                                listOfSuitableAccelerationModels.at( 0 ) );
                parameterSize = 3;
            }
        }
        break;
    }
    default:
        std::string errorMessage =
                "Error, did not recognize vector dependent variable type when making variable function: " +
                boost::lexical_cast< std::string >( dependentVariableSettings->variableType_ );
        throw std::runtime_error( errorMessage );
    }
    return std::make_pair( variableFunction, parameterSize );
}

//! Function to evaluate a set of double and vector-returning functions and concatenate the results.
/*!
 * Function to evaluate a set of double and vector-returning functions and concatenate the results. Results of double
 * function list are put in return vector first, followed by those in vector function list.
 * \param doubleFunctionList List of functions returning double variables
 * \param vectorFunctionList List of functions returning vector variables (pairs denote function and return vector size)
 * \param totalSize Total size of concatenated vector (used as input for efficiency.
 * \return Concatenated results from input functions.
 */
Eigen::VectorXd evaluateListOfFunctions(
        const std::vector< boost::function< double( ) > >& doubleFunctionList,
        const std::vector< std::pair< boost::function< Eigen::VectorXd( ) >, int > > vectorFunctionList,
        const int totalSize );

//! Function to create a function that evaluates a list of dependent variables and concatenates the results.
/*!
 *  Function to create a function that evaluates a list of dependent variables and concatenates the results.
 *  Dependent variables functions are created inside this function from a list of settings on their required
 *  types/properties.
 *  \param saveSettings Object containing types and other properties of dependent variables.
 *  \param bodyMap List of bodies to use in simulations (containing full environment).
 *  \param stateDerivativeModels List of state derivative models used in simulations (sorted by dynamics type as key)
 *  \return Function returning requested dependent variable values.
 *  NOTE: The environment and state derivative models need to
 *  be updated to current state and independent variable before computation is performed.
 */
template< typename TimeType = double, typename StateScalarType = double >
boost::function< Eigen::VectorXd( ) > createDependentVariableListFunction(
        const boost::shared_ptr< DependentVariableSaveSettings > saveSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::unordered_map< IntegratedStateType,
        std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels =
        std::unordered_map< IntegratedStateType,
        std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >( ) )
{
    // Retrieve list of save settings
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables =
            saveSettings->dependentVariables_;

    // create list of double and vector parameters
    std::vector< boost::function< double( ) > > doubleFunctionList;
    std::vector< std::pair< boost::function< Eigen::VectorXd( ) >, int > > vectorFunctionList;
    int totalVariableSize = 0;
    for( unsigned int i = 0; i < dependentVariables.size( ); i++ )
    {
        // Create double parameter
        if( getDependentVariableSize( dependentVariables.at( i )->variableType_ ) == 1 )
        {
            doubleFunctionList.push_back( getDoubleDependentVariableFunction(
                                              dependentVariables.at( i ),
                                              bodyMap, stateDerivativeModels ) );
            totalVariableSize++;
        }
        // Create vector parameter
        else
        {
            vectorFunctionList.push_back( getVectorDependentVariableFunction(
                                              dependentVariables.at( i ),
                                              bodyMap, stateDerivativeModels ) );
            totalVariableSize += vectorFunctionList.at( vectorFunctionList.size( ) - 1 ).second;
        }
    }

    // Create function conatenating function results.
    return boost::bind( &evaluateListOfFunctions, doubleFunctionList, vectorFunctionList, totalVariableSize );
}


} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATIONOUTPUT_H
