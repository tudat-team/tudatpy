/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <functional>

#include "tudat/basics/utilities.h"
#include "tudat/astro/basic_astro/astrodynamicsFunctions.h"
#include "tudat/astro/aerodynamics/aerodynamics.h"
#include "tudat/astro/ephemerides/frameManager.h"
#include "tudat/astro/propagators/dynamicsStateDerivativeModel.h"
#include "tudat/astro/propagators/rotationalMotionStateDerivative.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/environment_setup/createFlightConditions.h"
#include "tudat/math/basic/rotationRepresentations.h"
#include "tudat/astro/aerodynamics/nrlmsise00Atmosphere.h"

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
OutputType evaluateBivariateReferenceFunction(
        const std::function< OutputType( const InputType&, const InputType& ) > functionToEvaluate,
        const std::function< InputType( ) > firstInput,
        const std::function< InputType( ) > secondInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ) );
}

template< typename OutputType, typename InputType >
OutputType evaluateReferenceFunction(
        const std::function< OutputType( const InputType& ) > functionToEvaluate,
        const std::function< InputType( ) > firstInput )
{
    return functionToEvaluate( firstInput( ) );
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
OutputType evaluateBivariateFunction(
        const std::function< OutputType( const InputType, const InputType ) > functionToEvaluate,
        const std::function< InputType( ) > firstInput,
        const std::function< InputType( ) > secondInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ) );
}

//! Function to evaluate a function with three input variables from function pointers
/*!
 *  Function to evaluate a function with three input variables from function pointers that return these
 *  three input variables.
 *  \param functionToEvaluate Function that is to be evaluated with input from function pointers.
 *  \param firstInput Function returning first input to functionToEvaluate.
 *  \param secondInput Function returning second input to functionToEvaluate.
 *  \param thirdInput Function returning third input to functionToEvaluate.
 *  \return Output from functionToEvaluate, using functions firstInput, secondInput and thirdInput as input.
 */
template< typename OutputType, typename FirstInputType, typename SecondInputType, typename ThirdInputType >
OutputType evaluateTrivariateFunction(
        const std::function< OutputType( const FirstInputType&, const SecondInputType, const ThirdInputType ) >
        functionToEvaluate,
        const std::function< FirstInputType( ) > firstInput,
        const std::function< SecondInputType( ) > secondInput,
        const std::function< ThirdInputType( ) > thirdInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ), thirdInput( ) );
}

Eigen::VectorXd getNormsOfAccelerationDifferencesFromLists(
                       const std::function< Eigen::VectorXd( ) > firstAccelerationFunction,
                       const std::function< Eigen::VectorXd( ) > secondAccelerationFunction );


//! Funtion to get the size of a dependent variable save settings
/*!
 * Funtion to get the size of a dependent variable save settings.
 * \param singleDependentVariableSaveSettings Save settings for a dependent variable.
 * \return Size of requested dependent variable to save (equal to the size of the associated dependent variable,
 * or equal to 1 if the property `component_` is set).
 */
int getDependentVariableSaveSize(
        const std::shared_ptr< SingleDependentVariableSaveSettings >& singleDependentVariableSaveSettings,
        const simulation_setup::SystemOfBodies& bodies );

//! Funtion to get the size of a dependent variable
/*!
 * Funtion to get the size (i.e. number of values in variable: one for altitude, three for position, etc.)
 * of a dependent variable
 * \param dependentVariableSettings Dependent variable type for which size is to be determined.
 * \return Size of requested dependent variable.
 */
int getDependentVariableSize(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const simulation_setup::SystemOfBodies& bodies );

std::pair< int, int > getDependentVariableShape(
    const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
    const simulation_setup::SystemOfBodies& bodies );

bool isScalarDependentVariable(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const simulation_setup::SystemOfBodies& bodies );

bool isMatrixDependentVariable(
    const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings );

//! Get the vector representation of a rotation matrix.
/*!
 *  Get the vector representation of a rotation matrix.
 *  \param currentRotationMatrix Rotation matrix that is to be put into vector rerpesentation
 *  \return Column vector consisting of transpose of concatenated rows of currentRotationMatrix input.
 */
Eigen::VectorXd getVectorRepresentationForRotationMatrix(
        const Eigen::Matrix3d& currentRotationMatrix );

//! Get the vector representation of a rotation matrix.
/*!
 *  Get the vector representation of a rotation matrix.
 *  \param rotationFunction Function returning the rotation matrix that is to be put into vector rerpesentation
 *  \return Column vector consisting of transpose of concatenated rows of rotationFunction input.
 */
Eigen::VectorXd getVectorRepresentationForRotationMatrixFunction(
        const std::function< Eigen::Matrix3d( ) > rotationFunction );

//! Get the vector representation of a quaternion.
/*!
 *  Get the vector representation of a quaternion. Quaternion is converted to a rotation matrix, which is then put into
 *  a vector representation.
 *  \param rotationFunction Function returning the quaternion that is to be put inot vector rerpesentation
 *  \return Column vector consisting of transpose of concatenated rows of matrix representation of rotationFunction input.
 */
Eigen::VectorXd getVectorRepresentationForRotationQuaternion(
        const std::function< Eigen::Quaterniond( ) > rotationFunction );


//! Get the 3x3 matrix representation from a vector with 9 entries
/*!
 *  Get the matrix representation from a vector with 9 entries. The vector is the transpose of the concatenated rows of
 *  the associated matrix.
 *  \param vectorRepresentation Vector representation of 3x3 matrix (transpose of the concatenated rows of
 *  the associated matrix)
 *  \return Matrix from 3x3 vector representation.
 */
Eigen::Matrix3d getMatrixFromVectorRotationRepresentation(
        const Eigen::VectorXd vectorRepresentation );

//! Get the quaternion formulation of an orthonormal matrix, from input of a vector with 9 entries corresponding to matrix
//! entries.
/*!
 *  et the quaternion formulation of an orthonormal matrix, from input of a vector with 9 entries corresponding to matrix
 *  entries. The input vector is the transpose of the concatenated rows of the associated (orthonormal) matrix.
 *  \param vectorRepresentation Vector representation of 3x3 matrix (transpose of the concatenated rows of
 *  the associated matrix)
 *  \return Quaternion representation of orthonormal matrix obtained from 3x3 vector representation.
 */
Eigen::Quaterniond getQuaternionFromVectorRotationRepresentation(
        const Eigen::VectorXd vectorRepresentation );

//! Function to convert a matrix to the format used to save dependent variables
/*!
 * Function to convert a matrix to the format used to save dependent variables
 * \param matrix Matrix that is to be converted
 * \param vector Vector storage format of matrix
 */
void getMatrixInOutputVectorRepresentation(
        const Eigen::MatrixXd& matrix, Eigen::VectorXd& vector );

//! Function to convert a vector dependent variable output to its original matrix representation
/*!
 *  Function to convert a vector dependent variable output to its original matrix representation
 *  \param vector Vector dependent variable output
 *  \param matrix Original matrix representation
 *  \param rows Number of rows in matrix output
 *  \param columns Number of columns in matrix output
 */
void getOutputVectorInMatrixRepresentation(
        const Eigen::VectorXd& vector, Eigen::MatrixXd& matrix,
        const int rows, const int columns );

//! Function to retrieve matrix block function output in vector representation
/*!
 *  Function to retrieve matrix block function output in vector representation
 * \param blockFunction Function that returns (by reference) a matrix block
 * \param numberOfRows Number of rows in matrix output
 * \param numberOfColumns Number of columns in matrix output
 * \return Block-matrix in vector representation
 */
Eigen::VectorXd getVectorFunctionFromBlockFunction(
        const std::function< void( Eigen::Block< Eigen::MatrixXd > ) > blockFunction,
        const int numberOfRows, const int numberOfColumns );

//! Function to compute the Fay-Riddell equilibrium heat flux from body properties
/*!
 * Function to compute the Fay-Riddell equilibrium heat flux from body properties
 * \param flightConditions Object describing the current atmospheric flight conditions of the vehicle
 * \param vehicleSystems Object describing the physical properties of the vehicle
 * \return Equilibrium heat flux according to Fay-Riddell model
 */
double computeEquilibriumFayRiddellHeatFluxFromProperties(
        const std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions,
        const std::shared_ptr< system_models::VehicleSystems > vehicleSystems );

//! Function to retrieve relevant spherical harmonic acceleration model for dependent variable setting
/*!
 *  Function to retrieve relevant spherical harmonic acceleration model for dependent variable setting
 *  \param dependentVariableSettings Settings for dependent variable, associatedBody_ defines body undergoing acceleration,
 *  secondaryBody_ body exerting acceleration
 *  \param stateDerivativeModels List of state derivative models from which acceleration is to be retrieved
 *  \return Relevant spherical harmonic acceleration model for dependent variable setting
 */
template< typename StateScalarType, typename TimeType >
std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
getSphericalHarmonicAccelerationForDependentVariables(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels,
        const bool allowThirdBodyAcceleration = false )
{
    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > selectedAccelerationModel;

    // Retrieve list of suitable acceleration models (size should be one to avoid ambiguities)s
    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
            listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                dependentVariableSettings->associatedBody_,
                dependentVariableSettings->secondaryBody_, stateDerivativeModels,
                basic_astrodynamics::spherical_harmonic_gravity );


    // Check if third-body counterpart of acceleration is found
    if( listOfSuitableAccelerationModels.size( ) == 0 )
    {
        listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                    dependentVariableSettings->associatedBody_,
                    dependentVariableSettings->secondaryBody_,
                    stateDerivativeModels, basic_astrodynamics::getAssociatedThirdBodyAcceleration(
                        basic_astrodynamics::spherical_harmonic_gravity ) );
    }

    if( listOfSuitableAccelerationModels.size( ) != 1 )
    {
        std::string errorMessage = "Error when getting spherical harmonic acceleration components between bodies " +
                dependentVariableSettings->associatedBody_ + " and " +
                dependentVariableSettings->secondaryBody_ + " of type " +
                std::to_string(
                    basic_astrodynamics::spherical_harmonic_gravity ) +
                ", no such acceleration found";
        throw std::runtime_error( errorMessage );
    }
    else
    {
        std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration =
                std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                    listOfSuitableAccelerationModels.at( 0 ) );


        if( sphericalHarmonicAcceleration != nullptr )
        {
            selectedAccelerationModel = sphericalHarmonicAcceleration;
        }
        else if( !allowThirdBodyAcceleration )
        {
            std::string errorMessage = "Error when getting spherical harmonic acceleration for dependent variable " +
                    dependentVariableSettings->associatedBody_ + " and " +
                    dependentVariableSettings->secondaryBody_ + " type is ionconsistent (third body not checked)";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::shared_ptr< gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel >
                    thirdBodySphericalHarmonicAcceleration =
                    std::dynamic_pointer_cast< gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                        listOfSuitableAccelerationModels.at( 0 ) );

            if( thirdBodySphericalHarmonicAcceleration == nullptr )
            {
                std::string errorMessage = "Error when getting spherical harmonic acceleration for dependent variable " +
                        dependentVariableSettings->associatedBody_ + " and " +
                        dependentVariableSettings->secondaryBody_ + " type is ionconsistent (third body checked)";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                selectedAccelerationModel = thirdBodySphericalHarmonicAcceleration;
            }
        }
    }

    return selectedAccelerationModel;
}

//! Function to retrieve the gravitational acceleration model for dependent variable settings.
/*!
 *  Function to retrieve the gravitational acceleration model for dependent variable settings. Returns the gravitational
 *  acceleration model acting between the two relevant bodies.
 *  \param dependentVariableSettings Settings for dependent variable, associatedBody_ defines body undergoing acceleration,
 *  secondaryBody_ body exerting acceleration
 *  \param stateDerivativeModels List of state derivative models from which acceleration is to be retrieved
 *  \param selectedAccelerationModelType Acceleration type of the returned model
 *  \return Relevant gravitational acceleration model for dependent variable setting
 */
template< typename StateScalarType, typename TimeType >
std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
getGravitationalAccelerationForDependentVariables(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels,
        basic_astrodynamics::AvailableAcceleration& selectedAccelerationModelType )
{

    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > selectedAccelerationModel;

    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
            listOfSuitableAccelerationModels;

    for( basic_astrodynamics::AvailableAcceleration accelerationModelType:
        {basic_astrodynamics::spherical_harmonic_gravity, basic_astrodynamics::polyhedron_gravity,
         basic_astrodynamics::point_mass_gravity, basic_astrodynamics::mutual_spherical_harmonic_gravity} )
    {
        listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                dependentVariableSettings->associatedBody_, dependentVariableSettings->secondaryBody_,
                stateDerivativeModels, accelerationModelType );

        // Check if gravitational acceleration between the two bodies was found
        if( listOfSuitableAccelerationModels.size( ) == 1 )
        {
            selectedAccelerationModelType = accelerationModelType;
            selectedAccelerationModel = listOfSuitableAccelerationModels.at( 0 );
            break;
        }

        // Check if third-body counterpart of acceleration is found
        if( listOfSuitableAccelerationModels.size( ) == 0 )
        {
            listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                    dependentVariableSettings->associatedBody_, dependentVariableSettings->secondaryBody_,
                    stateDerivativeModels, basic_astrodynamics::getAssociatedThirdBodyAcceleration( accelerationModelType ) );

            // Check if gravitational acceleration between the two bodies was found
            if( listOfSuitableAccelerationModels.size( ) == 1 )
            {
                selectedAccelerationModelType = basic_astrodynamics::getAssociatedThirdBodyAcceleration( accelerationModelType );
                selectedAccelerationModel = listOfSuitableAccelerationModels.at( 0 );
                break;
            }
        }
    }

    if( listOfSuitableAccelerationModels.size( ) != 1 )
    {
        std::string errorMessage = "Error when getting acceleration between bodies " +
                dependentVariableSettings->associatedBody_ + " and " +
                dependentVariableSettings->secondaryBody_ + ": no gravitational acceleration found.";
        throw std::runtime_error( errorMessage );
    }

    return selectedAccelerationModel;
}

Eigen::VectorXd getConstellationMinimumDistance(
    const std::function< Eigen::Vector3d( ) >& mainBodyPositionFunction,
    const std::vector< std::function< Eigen::Vector3d( ) > >& bodiesToCheckPositionFunctions );

Eigen::VectorXd getConstellationMinimumVisibleDistance(
    const std::function< Eigen::Vector3d( ) >& mainBodyPositionFunction,
    const std::vector< std::function< Eigen::Vector3d( ) > >& bodiesToCheckPositionFunctions,
    const std::shared_ptr< ground_stations::PointingAnglesCalculator > stationPointingAngleCalculator,
    const double limitAngle,
    const double time );

//! Function to create a function returning a requested dependent variable value (of type VectorXd).
/*!
 *  Function to create a function returning a requested dependent variable value (of type VectorXd), retrieved from
 *  environment and/or state derivative models.
 *  \param dependentVariableSettings Settings for dependent variable that is to be returned by function created here.
 *  \param bodies List of bodies to use in simulations (containing full environment).
 *  \param stateDerivativeModels List of state derivative models used in simulations (sorted by dynamics type as key).
 *  \param stateDerivativePartials List of state derivative partials used in simulations (sorted by dynamics type as key).
 *  \return Function returning requested dependent variable. NOTE: The environment and state derivative models need to
 *  be updated to current state and independent variable before computation is performed.
 */
template< typename TimeType = double, typename StateScalarType = double >
std::pair< std::function< Eigen::VectorXd( ) >, int > getVectorDependentVariableFunction(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels =
        std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >( ),
        const std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >& stateDerivativePartials =
        std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >( ) )
{
    std::function< Eigen::VectorXd( ) > variableFunction;
    int parameterSize;

    // Retrieve base information on dependent variable
    PropagationDependentVariables dependentVariable = dependentVariableSettings->dependentVariableType_;
    const std::string& bodyWithProperty = dependentVariableSettings->associatedBody_;
    const std::string& secondaryBody = dependentVariableSettings->secondaryBody_;

    // Check dependent variable type and create function accordingly.
    switch( dependentVariable )
    {
    case relative_position_dependent_variable:
    {
        // Retrieve functions for positions of two bodies.
        std::function< Eigen::Vector3d( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                std::bind( &linear_algebra::computeVectorDifference< 3 >, std::placeholders::_1, std::placeholders::_2 );
        std::function< Eigen::Vector3d( ) > firstInput =
                std::bind( &simulation_setup::Body::getPosition, bodies.at( bodyWithProperty ) );

        std::function< Eigen::Vector3d( ) > secondInput;
        if( secondaryBody != "SSB" )
        {
            secondInput = std::bind( &simulation_setup::Body::getPosition, bodies.at( secondaryBody ) );
        }
        else if( simulation_setup::getGlobalFrameOrigin( bodies ) == "SSB" )
        {
            secondInput = []( ){ return Eigen::Vector3d::Zero( ); };
        }
        else
        {
            throw std::runtime_error( "Error, requested state of " + bodyWithProperty + " w.r.t. SSB, but SSB is not frame origin" );
        }
        variableFunction = std::bind(
                    &evaluateBivariateReferenceFunction< Eigen::Vector3d, Eigen::Vector3d >,
                    functionToEvaluate, firstInput, secondInput );
        parameterSize = 3;
        break;
    }
    case relative_velocity_dependent_variable:
    {
        // Retrieve functions for velocities of two bodies.
        std::function< Eigen::Vector3d( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                std::bind( &linear_algebra::computeVectorDifference< 3 >, std::placeholders::_1, std::placeholders::_2 );
        std::function< Eigen::Vector3d( ) > firstInput =
                std::bind( &simulation_setup::Body::getVelocity, bodies.at( bodyWithProperty ) );

        std::function< Eigen::Vector3d( ) > secondInput;
        if( secondaryBody != "SSB" )
        {
            secondInput = std::bind( &simulation_setup::Body::getVelocity, bodies.at( secondaryBody ) );
        }
        else if( simulation_setup::getGlobalFrameOrigin( bodies ) == "SSB" )
        {
            secondInput = []( ){ return Eigen::Vector3d::Zero( ); };
        }
        else
        {
            throw std::runtime_error( "Error, requested state of " + bodyWithProperty + " w.r.t. SSB, but SSB is not frame origin" );
        }

        variableFunction = std::bind(
                    &evaluateBivariateReferenceFunction< Eigen::Vector3d, Eigen::Vector3d >,
                    functionToEvaluate, firstInput, secondInput );
        parameterSize = 3;


        break;
    }
    case total_acceleration_dependent_variable:
    {
        // Retrieve model responsible for computing accelerations of requested bodies.
        std::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > nBodyModel =
                getTranslationalStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
        nBodyModel->setUpdateRemovedAcceleration( dependentVariableSettings->associatedBody_ );
        variableFunction =
                std::bind( &NBodyStateDerivative< StateScalarType, TimeType >::getTotalAccelerationForBody, nBodyModel,
                           bodyWithProperty );
        parameterSize = 3;


        break;
    }
    case single_acceleration_dependent_variable:
    {
        // Check input consistency.
        std::shared_ptr< SingleAccelerationDependentVariableSaveSettings > accelerationDependentVariableSettings =
                std::dynamic_pointer_cast< SingleAccelerationDependentVariableSaveSettings >( dependentVariableSettings );
        if( accelerationDependentVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type single_acceleration_dependent_variable";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            // Retrieve list of suitable acceleration models (size should be one to avoid ambiguities)
            std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                    listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                        accelerationDependentVariableSettings->associatedBody_,
                        accelerationDependentVariableSettings->secondaryBody_,
                        stateDerivativeModels, accelerationDependentVariableSettings->accelerationModelType_ );

            // Check if third-body counterpart of acceleration is found
            if( listOfSuitableAccelerationModels.size( ) == 0 && basic_astrodynamics::isAccelerationDirectGravitational(
                        accelerationDependentVariableSettings->accelerationModelType_ ) )
            {
                listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                            accelerationDependentVariableSettings->associatedBody_,
                            accelerationDependentVariableSettings->secondaryBody_,
                            stateDerivativeModels, basic_astrodynamics::getAssociatedThirdBodyAcceleration(
                                accelerationDependentVariableSettings->accelerationModelType_  ) );
            }

            if( listOfSuitableAccelerationModels.size( ) != 1 )
            {
                std::string errorMessage = "Error when getting acceleration between bodies " +
                        accelerationDependentVariableSettings->associatedBody_ + " and " +
                        accelerationDependentVariableSettings->secondaryBody_ + " of type " +
                        getAccelerationModelName(
                            accelerationDependentVariableSettings->accelerationModelType_ ) +
                        ", no such acceleration found";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                std::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > nBodyModel =
                        getTranslationalStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
                std::map< std::string, std::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d > > removedAcceleration =
                        nBodyModel->getRemovedCentralAcceleration( );
                if( removedAcceleration.count( accelerationDependentVariableSettings->associatedBody_ ) > 0 )
                {
                    if( listOfSuitableAccelerationModels.at( 0 ) ==
                            removedAcceleration.at( accelerationDependentVariableSettings->associatedBody_ ) )
                    {
                        nBodyModel->setUpdateRemovedAcceleration( accelerationDependentVariableSettings->associatedBody_ );
                    }
                }

                variableFunction = std::bind( &basic_astrodynamics::AccelerationModel3d::getAcceleration,
                                              listOfSuitableAccelerationModels.at( 0 ) );
                parameterSize = 3;
            }
        }
        break;
    }
    case spherical_harmonic_acceleration_norm_terms_dependent_variable:
    {
        // Check input consistency.
        std::shared_ptr< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings > accelerationComponentVariableSettings =
                std::dynamic_pointer_cast< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >( dependentVariableSettings );
        if( accelerationComponentVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type spherical_harmonic_acceleration_norm_terms_dependent_variable";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::shared_ptr< basic_astrodynamics::AccelerationModel3d > sphericalHarmonicAcceleration =
                    getSphericalHarmonicAccelerationForDependentVariables(
                        accelerationComponentVariableSettings, stateDerivativeModels, true );

            if( std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                        sphericalHarmonicAcceleration ) != nullptr )
            {
                std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > directSphericalHarmonicAcceleration =
                        std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                            sphericalHarmonicAcceleration );

                directSphericalHarmonicAcceleration->setSaveSphericalHarmonicTermsSeparately( true );
                variableFunction = std::bind(
                            &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getConcatenatedAccelerationComponentNorms,
                            directSphericalHarmonicAcceleration, accelerationComponentVariableSettings->componentIndices_ );
            }
            else if( std::dynamic_pointer_cast< gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                         sphericalHarmonicAcceleration ) != nullptr )
            {
                std::shared_ptr< gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel > thirdBodySphericalHarmonicAcceleration =
                        std::dynamic_pointer_cast< gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                            sphericalHarmonicAcceleration );
                thirdBodySphericalHarmonicAcceleration->getAccelerationModelForBodyUndergoingAcceleration( )->
                        setSaveSphericalHarmonicTermsSeparately( true );
                thirdBodySphericalHarmonicAcceleration->getAccelerationModelForCentralBody( )->
                        setSaveSphericalHarmonicTermsSeparately( true );

                std::function< Eigen::VectorXd( ) > directAccelerationsFunction = std::bind(
                            &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getConcatenatedAccelerationComponents,
                            thirdBodySphericalHarmonicAcceleration->getAccelerationModelForBodyUndergoingAcceleration( ),
                            accelerationComponentVariableSettings->componentIndices_ );
                std::function< Eigen::VectorXd( ) > centralBodyAccelerationsFunction = std::bind(
                            &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getConcatenatedAccelerationComponents,
                            thirdBodySphericalHarmonicAcceleration->getAccelerationModelForCentralBody( ),
                            accelerationComponentVariableSettings->componentIndices_ );
                variableFunction =
                        std::bind( &getNormsOfAccelerationDifferencesFromLists,
                                   directAccelerationsFunction, centralBodyAccelerationsFunction );
            }
            else
            {
                throw std::runtime_error( "Error when getting spherical_harmonic_acceleration_norm_terms_dependent_variable, did not recognize acceleration model." );
            }

            parameterSize = accelerationComponentVariableSettings->componentIndices_.size( );
        }
        break;
    }
    case spherical_harmonic_acceleration_terms_dependent_variable:
    {
        // Check input consistency.
        std::shared_ptr< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings > accelerationComponentVariableSettings =
                std::dynamic_pointer_cast< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >( dependentVariableSettings );
        if( accelerationComponentVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type spherical_harmonic_acceleration_terms_dependent_variable";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::shared_ptr< basic_astrodynamics::AccelerationModel3d > sphericalHarmonicAcceleration =
                    getSphericalHarmonicAccelerationForDependentVariables(
                        accelerationComponentVariableSettings, stateDerivativeModels, true );

            if( std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                        sphericalHarmonicAcceleration ) != nullptr )
            {
                std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > directSphericalHarmonicAcceleration =
                        std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                            sphericalHarmonicAcceleration );

                directSphericalHarmonicAcceleration->setSaveSphericalHarmonicTermsSeparately( true );
                variableFunction = std::bind(
                            &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getConcatenatedAccelerationComponents,
                            directSphericalHarmonicAcceleration, accelerationComponentVariableSettings->componentIndices_ );
            }
            else if( std::dynamic_pointer_cast< gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                         sphericalHarmonicAcceleration ) != nullptr )
            {
                std::shared_ptr< gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel > thirdBodySphericalHarmonicAcceleration =
                        std::dynamic_pointer_cast< gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                            sphericalHarmonicAcceleration );
                thirdBodySphericalHarmonicAcceleration->getAccelerationModelForBodyUndergoingAcceleration( )->
                        setSaveSphericalHarmonicTermsSeparately( true );
                thirdBodySphericalHarmonicAcceleration->getAccelerationModelForCentralBody( )->
                        setSaveSphericalHarmonicTermsSeparately( true );

                std::function< Eigen::VectorXd( ) > directAccelerationsFunction = std::bind(
                            &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getConcatenatedAccelerationComponents,
                            thirdBodySphericalHarmonicAcceleration->getAccelerationModelForBodyUndergoingAcceleration( ),
                            accelerationComponentVariableSettings->componentIndices_ );
                std::function< Eigen::VectorXd( ) > centralBodyAccelerationsFunction = std::bind(
                            &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getConcatenatedAccelerationComponents,
                            thirdBodySphericalHarmonicAcceleration->getAccelerationModelForCentralBody( ),
                            accelerationComponentVariableSettings->componentIndices_ );
                variableFunction =
                        std::bind( &utilities::subtractFunctionReturn< Eigen::VectorXd >,
                                   directAccelerationsFunction, centralBodyAccelerationsFunction );
            }
            else
            {
                throw std::runtime_error( "Error when getting spherical_harmonic_acceleration_terms_dependent_variable, did not recognize acceleration model." );
            }
            parameterSize = 3 * accelerationComponentVariableSettings->componentIndices_.size( );

        }
        break;
    }
    case total_gravity_field_variation_acceleration:
    {

        std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration =
                std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                    getSphericalHarmonicAccelerationForDependentVariables(
                        dependentVariableSettings, stateDerivativeModels, false ) );

        std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > timeDependentGravityField =
                std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                    bodies.at( dependentVariableSettings->secondaryBody_ )->getGravityFieldModel( ) );

        if( timeDependentGravityField == nullptr )
        {
            throw std::runtime_error( "Error when requesting save of gravity field variation acceleration, central body " +
                                      dependentVariableSettings->secondaryBody_ +
                                      " has no TimeDependentSphericalHarmonicsGravityField." );
        }
        else
        {
            std::function< Eigen::VectorXd( const Eigen::MatrixXd&, const Eigen::MatrixXd& ) > accelerationFunction =
                    std::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getAccelerationWithAlternativeCoefficients,
                               sphericalHarmonicAcceleration, std::placeholders::_1, std::placeholders::_2 );
            std::function< Eigen::MatrixXd( ) > cosineCorrectionFunction =
                    std::bind( &gravitation::TimeDependentSphericalHarmonicsGravityField::getTotalCosineCoefficientCorrection,
                               timeDependentGravityField,
                               sphericalHarmonicAcceleration->getMaximumDegree( ),
                               sphericalHarmonicAcceleration->getMaximumOrder( ) );
            std::function< Eigen::MatrixXd( ) > sineCorrectionFunction =
                    std::bind( &gravitation::TimeDependentSphericalHarmonicsGravityField::getTotalSineCoefficientCorrection,
                               timeDependentGravityField,
                               sphericalHarmonicAcceleration->getMaximumDegree( ),
                               sphericalHarmonicAcceleration->getMaximumOrder( ) );

            variableFunction = std::bind( evaluateBivariateReferenceFunction< Eigen::VectorXd, Eigen::MatrixXd >,
                                          accelerationFunction, cosineCorrectionFunction, sineCorrectionFunction );

            parameterSize = 3;
        }

        break;
    }
    case total_spherical_harmonic_cosine_coefficient_variation:
    {
        // Check input consistency.
        std::shared_ptr< TotalGravityFieldVariationSettings > totalGravityFieldVariationSettings =
                std::dynamic_pointer_cast< TotalGravityFieldVariationSettings >( dependentVariableSettings );
        if( totalGravityFieldVariationSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type total_spherical_harmonic_cosine_coefficient_variation";
            throw std::runtime_error( errorMessage );
        }
        else
        {

            std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > timeDependentGravityField =
                    std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                        bodies.at( dependentVariableSettings->associatedBody_ )->getGravityFieldModel( ) );

            if( timeDependentGravityField == nullptr )
            {
                throw std::runtime_error( "Error when requesting save of gravity field variation, central body " +
                                          dependentVariableSettings->secondaryBody_ +
                                          " has no TimeDependentSphericalHarmonicsGravityField." );
            }
            else
            {
                std::vector< std::pair< int, int > > componentIndices = totalGravityFieldVariationSettings->componentIndices_;
                unsigned int numberOfCoefficients =  componentIndices.size( );

                variableFunction = [=]( )
                {
                    Eigen::VectorXd coefficientCorrections =
                            Eigen::VectorXd::Zero( numberOfCoefficients );
                    for( unsigned int i = 0; i < numberOfCoefficients; i++ )
                    {
                        coefficientCorrections( i ) = timeDependentGravityField->getSingleCosineCoefficientCorrection(
                                    componentIndices.at( i ).first, componentIndices.at( i ).second );
                    }
                    return coefficientCorrections;
                };

                parameterSize = numberOfCoefficients;
            }

        }
        break;

    }
    case total_spherical_harmonic_sine_coefficient_variation:
    {
        // Check input consistency.
        std::shared_ptr< TotalGravityFieldVariationSettings > totalGravityFieldVariationSettings =
                std::dynamic_pointer_cast< TotalGravityFieldVariationSettings >( dependentVariableSettings );
        if( totalGravityFieldVariationSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type total_spherical_harmonic_sine_coefficient_variation";
            throw std::runtime_error( errorMessage );
        }
        else
        {

            std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > timeDependentGravityField =
                    std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                        bodies.at( dependentVariableSettings->associatedBody_ )->getGravityFieldModel( ) );

            if( timeDependentGravityField == nullptr )
            {
                throw std::runtime_error( "Error when requesting save of gravity field variation, central body " +
                                          dependentVariableSettings->secondaryBody_ +
                                          " has no TimeDependentSphericalHarmonicsGravityField." );
            }
            else
            {
                std::vector< std::pair< int, int > > componentIndices = totalGravityFieldVariationSettings->componentIndices_;
                int numberOfCoefficients =  componentIndices.size( );

                variableFunction = [=]( )
                {
                    Eigen::VectorXd coefficientCorrections =
                            Eigen::VectorXd::Zero( numberOfCoefficients );
                    for( unsigned int i = 0; i < componentIndices.size( ); i++ )
                    {
                        coefficientCorrections( i ) = timeDependentGravityField->getSingleSineCoefficientCorrection(
                                    componentIndices.at( i ).first, componentIndices.at( i ).second );
                    }
                    return coefficientCorrections;
                };

                parameterSize = numberOfCoefficients;
            }

        }
        break;

    }
    case single_gravity_field_variation_acceleration:
    {
        std::shared_ptr< SingleVariationSphericalHarmonicAccelerationSaveSettings > accelerationVariableSettings =
                std::dynamic_pointer_cast< SingleVariationSphericalHarmonicAccelerationSaveSettings >( dependentVariableSettings );
        if( accelerationVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type single_gravity_field_variation_acceleration";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration =
                    std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                        getSphericalHarmonicAccelerationForDependentVariables(
                            accelerationVariableSettings, stateDerivativeModels, false ) );

            std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > timeDependentGravityField =
                    std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                        bodies.at( dependentVariableSettings->secondaryBody_ )->getGravityFieldModel( ) );

            if( timeDependentGravityField == nullptr )
            {
                throw std::runtime_error( "Error when requesting save of gravity field variation acceleration, central body " +
                                          dependentVariableSettings->secondaryBody_ +
                                          " has no TimeDependentSphericalHarmonicsGravityField." );
            }
            else
            {
                std::function< Eigen::VectorXd( const Eigen::MatrixXd&, const Eigen::MatrixXd& ) > accelerationFunction =
                        std::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::getAccelerationWithAlternativeCoefficients,
                                   sphericalHarmonicAcceleration, std::placeholders::_1, std::placeholders::_2 );

                std::shared_ptr< gravitation::GravityFieldVariations > gravityFieldVatiation =
                        timeDependentGravityField->getGravityFieldVariationsSet( )->getGravityFieldVariation(
                            accelerationVariableSettings->deformationType_,
                            accelerationVariableSettings->identifier_ ).second;

                std::function< Eigen::MatrixXd( ) > cosineCorrectionFunction =
                        std::bind( &gravitation::GravityFieldVariations::getLastCosineCorrection,
                                   gravityFieldVatiation  );
                std::function< Eigen::MatrixXd( ) > sineCorrectionFunction =
                        std::bind( &gravitation::GravityFieldVariations::getLastSineCorrection,
                                   gravityFieldVatiation );

                variableFunction = std::bind( evaluateBivariateReferenceFunction< Eigen::VectorXd, Eigen::MatrixXd >,
                                              accelerationFunction, cosineCorrectionFunction, sineCorrectionFunction );

            }
        }
        parameterSize = 3;

        break;
    }
    case single_gravity_field_variation_acceleration_terms:
    {
        std::shared_ptr< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings > accelerationVariableSettings =
                std::dynamic_pointer_cast< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >( dependentVariableSettings );
        if( accelerationVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type single_gravity_field_variation_acceleration_terms";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration =
                    std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                        getSphericalHarmonicAccelerationForDependentVariables(
                            accelerationVariableSettings, stateDerivativeModels, false ) );

            std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > timeDependentGravityField =
                    std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                        bodies.at( dependentVariableSettings->secondaryBody_ )->getGravityFieldModel( ) );

            if( timeDependentGravityField == nullptr )
            {
                throw std::runtime_error( "Error when requesting save of gravity field variation acceleration, central body " +
                                          dependentVariableSettings->secondaryBody_ +
                                          " has no TimeDependentSphericalHarmonicsGravityField." );
            }
            else
            {
                std::function< Eigen::VectorXd( const Eigen::MatrixXd&, const Eigen::MatrixXd& ) > accelerationFunction =
                        std::bind( &gravitation::SphericalHarmonicsGravitationalAccelerationModel::
                                   getAccelerationComponentsWithAlternativeCoefficients,
                                   sphericalHarmonicAcceleration, std::placeholders::_1, std::placeholders::_2, accelerationVariableSettings->componentIndices_ );

                std::shared_ptr< gravitation::GravityFieldVariations > gravityFieldVatiation =
                        timeDependentGravityField->getGravityFieldVariationsSet( )->getGravityFieldVariation(
                            accelerationVariableSettings->deformationType_,
                            accelerationVariableSettings->identifier_ ).second;

                std::function< Eigen::MatrixXd( ) > cosineCorrectionFunction =
                        std::bind( &gravitation::GravityFieldVariations::getLastCosineCorrection,
                                   gravityFieldVatiation );
                std::function< Eigen::MatrixXd( ) > sineCorrectionFunction =
                        std::bind( &gravitation::GravityFieldVariations::getLastSineCorrection,
                                   gravityFieldVatiation );

                variableFunction = std::bind( evaluateBivariateReferenceFunction< Eigen::VectorXd, Eigen::MatrixXd >,
                                              accelerationFunction, cosineCorrectionFunction, sineCorrectionFunction );

            }
            parameterSize = 3 * accelerationVariableSettings->componentIndices_.size( );

        }
        break;
    }
    case aerodynamic_force_coefficients_dependent_variable:
    {
        if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                    bodies.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
        {
            simulation_setup::addAtmosphericFlightConditions(
                        bodies, bodyWithProperty, secondaryBody );
        }

        variableFunction = std::bind(
                    &aerodynamics::AerodynamicCoefficientInterface::getCurrentForceCoefficients,
                    std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodies.at( bodyWithProperty )->getFlightConditions( ) )->getAerodynamicCoefficientInterface( ) );
        parameterSize = 3;

        break;
    }
    case aerodynamic_moment_coefficients_dependent_variable:
    {
        if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                    bodies.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
        {
            simulation_setup::addAtmosphericFlightConditions(
                        bodies, bodyWithProperty, secondaryBody );
        }

        variableFunction = std::bind(
                    &aerodynamics::AerodynamicCoefficientInterface::getCurrentMomentCoefficients,
                    std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodies.at( bodyWithProperty )->getFlightConditions( ) )->getAerodynamicCoefficientInterface( ) );
        parameterSize = 3;

        break;
    }
    case aerodynamic_control_surface_free_force_coefficients_dependent_variable:
    {
        if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                bodies.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
        {
            simulation_setup::addAtmosphericFlightConditions(
                    bodies, bodyWithProperty, secondaryBody );
        }

        variableFunction = std::bind(
                &aerodynamics::AerodynamicCoefficientInterface::getCurrentControlSurfaceFreeForceCoefficients,
                std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodies.at( bodyWithProperty )->getFlightConditions( ) )->getAerodynamicCoefficientInterface( ) );
        parameterSize = 3;

        break;
    }
    case aerodynamic_control_surface_free_moment_coefficients_dependent_variable:
    {
        if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                bodies.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
        {
            simulation_setup::addAtmosphericFlightConditions(
                    bodies, bodyWithProperty, secondaryBody );
        }

        variableFunction = std::bind(
                &aerodynamics::AerodynamicCoefficientInterface::getCurrentControlSurfaceFreeMomentCoefficients,
                std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodies.at( bodyWithProperty )->getFlightConditions( ) )->getAerodynamicCoefficientInterface( ) );
        parameterSize = 3;

        break;
    }
    case aerodynamic_control_surface_force_coefficients_increment_dependent_variable:
    {
        std::shared_ptr< ControlSurfaceCoefficientDependentVariableSettings > controlSurfaceVariabelSettings =
                std::dynamic_pointer_cast< ControlSurfaceCoefficientDependentVariableSettings >( dependentVariableSettings );
        if( controlSurfaceVariabelSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type aerodynamic_control_surface_force_coefficients_increment_dependent_variable";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                    bodies.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                simulation_setup::addAtmosphericFlightConditions(
                        bodies, bodyWithProperty, secondaryBody );
            }

            variableFunction = std::bind(
                    &aerodynamics::AerodynamicCoefficientInterface::getCurrentForceCoefficientIncrement,
                    std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                            bodies.at( bodyWithProperty )->getFlightConditions( ) )->getAerodynamicCoefficientInterface( ),
                            controlSurfaceVariabelSettings->controlSurfaceName_ );
            parameterSize = 3;
        }
        break;
    }
    case aerodynamic_control_surface_moment_coefficients_increment_dependent_variable:
    {
        std::shared_ptr< ControlSurfaceCoefficientDependentVariableSettings > controlSurfaceVariabelSettings =
                std::dynamic_pointer_cast< ControlSurfaceCoefficientDependentVariableSettings >( dependentVariableSettings );
        if( controlSurfaceVariabelSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type aerodynamic_control_surface_moment_coefficients_increment_dependent_variable";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                    bodies.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                simulation_setup::addAtmosphericFlightConditions(
                        bodies, bodyWithProperty, secondaryBody );
            }

            variableFunction = std::bind(
                    &aerodynamics::AerodynamicCoefficientInterface::getCurrentMomentCoefficientIncrement,
                    std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                            bodies.at( bodyWithProperty )->getFlightConditions( ) )->getAerodynamicCoefficientInterface( ),
                    controlSurfaceVariabelSettings->controlSurfaceName_ );
            parameterSize = 3;
        }
        break;
    }
    case inertial_to_body_fixed_rotation_matrix_variable:
    {
        std::function< Eigen::Quaterniond( ) > rotationFunction =
                std::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame, bodies.at( bodyWithProperty ) );
        variableFunction = std::bind( &getVectorRepresentationForRotationQuaternion, rotationFunction );
        parameterSize = 9;
        break;
    }
    case intermediate_aerodynamic_rotation_matrix_variable:
    {
        if( bodies.at( bodyWithProperty )->getFlightConditions( ) == nullptr )
        {
            simulation_setup::addFlightConditions(
                        bodies, bodyWithProperty, secondaryBody );
        }

        // Check input consistency.
        std::shared_ptr< IntermediateAerodynamicRotationVariableSaveSettings >
                intermediateAerodynamicRotationVariableSaveSettings =
                std::dynamic_pointer_cast< IntermediateAerodynamicRotationVariableSaveSettings >(
                    dependentVariableSettings );
        if( intermediateAerodynamicRotationVariableSaveSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type intermediate_aerodynamic_rotation_matrix_variable";
            throw std::runtime_error( errorMessage );
        }

        std::function< Eigen::Quaterniond( ) > rotationFunction =
                std::bind( &reference_frames::AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                           bodies.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                           intermediateAerodynamicRotationVariableSaveSettings->baseFrame_,
                           intermediateAerodynamicRotationVariableSaveSettings->targetFrame_ );

        variableFunction = std::bind( &getVectorRepresentationForRotationQuaternion, rotationFunction );
        parameterSize = 9;
        break;
    }
    case body_fixed_airspeed_based_velocity_variable:
    {
        if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                    bodies.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
        {
            simulation_setup::addAtmosphericFlightConditions(
                        bodies, bodyWithProperty, secondaryBody );
        }

        variableFunction = std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAirspeedBasedVelocity,
                                      std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                          bodies.at( bodyWithProperty )->getFlightConditions( ) ) );
        parameterSize = 3;
        break;
    }
    case body_fixed_groundspeed_based_velocity_variable:
    {
        if( bodies.at( bodyWithProperty )->getFlightConditions( ) == nullptr )
        {
            simulation_setup::addFlightConditions(
                        bodies, bodyWithProperty, secondaryBody );
        }

        if(  bodies.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ) == nullptr )
        {
            std::string errorMessage= "Error, no aerodynamic angle calculator when creating dependent variable function of type body_fixed_groundspeed_based_velocity_variable";
            throw std::runtime_error( errorMessage );
        }

        variableFunction = std::bind( &reference_frames::AerodynamicAngleCalculator::getCurrentGroundspeedBasedBodyFixedVelocity,
                                      bodies.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ) );
        parameterSize = 3;
        break;
    }
    case tnw_to_inertial_frame_rotation_dependent_variable:
    {
        std::function< Eigen::Vector6d( ) > vehicleStateFunction =
                std::bind( &simulation_setup::Body::getState, bodies.at( dependentVariableSettings->associatedBody_ ) );
        std::function< Eigen::Vector6d( ) > centralBodyStateFunction;

        if( ephemerides::isFrameInertial( dependentVariableSettings->secondaryBody_ ) )
        {
            centralBodyStateFunction =  [ ]( ){ return Eigen::Vector6d::Zero( ); };
        }
        else
        {
            centralBodyStateFunction =
                    std::bind( &simulation_setup::Body::getState, bodies.at( dependentVariableSettings->secondaryBody_ ) );
        }

        std::function< Eigen::Matrix3d( ) > rotationFunction =
                std::bind( &reference_frames::getTnwToInertialRotationFromFunctions,
                           vehicleStateFunction, centralBodyStateFunction, true );
        variableFunction = std::bind(
                    &getVectorRepresentationForRotationMatrixFunction, rotationFunction );

        parameterSize = 9;

        break;
    }
    case rsw_to_inertial_frame_rotation_dependent_variable:
    {
        std::function< Eigen::Vector6d( ) > vehicleStateFunction =
                std::bind( &simulation_setup::Body::getState, bodies.at( dependentVariableSettings->associatedBody_ ) );
        std::function< Eigen::Vector6d( ) > centralBodyStateFunction;

        if( ephemerides::isFrameInertial( dependentVariableSettings->secondaryBody_ ) )
        {
            centralBodyStateFunction =  [ ]( ){ return Eigen::Vector6d::Zero( ); };
        }
        else
        {
            centralBodyStateFunction =
                    std::bind( &simulation_setup::Body::getState, bodies.at( dependentVariableSettings->secondaryBody_ ) );
        }

        std::function< Eigen::Matrix3d( ) > rotationFunction =
                [=]( ){ return reference_frames::getRswSatelliteCenteredToInertialFrameRotationMatrix(
                        vehicleStateFunction( ) - centralBodyStateFunction( ) ); };
        variableFunction = std::bind(
                    &getVectorRepresentationForRotationMatrixFunction, rotationFunction );

        parameterSize = 9;

        break;
    }
    case total_torque_dependent_variable:
    {

        // Retrieve model responsible for computing accelerations of requested bodies.
        std::shared_ptr< RotationalMotionStateDerivative< StateScalarType, TimeType > > rotationalDynamicsModel =
                getRotationalStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
        variableFunction = std::bind( &RotationalMotionStateDerivative< StateScalarType, TimeType >::getTotalTorqueForBody,
                                      rotationalDynamicsModel, bodyWithProperty );
        parameterSize = 3;


        break;
    }
    case single_torque_dependent_variable:
    {
        // Check input consistency.
        std::shared_ptr< SingleTorqueDependentVariableSaveSettings > torqueDependentVariableSettings =
                std::dynamic_pointer_cast< SingleTorqueDependentVariableSaveSettings >( dependentVariableSettings );
        if( torqueDependentVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type single_torque_dependent_variable";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            // Retrieve list of suitable torque models (size should be one to avoid ambiguities)
            std::vector< std::shared_ptr< basic_astrodynamics::TorqueModel > >
                    listOfSuitableTorqueModels = getTorqueBetweenBodies(
                        torqueDependentVariableSettings->associatedBody_,
                        torqueDependentVariableSettings->secondaryBody_,
                        stateDerivativeModels, torqueDependentVariableSettings->torqueModelType_ );


            if( listOfSuitableTorqueModels.size( ) != 1 )
            {
                std::string errorMessage = "Error when getting torque between bodies " +
                        torqueDependentVariableSettings->associatedBody_ + " and " +
                        torqueDependentVariableSettings->secondaryBody_ + " of type " +
                        std::to_string(
                            torqueDependentVariableSettings->torqueModelType_ ) +
                        ", no such torque found";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                //std::function< Eigen::Vector3d( ) > vectorFunction =
                variableFunction = std::bind( &basic_astrodynamics::TorqueModel::getTorque,
                                              listOfSuitableTorqueModels.at( 0 ) );
                parameterSize = 3;
            }
        }
        break;
    }
    case keplerian_state_dependent_variable:
    {
        // Retrieve model responsible for computing accelerations of requested bodies.
        std::string orbitingBody = dependentVariableSettings->associatedBody_;
        std::string centralBody = dependentVariableSettings->secondaryBody_;

        if( bodies.count( orbitingBody ) == 0 )
        {
            throw std::runtime_error( "Error when requesting kepler elements as dependent variables, orbiting body: '" + orbitingBody + "' not found" );
        }
        else if( bodies.count( centralBody ) == 0 )
        {
            throw std::runtime_error( "Error when requesting kepler elements as dependent variables, central body: '" + centralBody + "' not found" );
        }

        std::function< double( ) > centralBodyGravitationalParameter;
        if( bodies.at( centralBody )->getGravityFieldModel( ) == nullptr )
        {
            throw std::runtime_error( "Error when requesting kepler elements as dependent variables, central body: '" + centralBody + "' has no gravity field" );
        }
        else
        {
            centralBodyGravitationalParameter = std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                           bodies.at( centralBody )->getGravityFieldModel( ) );
        }

        std::function< double( ) > orbitingBodyGravitationalParameter;
        std::function< double( ) > effectiveGravitationalParameter;
        if( bodies.at( orbitingBody )->getGravityFieldModel( ) != nullptr )
        {
            orbitingBodyGravitationalParameter = std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                            bodies.at( orbitingBody )->getGravityFieldModel( ) );
            effectiveGravitationalParameter = std::bind( &utilities::sumFunctionReturn< double >,
                                                         orbitingBodyGravitationalParameter,
                                                         centralBodyGravitationalParameter );

        }
        else
        {
            effectiveGravitationalParameter = centralBodyGravitationalParameter;
        }

        // Retrieve functions for positions of two bodies.
        std::function< Eigen::Vector6d( const Eigen::Vector6d&, const Eigen::Vector6d& ) > functionToEvaluate =
                std::bind( &linear_algebra::computeVectorDifference< 6 >, std::placeholders::_1, std::placeholders::_2 );
        std::function< Eigen::Vector6d( ) > firstInput =
                std::bind( &simulation_setup::Body::getState, bodies.at( orbitingBody ) );
        std::function< Eigen::Vector6d( ) > secondInput =
                std::bind( &simulation_setup::Body::getState, bodies.at( centralBody ) );
        std::function< Eigen::Vector6d( ) > relativeStateFunction = std::bind(
                    &evaluateBivariateReferenceFunction< Eigen::Vector6d, Eigen::Vector6d >,
                    functionToEvaluate, firstInput, secondInput );

        variableFunction = std::bind( &orbital_element_conversions::convertCartesianToKeplerianElementsFromFunctions< double >,
                                      relativeStateFunction, effectiveGravitationalParameter );
        parameterSize = 6;


        break;
    }
    case modified_equinocial_state_dependent_variable:
    {
        // Retrieve model responsible for computing accelerations of requested bodies.
        std::string orbitingBody = dependentVariableSettings->associatedBody_;
        std::string centralBody = dependentVariableSettings->secondaryBody_;

        if( bodies.count( orbitingBody ) == 0 )
        {
            throw std::runtime_error( "Error when requesting modified equinoctial elements as dependent variables, orbiting body: '" + orbitingBody + "' not found" );
        }
        else if( bodies.count( centralBody ) == 0 )
        {
            throw std::runtime_error( "Error when requesting modified equinoctial elements as dependent variables, central body: '" + centralBody + "' not found" );
        }

        std::function< double( ) > centralBodyGravitationalParameter;
        if( bodies.at( centralBody )->getGravityFieldModel( ) == nullptr )
        {
            throw std::runtime_error( "Error when requesting modified equinoctial elements as dependent variables, central body: '" + centralBody + "' has no gravity field" );
        }
        else
        {
            centralBodyGravitationalParameter = std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                           bodies.at( centralBody )->getGravityFieldModel( ) );
        }

        std::function< double( ) > orbitingBodyGravitationalParameter;
        std::function< double( ) > effectiveGravitationalParameter;
        if( bodies.at( orbitingBody )->getGravityFieldModel( ) != nullptr )
        {
            orbitingBodyGravitationalParameter = std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                            bodies.at( orbitingBody )->getGravityFieldModel( ) );
            effectiveGravitationalParameter = std::bind( &utilities::sumFunctionReturn< double >,
                                                         orbitingBodyGravitationalParameter,
                                                         centralBodyGravitationalParameter );

        }
        else
        {
            effectiveGravitationalParameter = centralBodyGravitationalParameter;
        }

        // Retrieve functions for positions of two bodies.
        std::function< Eigen::Vector6d( const Eigen::Vector6d&, const Eigen::Vector6d& ) > functionToEvaluate =
                std::bind( &linear_algebra::computeVectorDifference< 6 >, std::placeholders::_1, std::placeholders::_2 );
        std::function< Eigen::Vector6d( ) > firstInput =
                std::bind( &simulation_setup::Body::getState, bodies.at( orbitingBody ) );
        std::function< Eigen::Vector6d( ) > secondInput =
                std::bind( &simulation_setup::Body::getState, bodies.at( centralBody ) );
        std::function< Eigen::Vector6d( ) > relativeStateFunction = std::bind(
                    &evaluateBivariateReferenceFunction< Eigen::Vector6d, Eigen::Vector6d >,
                    functionToEvaluate, firstInput, secondInput );

        variableFunction = std::bind(
                    &orbital_element_conversions::convertCartesianToModifiedEquinoctialElementsFromStateFunction< double >,
                    relativeStateFunction, effectiveGravitationalParameter );
        parameterSize = 6;


        break;
    }
    case body_fixed_relative_cartesian_position:
    {
        std::function< Eigen::Vector3d( ) > positionFunctionOfRelativeBody =
                std::bind( &simulation_setup::Body::getPosition, bodies.at( bodyWithProperty ) );
        std::function< Eigen::Vector3d( ) > positionFunctionOfCentralBody =
                std::bind( &simulation_setup::Body::getPosition, bodies.at( secondaryBody ) );
        std::function< Eigen::Quaterniond( ) > orientationFunctionOfCentralBody =
                std::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame, bodies.at( secondaryBody ) );

        variableFunction = std::bind(
                    &reference_frames::getBodyFixedCartesianPosition, positionFunctionOfCentralBody,
                    positionFunctionOfRelativeBody, orientationFunctionOfCentralBody );
        parameterSize = 3;
        break;
    }
    case body_fixed_relative_spherical_position:
    {
        std::function< Eigen::Vector3d( ) > positionFunctionOfRelativeBody =
                std::bind( &simulation_setup::Body::getPosition, bodies.at( bodyWithProperty ) );
        std::function< Eigen::Vector3d( ) > positionFunctionOfCentralBody =
                std::bind( &simulation_setup::Body::getPosition, bodies.at( secondaryBody ) );
        std::function< Eigen::Quaterniond( ) > orientationFunctionOfCentralBody =
                std::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame, bodies.at( secondaryBody ) );

        variableFunction = std::bind(
                    &reference_frames::getBodyFixedSphericalPosition, positionFunctionOfCentralBody,
                    positionFunctionOfRelativeBody, orientationFunctionOfCentralBody );
        parameterSize = 3;
        break;
    }
    case euler_angles_to_body_fixed_313:
    {
        std::function< Eigen::Quaterniond( ) > orientationFunctionOfBody =
                std::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame, bodies.at( bodyWithProperty ) );

        std::function< Eigen::Vector3d( const Eigen::Quaterniond ) > eulerAngleFunction =
                std::bind( &basic_mathematics::get313EulerAnglesFromQuaternion, std::placeholders::_1 );

        variableFunction = std::bind(
                    &evaluateReferenceFunction< Eigen::Vector3d, Eigen::Quaterniond >,
                    eulerAngleFunction, orientationFunctionOfBody );
        parameterSize = 3;
        break;
    }
#if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
    case acceleration_partial_wrt_body_translational_state:
    {
        std::shared_ptr< AccelerationPartialWrtStateSaveSettings > accelerationPartialVariableSettings =
                std::dynamic_pointer_cast< AccelerationPartialWrtStateSaveSettings >( dependentVariableSettings );
        if( accelerationPartialVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type acceleration_partial_wrt_body_translational_state";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            if( stateDerivativePartials.count( translational_state ) == 0 )
            {
                throw std::runtime_error( "Error when requesting acceleration_partial_wrt_body_translational_state dependent variable, no translational state partials found." );
            }

            std::shared_ptr< acceleration_partials::AccelerationPartial > partialToUse =
                    getAccelerationPartialForBody(
                        stateDerivativePartials.at( translational_state ),
                        accelerationPartialVariableSettings->accelerationModelType_,
                        accelerationPartialVariableSettings->associatedBody_,
                        accelerationPartialVariableSettings->secondaryBody_ );

            std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int > partialFunction =
                    partialToUse->getDerivativeFunctionWrtStateOfIntegratedBody(
                        std::make_pair( accelerationPartialVariableSettings->derivativeWrtBody_, "" ),
                        propagators::translational_state );

            if( partialFunction.second == 0 )
            {
                variableFunction = [ = ]( ){ return Eigen::VectorXd::Zero( 18 ); };
            }
            else
            {
                variableFunction = std::bind( &getVectorFunctionFromBlockFunction, partialFunction.first, 3, 6 );
            }

            parameterSize = 18;
        }
        break;
    }
    case total_acceleration_partial_wrt_body_translational_state:
    {
        std::shared_ptr< TotalAccelerationPartialWrtStateSaveSettings > totalAccelerationPartialVariableSettings =
                std::dynamic_pointer_cast< TotalAccelerationPartialWrtStateSaveSettings >( dependentVariableSettings );
        if( totalAccelerationPartialVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type total_acceleration_partial_wrt_body_translational_state";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            if( stateDerivativePartials.count( translational_state ) == 0 )
            {
                throw std::runtime_error( "Error when requesting total_acceleration_partial_wrt_body_translational_state dependent variable, no translational state partials found." );
            }

            // Retrieve model responsible for computing accelerations of requested bodies.
            std::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > nBodyModel =
                    getTranslationalStateDerivativeModelForBody( totalAccelerationPartialVariableSettings->associatedBody_, stateDerivativeModels );

            // Retrieve acceleration models and create partials
            basic_astrodynamics::SingleBodyAccelerationMap accelerationModelList =
                    nBodyModel->getFullAccelerationsMap( ).at( totalAccelerationPartialVariableSettings->associatedBody_ );

            std::vector< std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int > > partialFunctions;

            // Parse all accelerations.
            for ( basic_astrodynamics::SingleBodyAccelerationMap::iterator itr = accelerationModelList.begin( ) ;
                  itr != accelerationModelList.end( ) ; itr++ )
            {
                std::string bodyExertingAcceleration = itr->first;

                for ( unsigned int currentAcceleration = 0 ; currentAcceleration < itr->second.size( ) ; currentAcceleration++ )
                {
                    std::shared_ptr< acceleration_partials::AccelerationPartial > partialToUse =
                            getAccelerationPartialForBody(
                                stateDerivativePartials.at( translational_state ),
                                basic_astrodynamics::getAccelerationModelType( itr->second[ currentAcceleration ] ),
                                totalAccelerationPartialVariableSettings->associatedBody_,
                                bodyExertingAcceleration );

                    std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int > partialFunction =
                            partialToUse->getDerivativeFunctionWrtStateOfIntegratedBody(
                                std::make_pair( totalAccelerationPartialVariableSettings->derivativeWrtBody_, "" ),
                                propagators::translational_state );

                    partialFunctions.push_back( partialFunction );
                }
            }

            variableFunction = [ = ]( )
            {
                Eigen::VectorXd variable = Eigen::VectorXd::Zero( 18 );

                for ( unsigned int i = 0 ; i < partialFunctions.size( ) ; i++ )
                {
                    if( partialFunctions[ i ].second != 0 )
                    {
                        variable += getVectorFunctionFromBlockFunction( partialFunctions[ i ].first, 3, 6 );
                    }
                }

                return variable;
            };

            parameterSize = 18;
        }
        break;
    }
#endif
    case minimum_constellation_distance:
    {
        std::shared_ptr< MinimumConstellationDistanceDependentVariableSaveSettings > minimumDistanceDependentVariable =
                std::dynamic_pointer_cast< MinimumConstellationDistanceDependentVariableSaveSettings >( dependentVariableSettings );
        if( minimumDistanceDependentVariable == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type minimum_constellation_distance";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::function< Eigen::Vector3d( ) > mainBodyPositionFunction =
                    std::bind( &simulation_setup::Body::getPosition, bodies.at( bodyWithProperty ) );
            std::vector< std::function< Eigen::Vector3d( ) > > bodiesToCheckPositionFunctions;
            for( unsigned int i = 0; i < minimumDistanceDependentVariable->bodiesToCheck_.size( ); i++ )
            {
                bodiesToCheckPositionFunctions.push_back(
                            std::bind( &simulation_setup::Body::getPosition, bodies.at( minimumDistanceDependentVariable->bodiesToCheck_.at( i ) ) ) );
            }

            variableFunction = std::bind( &getConstellationMinimumDistance, mainBodyPositionFunction, bodiesToCheckPositionFunctions );
            parameterSize = 2;
        }
        break;
    }
    case minimum_constellation_ground_station_distance:
    {
        std::shared_ptr< MinimumConstellationStationDistanceDependentVariableSaveSettings > minimumDistanceDependentVariable =
                std::dynamic_pointer_cast< MinimumConstellationStationDistanceDependentVariableSaveSettings >( dependentVariableSettings );
        if( minimumDistanceDependentVariable == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type minimum_constellation_ground_station_distance";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            std::function< Eigen::Vector3d( ) > stationPositionFunction =
                    std::bind( &simulation_setup::getGroundStationPositionDuringPropagation< double >,
                                   bodies.at( bodyWithProperty ), secondaryBody, bodies );
            std::shared_ptr< ground_stations::PointingAnglesCalculator > stationPointingAngleCalculator =
                    bodies.at( bodyWithProperty )->getGroundStation( secondaryBody )->getPointingAnglesCalculator( );

            std::vector< std::function< Eigen::Vector3d( ) > > bodiesToCheckPositionFunctions;
            for( unsigned int i = 0; i < minimumDistanceDependentVariable->bodiesToCheck_.size( ); i++ )
            {
                bodiesToCheckPositionFunctions.push_back(
                            std::bind( &simulation_setup::Body::getPosition,
                                       bodies.at( minimumDistanceDependentVariable->bodiesToCheck_.at( i ) ) ) );
            }

            variableFunction = [=]( ){ return getConstellationMinimumVisibleDistance(
                                          stationPositionFunction, bodiesToCheckPositionFunctions,
                                          stationPointingAngleCalculator, minimumDistanceDependentVariable->elevationAngleLimit_,
                            bodies.at( bodyWithProperty )->getDoubleTimeOfCurrentState( ) ); };
            parameterSize = 3;
        }
        break;
    }
    case body_center_of_mass:
    {
        if( bodies.at( bodyWithProperty )->getMassProperties( ) == nullptr )
        {
            throw std::runtime_error( "Error when saving body center of mass for " + bodyWithProperty + ", body has no mass properties." );
        }
        else
        {
            auto massProperties = bodies.at( bodyWithProperty )->getMassProperties( );
            variableFunction = [=]( ){ return massProperties->getCurrentCenterOfMass( ); };
            parameterSize = 3;
        }
        break;
    }
    case body_inertia_tensor:
    {
        if( bodies.at( bodyWithProperty )->getMassProperties( ) == nullptr )
        {
            throw std::runtime_error( "Error when saving body inertia tensor for " + bodyWithProperty + ", body has no mass properties." );
        }
        else
        {
            auto massProperties = bodies.at( bodyWithProperty )->getMassProperties( );
            variableFunction = [=]( ){ return getVectorRepresentationForRotationMatrix( massProperties->getCurrentInertiaTensor( ) ); };
            parameterSize = 9;
        }
        break;
    }
    case vehicle_panel_inertial_surface_normals:
    case vehicle_panel_body_fixed_surface_normals:
    {
        if( bodies.at( bodyWithProperty )->getVehicleSystems( ) == nullptr )
        {
            throw std::runtime_error( "Error when saving vehicle panel orientation of " + bodyWithProperty + ", body has no system models." );
        }
        else if( bodies.at( bodyWithProperty )->getVehicleSystems( )->getVehicleExteriorPanels( ).size( ) == 0 )
        {
            throw std::runtime_error( "Error when saving vehicle panel orientation of " + bodyWithProperty + ", body has no surface panels." );
        }
        else if( bodies.at( bodyWithProperty )->getVehicleSystems( )->getVehicleExteriorPanels( ).count( secondaryBody ) == 0 )
        {
            throw std::runtime_error( "Error when saving vehicle panel orientation of " + bodyWithProperty + ", body has no surface panels for part:" + secondaryBody + "." );
        }
        else
        {
            std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > partPanels =
                bodies.at( bodyWithProperty )->getVehicleSystems( )->getVehicleExteriorPanels( ).at( secondaryBody );
            parameterSize = 3 * partPanels.size( );

            variableFunction = [=]( )
            {

                Eigen::VectorXd partSurfaceNormals = Eigen::VectorXd::Zero( parameterSize );
                Eigen::Quaterniond normalRotation = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
                if( secondaryBody != "" )
                {
                    normalRotation = bodies.at( bodyWithProperty )->getVehicleSystems( )->getPartRotationToBaseFrame( secondaryBody );
                }
                if( dependentVariable == vehicle_panel_inertial_surface_normals )
                {
                    normalRotation = bodies.at( bodyWithProperty )->getCurrentRotationToGlobalFrame( ) * normalRotation;
                }
                for( unsigned int i = 0; i < partPanels.size( ); i++ )
                {
                    partSurfaceNormals.segment( i * 3, 3 ) = normalRotation * partPanels.at( i )->getFrameFixedSurfaceNormal( )( );
                }

                return partSurfaceNormals;
            };

        }
        break;
    }
    case vehicle_surface_panel_radiation_pressure_force:
    {
        auto radiationPressureAccelerationList = getAccelerationBetweenBodies(
            dependentVariableSettings->associatedBody_,
            dependentVariableSettings->secondaryBody_,
            stateDerivativeModels, basic_astrodynamics::radiation_pressure );

        if (radiationPressureAccelerationList.empty())
        {
            std::string errorMessage = "Error, radiation pressure acceleration with target " +
                                       dependentVariableSettings->associatedBody_ + " and source " +
                                       dependentVariableSettings->secondaryBody_  +
                                       " not found";
            throw std::runtime_error(errorMessage);
        }

        auto radiationPressureAcceleration =
            std::dynamic_pointer_cast<electromagnetism::RadiationPressureAcceleration>(
                radiationPressureAccelerationList.front());
        if( radiationPressureAcceleration == nullptr )
        {
            std::string errorMessage = "Error, radiation pressure acceleration with target " +
                                       dependentVariableSettings->associatedBody_ + " and source " +
                                       dependentVariableSettings->secondaryBody_  +
                                       " not found with compatible type";
            throw std::runtime_error(errorMessage);
        }
        else if( std::dynamic_pointer_cast< electromagnetism::PaneledRadiationPressureTargetModel >( radiationPressureAcceleration->getTargetModel( ) ) == nullptr )
        {
            std::string errorMessage = "Error, radiation pressure acceleration with target " +
                                       dependentVariableSettings->associatedBody_ + " and source " +
                                       dependentVariableSettings->secondaryBody_  +
                                       " has no panelled target model when saving per-panel radiation pressure force";
            throw std::runtime_error(errorMessage);
        }
        else
        {
            std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > paneledTarget =
                std::dynamic_pointer_cast< electromagnetism::PaneledRadiationPressureTargetModel >( radiationPressureAcceleration->getTargetModel( ) );
            parameterSize = 3 * paneledTarget->getTotalNumberOfPanels( );

            variableFunction = [=]( )
            {
                Eigen::VectorXd panelForces = Eigen::VectorXd::Zero( parameterSize );
                std::vector< Eigen::Vector3d >& panelForceList = paneledTarget->getPanelForces( secondaryBody );
                for( unsigned int i = 0; i < panelForceList.size( ); i++ )
                {
                    panelForces.segment( 3 * i, 3 ) = panelForceList.at( i );
                }
                return panelForces;
            };
        }
        break;
    }
    case paneled_radiation_source_per_panel_irradiance:
    case paneled_radiation_source_geometry:
    {
        auto radiationPressureAccelerationList = getAccelerationBetweenBodies(
            dependentVariableSettings->associatedBody_,
            dependentVariableSettings->secondaryBody_,
            stateDerivativeModels, basic_astrodynamics::radiation_pressure );
        if (radiationPressureAccelerationList.empty())
        {
            std::string errorMessage = "Error, radiation pressure acceleration with target " +
                                       dependentVariableSettings->associatedBody_ + " and source " +
                                       dependentVariableSettings->secondaryBody_  +
                                       " not found";
            throw std::runtime_error(errorMessage);
        }

        auto radiationPressureAcceleration =
            std::dynamic_pointer_cast<electromagnetism::RadiationPressureAcceleration>(
                radiationPressureAccelerationList.front());
        if( radiationPressureAcceleration == nullptr )
        {
            std::string errorMessage = "Error, radiation pressure acceleration with target " +
                                       dependentVariableSettings->associatedBody_ + " and source " +
                                       dependentVariableSettings->secondaryBody_  +
                                       " not found with compatible type";
            throw std::runtime_error(errorMessage);
        }
        else if( std::dynamic_pointer_cast< electromagnetism::PaneledSourceRadiationPressureAcceleration >( radiationPressureAcceleration ) == nullptr )
        {
            std::string errorMessage = "Error, radiation pressure acceleration with target " +
                                       dependentVariableSettings->associatedBody_ + " and source " +
                                       dependentVariableSettings->secondaryBody_  +
                                       " has no panelled source model when saving per-source-panel properties";
            throw std::runtime_error(errorMessage);
        }
        else
        {
            std::shared_ptr< electromagnetism::PaneledSourceRadiationPressureAcceleration > paneledRadiationPressureAcceleration =
                std::dynamic_pointer_cast< electromagnetism::PaneledSourceRadiationPressureAcceleration >( radiationPressureAcceleration );
            if( dependentVariable == paneled_radiation_source_per_panel_irradiance )
            {
                paneledRadiationPressureAcceleration->enableSavePanellingIrradiance( );
                parameterSize = paneledRadiationPressureAcceleration->getPaneledSourceModel( )->getNumberOfPanels( );
                variableFunction = [=]( )
                {
                    Eigen::VectorXd panelIrradiances = Eigen::VectorXd::Zero( parameterSize );
                    std::vector< double >& panelIrradianceList = paneledRadiationPressureAcceleration->getSavedPanelIrradiances( );
                    for( unsigned int i = 0; i < panelIrradianceList.size( ); i++ )
                    {
                        panelIrradiances( i ) = panelIrradianceList.at( i );
                    }

                    return panelIrradiances;
                };
            }
            else if( dependentVariable == paneled_radiation_source_geometry )
            {
                paneledRadiationPressureAcceleration->enableSavePanellingGeometry( );
                parameterSize = 7 * paneledRadiationPressureAcceleration->getPaneledSourceModel( )->getNumberOfPanels( );
                variableFunction = [=]( )
                {
                    Eigen::VectorXd panelGeometries = Eigen::VectorXd::Zero( parameterSize );
                    std::vector< Eigen::Vector7d >& panelGeometriesList = paneledRadiationPressureAcceleration->getSavedPanelGeometries( );
                    for( unsigned int i = 0; i < panelGeometriesList.size( ); i++ )
                    {
                        panelGeometries.segment( 7 * i, 7 ) = panelGeometriesList.at( i );
                    }

                    return panelGeometries;
                };
            }
        }
        break;
    }
    case nrlmsise_input_data:
    {
        std::string acceleratedBody = dependentVariableSettings->associatedBody_;
        std::string centralBody = dependentVariableSettings->secondaryBody_;


        if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
            bodies.at( acceleratedBody )->getFlightConditions( ) ) == nullptr )
        {
            simulation_setup::addAtmosphericFlightConditions(
                bodies, acceleratedBody, secondaryBody );
        }

        std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions =
            std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                           bodies.at( acceleratedBody )->getFlightConditions( ) );

        if ( std::dynamic_pointer_cast<aerodynamics::NRLMSISE00Atmosphere>(
            bodies.at( centralBody )->getAtmosphereModel( )) == nullptr )
        {
            throw std::runtime_error(
                "Error when saving nrlmsise input of " + centralBody + ", body has no NRLMSISE atmosphere." );
        }
        auto atmosphereModel =
            std::dynamic_pointer_cast<aerodynamics::NRLMSISE00Atmosphere>(
            bodies.at( centralBody )->getAtmosphereModel( ) );
        variableFunction = [=]( )
        {
            atmosphereModel->setInputStruct(
                flightConditions->getCurrentAltitude( ),
                flightConditions->getCurrentLongitude( ),
                flightConditions->getCurrentLatitude( ),
                flightConditions->getCurrentTime( ) );
            return aerodynamics::getNrlmsiseInputAsVector( atmosphereModel->getNRLMSISE00InputStruct() );
        };
        parameterSize = 17;
        break;
    }
    case custom_dependent_variable:
    {
        std::shared_ptr< CustomDependentVariableSaveSettings > customVariableSettings =
                std::dynamic_pointer_cast< CustomDependentVariableSaveSettings >( dependentVariableSettings );

        if( customVariableSettings == nullptr )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type custom_dependent_variable";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            variableFunction = [=]( )
            {
                Eigen::VectorXd customVariables = customVariableSettings->customDependentVariableFunction_( );
                if( customVariables.rows( ) != customVariableSettings->dependentVariableSize_ )
                {
                    throw std::runtime_error( "Error when retrieving custom dependent variable, actual size is different from pre-defined size" );
                }
                return customVariables;
            };
            parameterSize = customVariableSettings->dependentVariableSize_;
        }
        break;
    }
    default:
        std::string errorMessage =
                "Error, did not recognize vector dependent variable type when making variable function: " +
                std::to_string( dependentVariableSettings->dependentVariableType_ );
        throw std::runtime_error( errorMessage );
    }
    return std::make_pair( variableFunction, parameterSize );
}

//! Acces element at index function
/*!
 * Acces element at index function
 * \param vectorFunction The function returning an Eigen vector
 * \param index The index to be accessed
 * \return The value of vector at index
 */
inline double elementAtIndexFunction( const std::function< Eigen::VectorXd( ) >& vectorFunction, const int index )
{
    return vectorFunction( )( index );
}

//! Function to create a function returning a requested dependent variable value (of type double).
/*!
 *  Function to create a function returning a requested dependent variable value (of type double), retrieved from
 *  environment and/or state derivative models.
 *  \param dependentVariableSettings Settings for dependent variable that is to be returned by function created here.
 *  \param bodies List of bodies to use in simulations (containing full environment).
 *  \param stateDerivativeModels List of state derivative models used in simulations (sorted by dynamics type as key).
 *  \param stateDerivativePartials List of state derivative partials used in simulations (sorted by dynamics type as key).
 *  \return Function returning requested dependent variable. NOTE: The environment and state derivative models need to
 *  be updated to current state and independent variable before computation is performed.
 */
template< typename TimeType = double, typename StateScalarType = double >
std::function< double( ) > getDoubleDependentVariableFunction(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels =
        std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >( ),
        const std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >& stateDerivativePartials =
        std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >( ) )
{
    const int componentIndex = dependentVariableSettings->componentIndex_;
    const int dependentVariableSize = getDependentVariableSize( dependentVariableSettings, bodies );
    if ( componentIndex >= 0 )
    {
        if ( dependentVariableSettings->componentIndex_ > dependentVariableSize - 1 )
        {
            throw std::runtime_error( "Error, cannot access component of variable because it exceeds its size" );
        }

        const std::pair< std::function< Eigen::VectorXd( ) >, int > vectorFunction =
                getVectorDependentVariableFunction< TimeType, StateScalarType >(
                    dependentVariableSettings, bodies, stateDerivativeModels, stateDerivativePartials );
        return std::bind( &elementAtIndexFunction, vectorFunction.first, componentIndex );
    }
    else
    {
        std::function< double( ) > variableFunction;

        // Retrieve base information on dependent variable
        PropagationDependentVariables dependentVariable = dependentVariableSettings->dependentVariableType_;
        const std::string& bodyWithProperty = dependentVariableSettings->associatedBody_;
        const std::string& secondaryBody = dependentVariableSettings->secondaryBody_;

        // Check dependent variable type and create function accordingly.
        switch( dependentVariable )
        {
        case mach_number_dependent_variable:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodies.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                simulation_setup::addAtmosphericFlightConditions(
                            bodies, bodyWithProperty, secondaryBody );
            }

            std::function< double( const double, const double ) > functionToEvaluate =
                    std::bind( &aerodynamics::computeMachNumber, std::placeholders::_1, std::placeholders::_2 );

            // Retrieve functions for airspeed and speed of sound.
            std::function< double( ) > firstInput =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAirspeed,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodies.at( bodyWithProperty )->getFlightConditions( ) ) );
            std::function< double( ) > secondInput =
                    std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentSpeedOfSound,
                               std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                   bodies.at( bodyWithProperty )->getFlightConditions( ) ) );


            variableFunction = std::bind( &evaluateBivariateFunction< double, double >,
                                          functionToEvaluate, firstInput, secondInput );
            break;
        }
        case altitude_dependent_variable:
            if( bodies.at( bodyWithProperty )->getFlightConditions( ) == nullptr )
            {
                simulation_setup::addFlightConditions(
                            bodies, bodyWithProperty, secondaryBody );
            }
            variableFunction = std::bind( &aerodynamics::FlightConditions::getCurrentAltitude,
                                          bodies.at( bodyWithProperty )->getFlightConditions( ) );
            break;
        case airspeed_dependent_variable:
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodies.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                simulation_setup::addAtmosphericFlightConditions(
                            bodies, bodyWithProperty, secondaryBody );
            }
            variableFunction = std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAirspeed,
                                          std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                              bodies.at( bodyWithProperty )->getFlightConditions( ) ) );
            break;
        case local_density_dependent_variable:
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodies.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                simulation_setup::addAtmosphericFlightConditions(
                            bodies, bodyWithProperty, secondaryBody );
            }
            variableFunction = std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentDensity,
                                          std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                              bodies.at( bodyWithProperty )->getFlightConditions( ) ) );
            break;
        case radiation_pressure_dependent_variable:
        {
            auto radiationPressureAccelerationList = getAccelerationBetweenBodies(
                dependentVariableSettings->associatedBody_,
                dependentVariableSettings->secondaryBody_,
                stateDerivativeModels, basic_astrodynamics::radiation_pressure );

            if (radiationPressureAccelerationList.empty())
            {
                std::string errorMessage = "Error, radiation pressure acceleration with target " +
                        dependentVariableSettings->associatedBody_ + " and source " +
                        dependentVariableSettings->secondaryBody_  +
                       " not found";
                throw std::runtime_error(errorMessage);
            }

            auto radiationPressureAcceleration =
                    std::dynamic_pointer_cast<electromagnetism::RadiationPressureAcceleration>(
                            radiationPressureAccelerationList.front());

            variableFunction = [=] () { return radiationPressureAcceleration->getCurrentRadiationPressure(); };

            break;
        }
        case relative_distance_dependent_variable:
        {
            // Retrieve functions for positions of two bodies.
            std::function< double( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                    std::bind( &linear_algebra::computeNormOfVectorDifference, std::placeholders::_1, std::placeholders::_2 );
            std::function< Eigen::Vector3d( ) > firstInput =
                    std::bind( &simulation_setup::Body::getPosition, bodies.at( bodyWithProperty ) );

            std::function< Eigen::Vector3d( ) > secondInput;
            if( secondaryBody != "SSB" )
            {
                secondInput = std::bind( &simulation_setup::Body::getPosition, bodies.at( secondaryBody ) );
            }
            else if( simulation_setup::getGlobalFrameOrigin( bodies ) == "SSB" )
            {
                secondInput = []( ){ return Eigen::Vector3d::Zero( ); };
            }
            else
            {
                throw std::runtime_error( "Error, requested state of " + bodyWithProperty + " w.r.t. SSB, but SSB is not frame origin" );
            }

            variableFunction = std::bind(
                        &evaluateBivariateReferenceFunction< double, Eigen::Vector3d >, functionToEvaluate, firstInput, secondInput );
            break;
        }
        case relative_speed_dependent_variable:
        {
            // Retrieve functions for velicoty of two bodies.
            std::function< double( const Eigen::Vector3d&, const Eigen::Vector3d& ) > functionToEvaluate =
                    std::bind( &linear_algebra::computeNormOfVectorDifference, std::placeholders::_1, std::placeholders::_2 );
            std::function< Eigen::Vector3d( ) > firstInput =
                    std::bind( &simulation_setup::Body::getVelocity, bodies.at( bodyWithProperty ) );

            std::function< Eigen::Vector3d( ) > secondInput;
            if( secondaryBody != "SSB" )
            {
                secondInput = std::bind( &simulation_setup::Body::getVelocity, bodies.at( secondaryBody ) );
            }
            else if( simulation_setup::getGlobalFrameOrigin( bodies ) == "SSB" )
            {
                secondInput = []( ){ return Eigen::Vector3d::Zero( ); };
            }
            else
            {
                throw std::runtime_error( "Error, requested state of " + bodyWithProperty + " w.r.t. SSB, but SSB is not frame origin" );
            }


            variableFunction = std::bind(
                        &evaluateBivariateReferenceFunction< double, Eigen::Vector3d >, functionToEvaluate, firstInput, secondInput );

            break;
        }
        case single_acceleration_norm_dependent_variable:
        {
            // Check input consistency
            std::shared_ptr< SingleAccelerationDependentVariableSaveSettings > accelerationDependentVariableSettings =
                    std::dynamic_pointer_cast< SingleAccelerationDependentVariableSaveSettings >(
                        dependentVariableSettings );
            if( accelerationDependentVariableSettings == nullptr )
            {
                std::string errorMessage = "Error, inconsistent inout when creating dependent variable function of type single_acceleration_norm_dependent_variable";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                // Retrieve list of suitable acceleration models (size should be one to avoid ambiguities)
                std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                        listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                            accelerationDependentVariableSettings->associatedBody_,
                            accelerationDependentVariableSettings->secondaryBody_,
                            stateDerivativeModels, accelerationDependentVariableSettings->accelerationModelType_ );

                // Check if third-body counterpart of acceleration is found
                if( listOfSuitableAccelerationModels.size( ) == 0 && basic_astrodynamics::isAccelerationDirectGravitational(
                            accelerationDependentVariableSettings->accelerationModelType_ ) )
                {
                    listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                                accelerationDependentVariableSettings->associatedBody_,
                                accelerationDependentVariableSettings->secondaryBody_,
                                stateDerivativeModels, basic_astrodynamics::getAssociatedThirdBodyAcceleration(
                                    accelerationDependentVariableSettings->accelerationModelType_  ) );
                }

                if( listOfSuitableAccelerationModels.size( ) != 1 )
                {
                    std::string errorMessage = "Error when getting acceleration between bodies " +
                            accelerationDependentVariableSettings->associatedBody_ + " and " +
                            accelerationDependentVariableSettings->secondaryBody_ + " of type " +
                            std::to_string(
                                accelerationDependentVariableSettings->accelerationModelType_ ) +
                            ", no such acceleration found";
                    throw std::runtime_error( errorMessage );
                }
                else
                {
                    std::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > nBodyModel =
                            getTranslationalStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
                    std::map< std::string, std::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d > > removedAcceleration =
                            nBodyModel->getRemovedCentralAcceleration( );
                    if( removedAcceleration.count( accelerationDependentVariableSettings->associatedBody_ ) > 0 )
                    {
                        if( listOfSuitableAccelerationModels.at( 0 ) ==
                                removedAcceleration.at( accelerationDependentVariableSettings->associatedBody_ ) )
                        {
                            nBodyModel->setUpdateRemovedAcceleration( accelerationDependentVariableSettings->associatedBody_ );
                        }
                    }

                    std::function< Eigen::Vector3d( ) > vectorFunction =
                            std::bind( &basic_astrodynamics::AccelerationModel3d::getAcceleration,
                                       listOfSuitableAccelerationModels.at( 0 ) );
                    variableFunction = std::bind( &linear_algebra::getVectorNormFromFunction, vectorFunction );
                }
            }
            break;
        }
        case total_acceleration_norm_dependent_variable:
        {
            // Retrieve model responsible for computing accelerations of requested bodies.
            std::shared_ptr< NBodyStateDerivative< StateScalarType, TimeType > > nBodyModel =
                    getTranslationalStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
            nBodyModel->setUpdateRemovedAcceleration( dependentVariableSettings->associatedBody_ );
            std::function< Eigen::Vector3d( ) > vectorFunction =
                    std::bind( &NBodyStateDerivative< StateScalarType, TimeType >::getTotalAccelerationForBody,
                               nBodyModel, bodyWithProperty );
            variableFunction = std::bind( &linear_algebra::getVectorNormFromFunction, vectorFunction );

            break;
        }
        case total_torque_norm_dependent_variable:
        {

            // Retrieve model responsible for computing accelerations of requested bodies.
            std::shared_ptr< RotationalMotionStateDerivative< StateScalarType, TimeType > > rotationalDynamicsModel =
                    getRotationalStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
            std::function< Eigen::Vector3d( ) > vectorFunction  =
                    std::bind( &RotationalMotionStateDerivative< StateScalarType, TimeType >::getTotalTorqueForBody,
                               rotationalDynamicsModel, bodyWithProperty );
            variableFunction = std::bind( &linear_algebra::getVectorNormFromFunction, vectorFunction );

            break;
        }
        case single_torque_norm_dependent_variable:
        {
            // Check input consistency.
            std::shared_ptr< SingleTorqueDependentVariableSaveSettings > torqueDependentVariableSettings =
                    std::dynamic_pointer_cast< SingleTorqueDependentVariableSaveSettings >( dependentVariableSettings );
            if( torqueDependentVariableSettings == nullptr )
            {
                std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type single_torque_norm_dependent_variable";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                // Retrieve list of suitable torque models (size should be one to avoid ambiguities)
                std::vector< std::shared_ptr< basic_astrodynamics::TorqueModel > >
                        listOfSuitableTorqueModels = getTorqueBetweenBodies(
                            torqueDependentVariableSettings->associatedBody_,
                            torqueDependentVariableSettings->secondaryBody_,
                            stateDerivativeModels, torqueDependentVariableSettings->torqueModelType_ );


                if( listOfSuitableTorqueModels.size( ) != 1 )
                {
                    std::string errorMessage = "Error when getting torque between bodies " +
                            torqueDependentVariableSettings->associatedBody_ + " and " +
                            torqueDependentVariableSettings->secondaryBody_ + " of type " +
                            std::to_string(
                                torqueDependentVariableSettings->torqueModelType_ ) +
                            ", no such torque found";
                    throw std::runtime_error( errorMessage );
                }
                else
                {
                    //std::function< Eigen::Vector3d( ) > vectorFunction =
                    std::function< Eigen::Vector3d( ) > vectorFunction =
                            std::bind( &basic_astrodynamics::TorqueModel::getTorque,
                                       listOfSuitableTorqueModels.at( 0 ) );
                    variableFunction = std::bind( &linear_algebra::getVectorNormFromFunction, vectorFunction );
                }
            }
            break;
        }
        case relative_body_aerodynamic_orientation_angle_variable:
        {
            if( bodies.at( bodyWithProperty )->getFlightConditions( ) == nullptr )
            {
                simulation_setup::addFlightConditions(
                            bodies, bodyWithProperty, secondaryBody );
            }

            // Check input consistency.
            std::shared_ptr< BodyAerodynamicAngleVariableSaveSettings > bodyAerodynamicAngleVariableSaveSettings =
                    std::dynamic_pointer_cast< BodyAerodynamicAngleVariableSaveSettings >( dependentVariableSettings );
            if( bodyAerodynamicAngleVariableSaveSettings == nullptr )
            {
                std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type relative_body_aerodynamic_orientation_angle_variable";
                throw std::runtime_error( errorMessage );
            }

            variableFunction = std::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                                          bodies.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                                          bodyAerodynamicAngleVariableSaveSettings->angle_ );
            break;
        }
        case total_aerodynamic_g_load_variable:
        {
            // Check input consistency
            std::shared_ptr< SingleAccelerationDependentVariableSaveSettings > aerodynamicAccelerationSettings =
                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::aerodynamic, bodyWithProperty, secondaryBody, 0 );
            std::function< Eigen::VectorXd( ) > aerodynamicAccelerationFunction =
                    getVectorDependentVariableFunction( aerodynamicAccelerationSettings, bodies, stateDerivativeModels ).first;


            std::function< double( Eigen::Vector3d )> functionToEvaluate =
                    std::bind( &aerodynamics::computeAerodynamicLoadFromAcceleration, std::placeholders::_1 );
            variableFunction = std::bind(
                        &evaluateReferenceFunction< double, Eigen::Vector3d >,
                        functionToEvaluate, aerodynamicAccelerationFunction );

            break;
        }
        case stagnation_point_heat_flux_dependent_variable:
        {

            std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions =
                    std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodies.at( bodyWithProperty )->getFlightConditions( ) );
            if( flightConditions == nullptr )
            {
                simulation_setup::addAtmosphericFlightConditions(
                            bodies, bodyWithProperty, secondaryBody );
            }

            if( bodies.at( bodyWithProperty )->getVehicleSystems( ) == nullptr )
            {
                std::string errorMessage = "Error, no vehicle systems available when requesting stagnation point heating output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            std::shared_ptr< system_models::VehicleSystems > vehicleSystems =
                    bodies.at( bodyWithProperty )->getVehicleSystems( );


            if( !( vehicleSystems->getNoseRadius( ) == vehicleSystems->getNoseRadius( ) ) )
            {
                std::string errorMessage = "Error, no nose radius available when requesting stagnation point heating output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            if( !( vehicleSystems->getWallEmissivity( ) == vehicleSystems->getWallEmissivity( ) ) )
            {
                std::string errorMessage = "Error, no wall emmisivityavailable when requesting stagnation point heating output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            variableFunction = std::bind(
                        &computeEquilibriumFayRiddellHeatFluxFromProperties,
                        flightConditions, vehicleSystems );

            break;
        }
        case local_temperature_dependent_variable:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodies.at( bodyWithProperty )->getFlightConditions( ) )== nullptr )
            {
                simulation_setup::addAtmosphericFlightConditions(
                            bodies, bodyWithProperty, secondaryBody );
            }
            variableFunction = std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentFreestreamTemperature,
                                          std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                              bodies.at( bodyWithProperty )->getFlightConditions( ) ) );
            break;
        }
        case local_dynamic_pressure_dependent_variable:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodies.at( bodyWithProperty )->getFlightConditions( ) ) == nullptr )
            {
                simulation_setup::addAtmosphericFlightConditions(
                            bodies, bodyWithProperty, secondaryBody );
            }
            variableFunction = std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentDynamicPressure,
                                          std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                              bodies.at( bodyWithProperty )->getFlightConditions( ) ) );
            break;
        }
        case local_aerodynamic_heat_rate_dependent_variable:
        {
            if( std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodies.at( bodyWithProperty )->getFlightConditions( ) ) == nullptr )
            {
                simulation_setup::addAtmosphericFlightConditions(
                            bodies, bodyWithProperty, secondaryBody );
            }
            variableFunction = std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentAerodynamicHeatRate,
                                          std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                                              bodies.at( bodyWithProperty )->getFlightConditions( ) ) );
            break;
        }
        case geodetic_latitude_dependent_variable:
        {
            if( bodies.at( bodyWithProperty )->getFlightConditions( ) == nullptr )
            {
                simulation_setup::addFlightConditions(
                            bodies, bodyWithProperty, secondaryBody );
            }

            variableFunction = std::bind( &aerodynamics::FlightConditions::getCurrentGeodeticLatitude,
                                          bodies.at( bodyWithProperty )->getFlightConditions( ) );
            break;
        }
        case control_surface_deflection_dependent_variable:
        {
            if( bodies.at( bodyWithProperty )->getVehicleSystems( ) == nullptr )
            {
                std::string errorMessage = "Error, no vehicle systems available when requesting control surface deflection output of " +
                        bodyWithProperty + "with surface" + secondaryBody;
                throw std::runtime_error( errorMessage );
            }

            variableFunction = std::bind( &system_models::VehicleSystems::getCurrentControlSurfaceDeflection,
                                          bodies.at( bodyWithProperty )->getVehicleSystems( ),
                                          dependentVariableSettings->secondaryBody_ );
            break;
        }
        case total_mass_rate_dependent_variables:
        {
            // Retrieve model responsible for computing mass rate of requested bodies.
            std::shared_ptr< BodyMassStateDerivative< StateScalarType, TimeType > > nBodyModel =
                    getBodyMassStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
            variableFunction =
                    std::bind( &BodyMassStateDerivative< StateScalarType, TimeType >::getTotalMassRateForBody, nBodyModel,
                               bodyWithProperty );

            break;
        }
        case periapsis_altitude_dependent_variable:
        {
            using namespace Eigen;
            std::function< double( const Vector6d&, const double, const double ) > functionToEvaluate =
                    std::bind( &basic_astrodynamics::computePeriapsisAltitudeFromCartesianState, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );

            // Retrieve function for propagated body's Cartesian state in the global reference frame.
            std::function< Vector6d( ) > propagatedBodyStateFunction =
                    std::bind( &simulation_setup::Body::getState, bodies.at( bodyWithProperty ) );

            // Retrieve function for central body's Cartesian state in the global reference frame.
            std::function< Vector6d( ) > centralBodyStateFunction =
                    std::bind( &simulation_setup::Body::getState, bodies.at( secondaryBody ) );

            // Retrieve function for propagated body's Cartesian state in the propagation reference frame.
            std::function< Vector6d( ) > firstInput =
                    std::bind( &utilities::subtractFunctionReturn< Vector6d >,
                               propagatedBodyStateFunction, centralBodyStateFunction );

            // Retrieve function for central body's gravitational parameter.
            std::function< double( ) > secondInput =
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                               bodies.at( secondaryBody )->getGravityFieldModel( ) );

            // Retrieve function for central body's average radius.
            std::function< double( ) > thirdInput =
                    std::bind( &basic_astrodynamics::BodyShapeModel::getAverageRadius,
                               bodies.at( secondaryBody )->getShapeModel( ) );

            variableFunction = std::bind( &evaluateTrivariateFunction< double, Vector6d, double, double >,
                                          functionToEvaluate, firstInput, secondInput, thirdInput );
            break;
        }
        case apoapsis_altitude_dependent_variable:
        {
            using namespace Eigen;
            std::function< double( const Vector6d&, const double, const double ) > functionToEvaluate =
                    std::bind( &basic_astrodynamics::computeApoapsisAltitudeFromCartesianState, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );

            // Retrieve function for propagated body's Cartesian state in the global reference frame.
            std::function< Vector6d( ) > propagatedBodyStateFunction =
                    std::bind( &simulation_setup::Body::getState, bodies.at( bodyWithProperty ) );

            // Retrieve function for central body's Cartesian state in the global reference frame.
            std::function< Vector6d( ) > centralBodyStateFunction =
                    std::bind( &simulation_setup::Body::getState, bodies.at( secondaryBody ) );

            // Retrieve function for propagated body's Cartesian state in the propagation reference frame.
            std::function< Vector6d( ) > firstInput =
                    std::bind( &utilities::subtractFunctionReturn< Vector6d >,
                               propagatedBodyStateFunction, centralBodyStateFunction );

            // Retrieve function for central body's gravitational parameter.
            std::function< double( ) > secondInput =
                    std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                               bodies.at( secondaryBody )->getGravityFieldModel( ) );

            // Retrieve function for central body's average radius.
            std::function< double( ) > thirdInput =
                    std::bind( &basic_astrodynamics::BodyShapeModel::getAverageRadius,
                               bodies.at( secondaryBody )->getShapeModel( ) );

            variableFunction = std::bind( &evaluateTrivariateFunction< double, Vector6d, double, double >,
                                          functionToEvaluate, firstInput, secondInput, thirdInput );
            break;
        }
        case current_body_mass_dependent_variable:
        {
            variableFunction = std::bind(
                        &simulation_setup::Body::getBodyMass, bodies.at( bodyWithProperty ) );
            break;
        }
        case radiation_pressure_coefficient_dependent_variable:
        {
            if( bodies.at( bodyWithProperty )->getRadiationPressureTargetModel( ) == nullptr )
            {
                std::string errorMessage = "Error, no radiation pressure interface when requesting radiation pressure output of " +
                        bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }
            else if( std::dynamic_pointer_cast< electromagnetism::CannonballRadiationPressureTargetModel >(
                bodies.at( bodyWithProperty )-> getRadiationPressureTargetModel( ) ) == nullptr )
            {
                std::string errorMessage = "Error, no cannonball radiation pressure interface when requesting radiation pressure output of " +
                                           bodyWithProperty + "w.r.t." + secondaryBody;
                throw std::runtime_error( errorMessage );
            }
            variableFunction = std::bind( &electromagnetism::CannonballRadiationPressureTargetModel::getCoefficient,
                                          std::dynamic_pointer_cast< electromagnetism::CannonballRadiationPressureTargetModel >(
                                              bodies.at( bodyWithProperty )-> getRadiationPressureTargetModel( ) ) );
            break;
        }
        case gravity_field_potential_dependent_variable:
        {

            // Retrieve list of suitable acceleration models (size should be one to avoid ambiguities): either
            // point mass gravity, spherical harmonics gravity or polyhedron gravity
            basic_astrodynamics::AvailableAcceleration selectedAccelerationModelType;
            std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > selectedAccelerationModel =
                    getGravitationalAccelerationForDependentVariables( dependentVariableSettings, stateDerivativeModels,
                                                                       selectedAccelerationModelType );

            // Position functions
            std::function< Eigen::Vector3d( ) > positionFunctionOfRelativeBody =
                    std::bind( &simulation_setup::Body::getPosition, bodies.at( bodyWithProperty ) );
            std::function< Eigen::Vector3d( ) > positionFunctionOfCentralBody =
                    std::bind( &simulation_setup::Body::getPosition, bodies.at( secondaryBody ) );

            // Retrieve orientation function depending on type of gravity field
            std::function< Eigen::Quaterniond( ) > orientationFunctionOfCentralBody;

            // If gravity field is point mass, don't need orientation function
            if ( selectedAccelerationModelType == basic_astrodynamics::point_mass_gravity ||
                selectedAccelerationModelType == basic_astrodynamics::third_body_point_mass_gravity )
            {
                orientationFunctionOfCentralBody = [=]( ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); };
            }
            // Else get orientation function
            else
            {
                orientationFunctionOfCentralBody = std::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame, bodies.at( secondaryBody ) );
            }

            // Retrieve function to get body fixed position of body
            std::function< Eigen::Vector3d() > bodyFixedPositionOfBodyWithProperty =
                    std::bind( &reference_frames::getBodyFixedCartesianPosition, positionFunctionOfCentralBody,
                               positionFunctionOfRelativeBody, orientationFunctionOfCentralBody );

            // Downcast gravity model and get variable function
            if ( selectedAccelerationModelType == basic_astrodynamics::spherical_harmonic_gravity )
            {
                std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel >
                        sphericalHarmonicsAccelerationModel =  std::dynamic_pointer_cast<
                                gravitation::SphericalHarmonicsGravitationalAccelerationModel >( selectedAccelerationModel );

                sphericalHarmonicsAccelerationModel->resetUpdatePotential( true );
                variableFunction = [=]( ){ return sphericalHarmonicsAccelerationModel->getCurrentPotential( ); };
            }
            else if ( selectedAccelerationModelType == basic_astrodynamics::third_body_spherical_harmonic_gravity )
            {
                std::shared_ptr< gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel >
                        thirdBodySphericalHarmonicsAccelerationModel =  std::dynamic_pointer_cast<
                                gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel >( selectedAccelerationModel );

                std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel >
                        sphericalHarmonicsAccelerationModel =  std::dynamic_pointer_cast<
                                gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                                        thirdBodySphericalHarmonicsAccelerationModel->getAccelerationModelForBodyUndergoingAcceleration( ) );

                sphericalHarmonicsAccelerationModel->resetUpdatePotential( true );
                variableFunction = [=]( ){ return sphericalHarmonicsAccelerationModel->getCurrentPotential( ); };
            }
            else if ( selectedAccelerationModelType == basic_astrodynamics::point_mass_gravity )
            {
                std::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d >
                        pointMassAccelerationModel =  std::dynamic_pointer_cast<
                                gravitation::CentralGravitationalAccelerationModel3d >( selectedAccelerationModel );
                pointMassAccelerationModel->resetUpdatePotential( true );
                variableFunction = [=]( ){ return pointMassAccelerationModel->getCurrentPotential( ); };
            }
            else if ( selectedAccelerationModelType == basic_astrodynamics::third_body_point_mass_gravity )
            {

                std::shared_ptr< gravitation::ThirdBodyCentralGravityAcceleration >
                        thirdBodyPointMassAccelerationModel =  std::dynamic_pointer_cast<
                                gravitation::ThirdBodyCentralGravityAcceleration >( selectedAccelerationModel );

                std::shared_ptr< gravitation::CentralGravitationalAccelerationModel3d >
                        pointMassAccelerationModel =  std::dynamic_pointer_cast<
                                gravitation::CentralGravitationalAccelerationModel3d >(
                                        thirdBodyPointMassAccelerationModel->getAccelerationModelForBodyUndergoingAcceleration( ) );

                pointMassAccelerationModel->resetUpdatePotential( true );
                variableFunction = [=]( ){ return pointMassAccelerationModel->getCurrentPotential( ); };
            }
            else if ( selectedAccelerationModelType == basic_astrodynamics::polyhedron_gravity )
            {
                std::shared_ptr< gravitation::PolyhedronGravitationalAccelerationModel >
                        polyhedronAccelerationModel =  std::dynamic_pointer_cast<
                                gravitation::PolyhedronGravitationalAccelerationModel >( selectedAccelerationModel );
                polyhedronAccelerationModel->resetUpdatePotential( true );
                variableFunction = [=]( ){ return polyhedronAccelerationModel->getCurrentPotential( ); };
            }
            else if ( selectedAccelerationModelType == basic_astrodynamics::third_body_polyhedron_gravity )
            {
                std::shared_ptr< gravitation::ThirdBodyPolyhedronGravitationalAccelerationModel >
                        thirdBodyPolyhedronAccelerationModel =  std::dynamic_pointer_cast<
                                gravitation::ThirdBodyPolyhedronGravitationalAccelerationModel >( selectedAccelerationModel );

                std::shared_ptr< gravitation::PolyhedronGravitationalAccelerationModel >
                        polyhedronAccelerationModel =  std::dynamic_pointer_cast<
                                gravitation::PolyhedronGravitationalAccelerationModel >(
                                        thirdBodyPolyhedronAccelerationModel->getAccelerationModelForBodyUndergoingAcceleration( ) );

                polyhedronAccelerationModel->resetUpdatePotential( true );
                variableFunction = [=]( ){ return polyhedronAccelerationModel->getCurrentPotential( ); };
            }
            else if ( selectedAccelerationModelType == basic_astrodynamics::ring_gravity )
            {
                std::shared_ptr< gravitation::RingGravitationalAccelerationModel >
                        ringAccelerationModel =  std::dynamic_pointer_cast<
                                gravitation::RingGravitationalAccelerationModel >( selectedAccelerationModel );
                ringAccelerationModel->resetUpdatePotential( true );
                variableFunction = [=]( ){ return ringAccelerationModel->getCurrentPotential( ); };
            }
            else if ( selectedAccelerationModelType == basic_astrodynamics::third_body_ring_gravity )
            {
                std::shared_ptr< gravitation::ThirdBodyRingGravitationalAccelerationModel >
                        thirdBodyRingAccelerationModel =  std::dynamic_pointer_cast<
                                gravitation::ThirdBodyRingGravitationalAccelerationModel >( selectedAccelerationModel );

                std::shared_ptr< gravitation::RingGravitationalAccelerationModel >
                        ringAccelerationModel =  std::dynamic_pointer_cast<
                                gravitation::RingGravitationalAccelerationModel >(
                                        thirdBodyRingAccelerationModel->getAccelerationModelForBodyUndergoingAcceleration( ) );

                ringAccelerationModel->resetUpdatePotential( true );
                variableFunction = [=]( ){ return ringAccelerationModel->getCurrentPotential( ); };
            }
            else
            {
                std::string errorMessage = "Error, when setting up gravitational potential as dependent variable, for" +
                        std::to_string( selectedAccelerationModelType ) + " gravity model of " + bodyWithProperty +
                        "w.r.t." + secondaryBody + "acceleration model does not have the computation of potential implemented.";
                throw std::runtime_error( errorMessage );
            }

            break;
        }
        case gravity_field_laplacian_of_potential_dependent_variable:
        {
            // Retrieve list of suitable acceleration models (size should be one to avoid ambiguities): either
            // point mass gravity, spherical harmonics gravity or polyhedron gravity
            basic_astrodynamics::AvailableAcceleration selectedAccelerationModelType;
            std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > selectedAccelerationModel =
                    getGravitationalAccelerationForDependentVariables( dependentVariableSettings, stateDerivativeModels,
                                                                       selectedAccelerationModelType );

            // Position functions
            std::function< Eigen::Vector3d ( ) > positionFunctionOfRelativeBody =
                    std::bind( &simulation_setup::Body::getPosition, bodies.at( bodyWithProperty ) );
            std::function< Eigen::Vector3d ( ) > positionFunctionOfCentralBody =
                    std::bind( &simulation_setup::Body::getPosition, bodies.at( secondaryBody ) );

            // Retrieve orientation function depending on type of gravity field
            std::function< Eigen::Quaterniond ( ) > orientationFunctionOfCentralBody;

            // If gravity field is point mass, don't need orientation function
            if ( selectedAccelerationModelType == basic_astrodynamics::point_mass_gravity ||
                 selectedAccelerationModelType == basic_astrodynamics::third_body_point_mass_gravity ) {
                orientationFunctionOfCentralBody = [ = ] ( ) {
                    return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
                };
            }
                // Else get orientation function
            else {
                orientationFunctionOfCentralBody = std::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame,
                                                              bodies.at( secondaryBody ) );
            }

            // Retrieve function to get body fixed position of body
            std::function< Eigen::Vector3d ( ) > bodyFixedPositionOfBodyWithProperty =
                    std::bind( &reference_frames::getBodyFixedCartesianPosition, positionFunctionOfCentralBody,
                               positionFunctionOfRelativeBody, orientationFunctionOfCentralBody );

            // Downcast gravity model and get variable function
            if ( selectedAccelerationModelType == basic_astrodynamics::polyhedron_gravity ) {
                std::shared_ptr< gravitation::PolyhedronGravitationalAccelerationModel >
                        polyhedronAccelerationModel = std::dynamic_pointer_cast<
                        gravitation::PolyhedronGravitationalAccelerationModel >( selectedAccelerationModel );
                polyhedronAccelerationModel->resetUpdateLaplacianOfPotential( true );
                variableFunction = [ = ] ( ) { return polyhedronAccelerationModel->getCurrentLaplacianOfPotential( ); };
            }
            else if ( selectedAccelerationModelType == basic_astrodynamics::third_body_polyhedron_gravity ) {
                std::shared_ptr< gravitation::ThirdBodyPolyhedronGravitationalAccelerationModel >
                        thirdBodyPolyhedronAccelerationModel = std::dynamic_pointer_cast<
                        gravitation::ThirdBodyPolyhedronGravitationalAccelerationModel >( selectedAccelerationModel );

                std::shared_ptr< gravitation::PolyhedronGravitationalAccelerationModel >
                        polyhedronAccelerationModel = std::dynamic_pointer_cast<
                        gravitation::PolyhedronGravitationalAccelerationModel >(
                        thirdBodyPolyhedronAccelerationModel->getAccelerationModelForBodyUndergoingAcceleration( ) );

                polyhedronAccelerationModel->resetUpdateLaplacianOfPotential( true );
                variableFunction = [ = ] ( ) { return polyhedronAccelerationModel->getCurrentLaplacianOfPotential( ); };
            }
            else {
                std::string errorMessage = "Error, when setting up gravitational potential as dependent variable, for " +
                                           std::to_string( selectedAccelerationModelType ) + " gravity model of " +
                                           bodyWithProperty +
                                           " w.r.t. " + secondaryBody +
                                           ": acceleration model does not have the computation of potential implemented.";
                throw std::runtime_error( errorMessage );
            }

            break;
        }
        case received_irradiance:
        {
            auto radiationPressureAccelerationList = getAccelerationBetweenBodies(
                dependentVariableSettings->associatedBody_,
                dependentVariableSettings->secondaryBody_,
                stateDerivativeModels, basic_astrodynamics::radiation_pressure );

            if (radiationPressureAccelerationList.empty())
            {
                std::string errorMessage = "Error, radiation pressure acceleration with target " +
                        dependentVariableSettings->associatedBody_ + " and source " +
                        dependentVariableSettings->secondaryBody_  +
                       " not found";
                throw std::runtime_error(errorMessage);
            }

            auto radiationPressureAcceleration =
                    std::dynamic_pointer_cast<electromagnetism::RadiationPressureAcceleration>(
                            radiationPressureAccelerationList.front());

            variableFunction = [=] () { return radiationPressureAcceleration->getReceivedIrradiance(); };

            break;
        }
        case received_fraction:
        {
            auto radiationPressureAccelerationList = getAccelerationBetweenBodies(
                dependentVariableSettings->associatedBody_,
                dependentVariableSettings->secondaryBody_,
                stateDerivativeModels, basic_astrodynamics::radiation_pressure );

            if (radiationPressureAccelerationList.empty())
            {
                std::string errorMessage = "Error, radiation pressure acceleration with target " +
                        dependentVariableSettings->associatedBody_ + " and source " +
                        dependentVariableSettings->secondaryBody_  +
                       " not found";
                throw std::runtime_error(errorMessage);
            }

            auto radiationPressureAcceleration =
                    std::dynamic_pointer_cast<electromagnetism::IsotropicPointSourceRadiationPressureAcceleration>(
                            radiationPressureAccelerationList.front());
            if (radiationPressureAcceleration == nullptr)
            {
                throw std::runtime_error("Error, source body " + dependentVariableSettings->secondaryBody_ +
                    " does not have an isotropic point source, which is required for the received fraction dependent " +
                    "variable");
            }

            variableFunction = [=] () { return radiationPressureAcceleration->getSourceToTargetReceivedFraction(); };

            break;
        }
        case visible_and_emitting_source_panel_count:
        case visible_source_area:
        {
            auto radiationPressureAccelerationList = getAccelerationBetweenBodies(
                dependentVariableSettings->associatedBody_,
                dependentVariableSettings->secondaryBody_,
                stateDerivativeModels, basic_astrodynamics::radiation_pressure );

            if (radiationPressureAccelerationList.empty())
            {
                std::string errorMessage = "Error, radiation pressure acceleration with target " +
                        dependentVariableSettings->associatedBody_ + " and source " +
                        dependentVariableSettings->secondaryBody_  +
                       " not found";
                throw std::runtime_error(errorMessage);
            }

            auto radiationPressureAcceleration =
                    std::dynamic_pointer_cast<electromagnetism::PaneledSourceRadiationPressureAcceleration>(
                            radiationPressureAccelerationList.front());
            if (radiationPressureAcceleration == nullptr)
            {
                if (dependentVariable == visible_source_area) {
                    throw std::runtime_error("Error, source body " + dependentVariableSettings->secondaryBody_ +
                         " does not have a paneled source, which is required for the visible source area "+
                         "dependent variable");
                }
                else
                {
                    throw std::runtime_error("Error, source body " + dependentVariableSettings->secondaryBody_ +
                         " does not have a paneled source, which is required for the visible/illuminated source "+
                         "panel count dependent variable");
                }
            }

            switch (dependentVariable)
            {
                case visible_and_emitting_source_panel_count:
                    variableFunction = [=] () { return radiationPressureAcceleration->getVisibleAndEmittingSourcePanelCount(); };
                    break;
                case visible_source_area:
                    variableFunction = [=] () { return radiationPressureAcceleration->getVisibleSourceArea(); };
                    break;
                default:
                    break;
            }

            break;
        }
        case custom_dependent_variable:
        {
            std::shared_ptr< CustomDependentVariableSaveSettings > customVariableSettings =
                    std::dynamic_pointer_cast< CustomDependentVariableSaveSettings >( dependentVariableSettings );

            if( customVariableSettings == nullptr )
            {
                std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type custom_dependent_variable";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                variableFunction = [=]( )
                {
                    Eigen::VectorXd customVariables = customVariableSettings->customDependentVariableFunction_( );
                    if( customVariables.rows( ) != 1 )
                    {
                        throw std::runtime_error( "Error when retrieving size-1 custom dependent variable, actual size is " + std::to_string( customVariables.rows( ) ) );
                    }
                    return customVariables( 0 );
                };
            }
            break;
        }
        default:
            std::string errorMessage =
                    "Error, did not recognize double dependent variable type when making variable function: " +
                    std::to_string( dependentVariableSettings->dependentVariableType_ );
            throw std::runtime_error( errorMessage );
        }

        return variableFunction;
    }
}


//! Function to return a vector containing only one value given by doubleFunction
/*!
 *
 * \param doubleFunction Function returning the double value
 * \return The vector containing only one element
 */
Eigen::VectorXd getVectorFromDoubleFunction( const std::function< double( ) >& doubleFunction );

//! Function to evaluate a set of vector-returning functions and concatenate the results.
/*!
 * Function to evaluate a set of vector-returning functions and concatenate the results.
 * \param vectorFunctionList List of functions returning vector variables (pairs denote function and return vector size)
 * \param totalSize Total size of concatenated vector (used as input for efficiency.
 * \return Concatenated results from input functions.
 */
Eigen::VectorXd evaluateListOfVectorFunctions(
        const std::vector< std::pair< std::function< Eigen::VectorXd( ) >, int > > vectorFunctionList,
        const int totalSize );

//! Function to create a function that evaluates a list of dependent variables and concatenates the results.
/*!
 *  Function to create a function that evaluates a list of dependent variables and concatenates the results.
 *  Dependent variables functions are created inside this function from a list of settings on their required
 *  types/properties.
 *  \param saveSettings Object containing types and other properties of dependent variables.
 *  \param bodies List of bodies to use in simulations (containing full environment).
 *  \param stateDerivativeModels List of state derivative models used in simulations (sorted by dynamics type as key)
 *  \return Pair with function returning requested dependent variable values, and list variable names with start entries.
 *  NOTE: The environment and state derivative models need to
 *  be updated to current state and independent variable before computation is performed.
 */
template< typename TimeType = double, typename StateScalarType = double >
std::pair< std::function< Eigen::VectorXd( ) >, std::map< std::pair< int, int >, std::string > > createDependentVariableListFunction(
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables,
        const simulation_setup::SystemOfBodies& bodies,
        std::map< std::pair< int, int >, std::shared_ptr< SingleDependentVariableSaveSettings > >& orderedDependentVariables,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > >& stateDerivativeModels,
        const std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >& stateDerivativePartials =
        std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >( ) )
{
    // create list of vector parameters
    std::vector< std::pair< std::function< Eigen::VectorXd( ) >, int > > vectorFunctionList;
    std::vector< std::pair< std::string, int > > vectorVariableList;

    for( std::shared_ptr< SingleDependentVariableSaveSettings > variable: dependentVariables )
    {
        std::pair< std::function< Eigen::VectorXd( ) >, int > vectorFunction;
        // Create double parameter
        if( isScalarDependentVariable( variable, bodies ) )
        {
#if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
            std::function< double( ) > doubleFunction =
                    getDoubleDependentVariableFunction( variable, bodies, stateDerivativeModels, stateDerivativePartials );
#else
            std::function< double( ) > doubleFunction =
                    getDoubleDependentVariableFunction( variable, bodies, stateDerivativeModels );
#endif
            vectorFunction = std::make_pair( std::bind( &getVectorFromDoubleFunction, doubleFunction ), 1 );
        }
        // Create vector parameter
        else
        {
#if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
            vectorFunction = getVectorDependentVariableFunction(
                        variable, bodies, stateDerivativeModels, stateDerivativePartials );
#else
            vectorFunction =
                    getVectorDependentVariableFunction( variable, bodies, stateDerivativeModels );
#endif
        }
        vectorFunctionList.push_back( vectorFunction );
        vectorVariableList.push_back( std::make_pair( getDependentVariableId( variable ), vectorFunction.second ) );
    }

    // Set list of variable ids/indices in correc otder.
    int totalVariableSize = 0;
    std::map< std::pair< int, int >, std::string > dependentVariableIds;

    int variableCounter = 0;
    for( std::pair< std::string, int > vectorVariable: vectorVariableList )
    {
        dependentVariableIds[ { totalVariableSize, vectorVariable.second } ] = vectorVariable.first;
        orderedDependentVariables[ { totalVariableSize, vectorVariable.second } ] = dependentVariables.at( variableCounter );
        variableCounter++;
        totalVariableSize += vectorVariable.second;
    }

    // Create function conatenating function results.
    return std::make_pair( std::bind( &evaluateListOfVectorFunctions, vectorFunctionList, totalVariableSize ),
                           dependentVariableIds );
}

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATIONOUTPUT_H
