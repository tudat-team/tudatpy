/*    Copyright (c) 2010-2017, Delft University of Technology
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
#include "Tudat/Astrodynamics/Ephemerides/frameManager.h"
#include "Tudat/Astrodynamics/Propagators/dynamicsStateDerivativeModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationOutputSettings.h"

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
        const boost::function< OutputType( const InputType&, const InputType& ) > functionToEvaluate,
        const boost::function< InputType( ) > firstInput,
        const boost::function< InputType( ) > secondInput )
{
    return functionToEvaluate( firstInput( ), secondInput( ) );
}

template< typename OutputType, typename InputType >
OutputType evaluateReferenceFunction(
        const boost::function< OutputType( const InputType& ) > functionToEvaluate,
        const boost::function< InputType( ) > firstInput )
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
        const boost::function< Eigen::Matrix3d( ) > rotationFunction );

//! Get the vector representation of a quaternion.
/*!
 *  Get the vector representation of a quaternion. Quaternion is converted to a rotation matrix, which is then put into
 *  a vector representation.
 *  \param rotationFunction Function returning the quaternion that is to be put inot vector rerpesentation
 *  \return Column vector consisting of transpose of concatenated rows of matrix representation of rotationFunction input.
 */
Eigen::VectorXd getVectorRepresentationForRotationQuaternion(
        const boost::function< Eigen::Quaterniond( ) > rotationFunction );


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

//! Function to compute the Fay-Riddell equilibrium heat flux from body properties
/*!
 * Function to compute the Fay-Riddell equilibrium heat flux from body properties
 * \param flightConditions Object describing the current atmospheric flight conditions of the vehicle
 * \param vehicleSystems Object describing the physical properties of the vehicle
 * \return Equilibrium heat flux according to Fay-Riddell model
 */
double computeEquilibriumFayRiddellHeatFluxFromProperties(
        const boost::shared_ptr< aerodynamics::FlightConditions > flightConditions,
        const boost::shared_ptr< system_models::VehicleSystems > vehicleSystems );


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
            std::string errorMessage = "Error, no flight conditions available when requesting Mach number output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
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


        variableFunction = boost::bind( &evaluateBivariateFunction< double, double >,
                                        functionToEvaluate, firstInput, secondInput );
        break;
    }
    case altitude_dependent_variable:
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {
            std::string errorMessage = "Error, no flight conditions available when requesting altitude output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }
        variableFunction = boost::bind( &aerodynamics::FlightConditions::getCurrentAltitude,
                                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
        break;
    case airspeed_dependent_variable:
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {
            std::string errorMessage = "Error, no flight conditions available when requesting airspeed output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }
        variableFunction = boost::bind( &aerodynamics::FlightConditions::getCurrentAirspeed,
                                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
        break;
    case local_density_dependent_variable:
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {
            std::string errorMessage = "Error, no flight conditions available when requesting density output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }
        variableFunction = boost::bind( &aerodynamics::FlightConditions::getCurrentDensity,
                                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
        break;
    case radiation_pressure_dependent_variable:
        if( bodyMap.at( bodyWithProperty )->getRadiationPressureInterfaces( ).count( secondaryBody ) == 0 )
        {
            std::string errorMessage = "Error, no radiation pressure interfaces when requesting radiation pressure output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
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
                    &evaluateBivariateReferenceFunction< double, Eigen::Vector3d >, functionToEvaluate, firstInput, secondInput );
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
                    &evaluateBivariateReferenceFunction< double, Eigen::Vector3d >, functionToEvaluate, firstInput, secondInput );

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

            // Check if thirfd-body counterpart of acceleration is found
            if( listOfSuitableAccelerationModels.size( ) == 0 && basic_astrodynamics::isAccelerationDirectGravitational(
                        accelerationDependentVariableSettings->accelerationModeType_ ) )
            {
                listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                    accelerationDependentVariableSettings->associatedBody_,
                    accelerationDependentVariableSettings->secondaryBody_,
                    stateDerivativeModels, basic_astrodynamics::getAssociatedThirdBodyAcceleration(
                                accelerationDependentVariableSettings->accelerationModeType_  ) );
            }

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
    case relative_body_aerodynamic_orientation_angle_variable:
    {
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {
            std::string errorMessage = "Error when flight conditions for relative_body_aerodynamic_orientation_angle_variable output " +
                    bodyWithProperty + " has no flight conditions";
            throw std::runtime_error( errorMessage );
        }

        // Check input consistency.
        boost::shared_ptr< BodyAerodynamicAngleVariableSaveSettings > bodyAerodynamicAngleVariableSaveSettings =
                boost::dynamic_pointer_cast< BodyAerodynamicAngleVariableSaveSettings >( dependentVariableSettings );
        if( bodyAerodynamicAngleVariableSaveSettings == NULL )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type relative_body_aerodynamic_orientation_angle_variable";
            throw std::runtime_error( errorMessage );
        }

        variableFunction = boost::bind( &reference_frames::AerodynamicAngleCalculator::getAerodynamicAngle,
                                        bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                                        bodyAerodynamicAngleVariableSaveSettings->angle_ );
        break;
    }
    case total_aerodynamic_g_load_variable:
    {
        // Check input consistency
        boost::shared_ptr< SingleAccelerationDependentVariableSaveSettings > aerodynamicAccelerationSettings =
                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    basic_astrodynamics::aerodynamic, bodyWithProperty, secondaryBody, 0 );
        boost::function< Eigen::VectorXd( ) > aerodynamicAccelerationFunction =
                getVectorDependentVariableFunction( aerodynamicAccelerationSettings, bodyMap, stateDerivativeModels ).first;


        boost::function< double( Eigen::Vector3d )> functionToEvaluate =
                boost::bind( &aerodynamics::computeAerodynamicLoadFromAcceleration, _1 );
        variableFunction = boost::bind(
                    &evaluateReferenceFunction< double, Eigen::Vector3d >,
                    functionToEvaluate, aerodynamicAccelerationFunction );

        break;
    }
    case stagnation_point_heat_flux_dependent_variable:
    {

        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {
            std::string errorMessage = "Error no flight conditions available when requesting stagnation point heating output of" +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        boost::shared_ptr< aerodynamics::FlightConditions > flightConditions =
                bodyMap.at( bodyWithProperty )->getFlightConditions( );

        if( bodyMap.at( bodyWithProperty )->getVehicleSystems( ) == NULL )
        {
            std::string errorMessage = "Error, no vehicle systems available when requesting stagnation point heating output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        boost::shared_ptr< system_models::VehicleSystems > vehicleSystems =
                bodyMap.at( bodyWithProperty )->getVehicleSystems( );


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

        variableFunction = boost::bind(
                    &computeEquilibriumFayRiddellHeatFluxFromProperties,
                    flightConditions, vehicleSystems );

        break;
    }
    case local_temperature_dependent_variable:
    {
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {
            std::string errorMessage = "Error, no flight conditions available when requesting temperature output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }
        variableFunction = boost::bind( &aerodynamics::FlightConditions::getCurrentFreestreamTemperature,
                                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
        break;
    }
    case geodetic_latitude_dependent_variable:
    {
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {
            std::string errorMessage = "Error, no flight conditions available when requesting geodetic latitude output of " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        variableFunction = boost::bind( &aerodynamics::FlightConditions::getCurrentGeodeticLatitude,
                                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
        break;
    }
    case control_surface_deflection_dependent_variable:
    {
        if( bodyMap.at( bodyWithProperty )->getVehicleSystems( ) == NULL )
        {
            std::string errorMessage = "Error, no vehicle systems available when requesting control surface deflection output of " +
                    bodyWithProperty + "with surface" + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        variableFunction = boost::bind( &system_models::VehicleSystems::getCurrentControlSurfaceDeflection,
                                        bodyMap.at( bodyWithProperty )->getVehicleSystems( ),
                                        dependentVariableSettings->secondaryBody_ );
        break;
    }
    case total_mass_rate_dependent_variables:
    {
        // Retrieve model responsible for computing mass rate of requested bodies.
        boost::shared_ptr< BodyMassStateDerivative< StateScalarType, TimeType > > nBodyModel =
                getBodyMassStateDerivativeModelForBody( bodyWithProperty, stateDerivativeModels );
        variableFunction =
                boost::bind( &BodyMassStateDerivative< StateScalarType, TimeType >::getTotalMassRateForBody, nBodyModel,
                             bodyWithProperty );

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
                    &evaluateBivariateReferenceFunction< Eigen::Vector3d, Eigen::Vector3d >,
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
                    &evaluateBivariateReferenceFunction< Eigen::Vector3d, Eigen::Vector3d >,
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

            // Check if thirfd-body counterpart of acceleration is found
            if( listOfSuitableAccelerationModels.size( ) == 0 && basic_astrodynamics::isAccelerationDirectGravitational(
                        accelerationDependentVariableSettings->accelerationModeType_ ) )
            {
                listOfSuitableAccelerationModels = getAccelerationBetweenBodies(
                    accelerationDependentVariableSettings->associatedBody_,
                    accelerationDependentVariableSettings->secondaryBody_,
                    stateDerivativeModels, basic_astrodynamics::getAssociatedThirdBodyAcceleration(
                                accelerationDependentVariableSettings->accelerationModeType_  ) );
            }

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
    case aerodynamic_force_coefficients_dependent_variable:
    {
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {
            std::string errorMessage = "Error, no flight conditions available when requesting density output of aerodynamic force coefficients " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        variableFunction = boost::bind(
                    &aerodynamics::AerodynamicCoefficientInterface::getCurrentForceCoefficients,
                    bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicCoefficientInterface( ) );
        parameterSize = 3;

        break;
    }
    case aerodynamic_moment_coefficients_dependent_variable:
    {
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {

            std::string errorMessage = "Error, no flight conditions available when requesting density output of aerodynamic moment coefficients " +
                    bodyWithProperty + "w.r.t." + secondaryBody;
            throw std::runtime_error( errorMessage );
        }

        variableFunction = boost::bind(
                    &aerodynamics::AerodynamicCoefficientInterface::getCurrentMomentCoefficients,
                    bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicCoefficientInterface( ) );
        parameterSize = 3;

        break;
    }
    case rotation_matrix_to_body_fixed_frame_variable:
    {
        boost::function< Eigen::Quaterniond( ) > rotationFunction =
                boost::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame, bodyMap.at( bodyWithProperty ) );
        variableFunction = boost::bind( &getVectorRepresentationForRotationQuaternion, rotationFunction );
        parameterSize = 9;
        break;
    }
    case intermediate_aerodynamic_rotation_matrix_variable:
    {
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {
            std::string errorMessage= "Error, no flight conditions when creating dependent variable function of type intermediate_aerodynamic_rotation_matrix_variable";
            throw std::runtime_error( errorMessage );
        }

        // Check input consistency.
        boost::shared_ptr< IntermediateAerodynamicRotationVariableSaveSettings >
                intermediateAerodynamicRotationVariableSaveSettings =
                boost::dynamic_pointer_cast< IntermediateAerodynamicRotationVariableSaveSettings >(
                    dependentVariableSettings );
        if( intermediateAerodynamicRotationVariableSaveSettings == NULL )
        {
            std::string errorMessage= "Error, inconsistent inout when creating dependent variable function of type intermediate_aerodynamic_rotation_matrix_variable";
            throw std::runtime_error( errorMessage );
        }

        boost::function< Eigen::Quaterniond( ) > rotationFunction =
                boost::bind( &reference_frames::AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                             bodyMap.at( bodyWithProperty )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                             intermediateAerodynamicRotationVariableSaveSettings->baseFrame_,
                             intermediateAerodynamicRotationVariableSaveSettings->targetFrame_ );

        variableFunction = boost::bind( &getVectorRepresentationForRotationQuaternion, rotationFunction );
        parameterSize = 9;
        break;
    }
    case body_fixed_airspeed_based_velocity_variable:
    {
        if( bodyMap.at( bodyWithProperty )->getFlightConditions( ) == NULL )
        {
            std::string errorMessage= "Error, no flight conditions when creating dependent variable function of type current_airpeed_based_velocity_variable_variable";
            throw std::runtime_error( errorMessage );
        }

        variableFunction = boost::bind( &aerodynamics::FlightConditions::getCurrentAirspeedBasedVelocity,
                                        bodyMap.at( bodyWithProperty )->getFlightConditions( ) );
        parameterSize = 3;
        break;
    }
    case lvlh_to_inertial_frame_rotation_dependent_variable:
    {
        boost::function< Eigen::Vector6d( ) > vehicleStateFunction =
                boost::bind( &simulation_setup::Body::getState, bodyMap.at( dependentVariableSettings->associatedBody_ ) );
        boost::function< Eigen::Vector6d( ) > centralBodyStateFunction;

        if( ephemerides::isFrameInertial( dependentVariableSettings->secondaryBody_ ) )
        {
            centralBodyStateFunction =  boost::lambda::constant( Eigen::Vector6d::Zero( ) );
        }
        else
        {
            centralBodyStateFunction =
                    boost::bind( &simulation_setup::Body::getState, bodyMap.at( dependentVariableSettings->secondaryBody_ ) );
        }

        boost::function< Eigen::Matrix3d( ) > rotationFunction =
                boost::bind( &reference_frames::getVelocityBasedLvlhToInertialRotationFromFunctions,
                             vehicleStateFunction, centralBodyStateFunction, true );
        variableFunction = boost::bind(
                    &getVectorRepresentationForRotationMatrixFunction, rotationFunction );

        parameterSize = 9;

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
 *  \return Pair with function returning requested dependent variable values, and list variable names with start entries.
 *  NOTE: The environment and state derivative models need to
 *  be updated to current state and independent variable before computation is performed.
 */
template< typename TimeType = double, typename StateScalarType = double >
std::pair< boost::function< Eigen::VectorXd( ) >, std::map< int, std::string > > createDependentVariableListFunction(
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

    std::vector< std::string > doubleVariableList;
    std::vector< std::pair< std::string, int > > vectorVariableList;

    for( unsigned int i = 0; i < dependentVariables.size( ); i++ )
    {
        // Create double parameter
        if( getDependentVariableSize( dependentVariables.at( i )->variableType_ ) == 1 )
        {
            doubleFunctionList.push_back( getDoubleDependentVariableFunction(
                                              dependentVariables.at( i ),
                                              bodyMap, stateDerivativeModels ) );
            doubleVariableList.push_back( getDependentVariableId(
                        dependentVariables.at( i ) ) );
        }
        // Create vector parameter
        else
        {
            vectorFunctionList.push_back( getVectorDependentVariableFunction(
                                              dependentVariables.at( i ),
                                              bodyMap, stateDerivativeModels ) );
            vectorVariableList.push_back( std::make_pair( getDependentVariableId(
                        dependentVariables.at( i ) ), vectorFunctionList.at( vectorFunctionList.size( ) - 1 ).second ) );
        }
    }

    // Set list of variable ids/indices in correc otder.
    int totalVariableSize = 0;
    std::map< int, std::string > dependentVariableId;
    for( unsigned int i = 0; i < doubleVariableList.size( ); i++ )
    {
        dependentVariableId[ totalVariableSize ] = doubleVariableList.at( i );
        totalVariableSize++;
    }

    for( unsigned int i = 0; i < vectorFunctionList.size( ); i++ )
    {
        dependentVariableId[ totalVariableSize ] = vectorVariableList.at( i ).first;
        totalVariableSize += vectorVariableList.at( i ).second;
    }

    // Create function conatenating function results.
    return std::make_pair(
                boost::bind( &evaluateListOfFunctions, doubleFunctionList, vectorFunctionList, totalVariableSize ),
                dependentVariableId );
}


} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATIONOUTPUT_H
