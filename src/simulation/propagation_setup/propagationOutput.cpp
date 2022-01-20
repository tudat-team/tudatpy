/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/aerodynamics/aerodynamics.h"
#include "tudat/simulation/propagation_setup/propagationOutput.h"

namespace tudat
{

namespace propagators
{

//! Get the vector representation of a rotation matrix.
Eigen::VectorXd getVectorRepresentationForRotationMatrix(
        const Eigen::Matrix3d& currentRotationMatrix )
{
    Eigen::VectorXd vectorRepresentation = Eigen::VectorXd( 9 );
    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            vectorRepresentation( i * 3 + j ) = currentRotationMatrix( i, j );
        }
    }
    return vectorRepresentation;
}

//! Get the vector representation of a rotation matrix.
Eigen::VectorXd getVectorRepresentationForRotationMatrixFunction(
        const std::function< Eigen::Matrix3d( ) > rotationFunction )
{
    return getVectorRepresentationForRotationMatrix( rotationFunction( ) );
}

//! Get the vector representation of a quaternion.
Eigen::VectorXd getVectorRepresentationForRotationQuaternion(
        const std::function< Eigen::Quaterniond( ) > rotationFunction )
{
    return getVectorRepresentationForRotationMatrix( rotationFunction( ).toRotationMatrix( ) );
}

//! Get the 3x3 matrix representation from a vector with 9 entries
Eigen::Matrix3d getMatrixFromVectorRotationRepresentation(
        const Eigen::VectorXd vectorRepresentation )
{
    if( vectorRepresentation.rows( ) != 9 )
    {
        throw std::runtime_error( "Error when putting vector in matrix representation, size is incompatible" );
    }
    Eigen::Matrix3d currentRotationMatrix;
    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            currentRotationMatrix( i, j ) = vectorRepresentation( i * 3 + j );
        }
    }
    return currentRotationMatrix;
}


//! Get the quaternion formulation of an orthonormal matrix, from input of a vector with 9 entries corresponding to matrix
//! entries.
Eigen::Quaterniond getQuaternionFromVectorRotationRepresentation(
        const Eigen::VectorXd vectorRepresentation )
{
    return Eigen::Quaterniond( getMatrixFromVectorRotationRepresentation( vectorRepresentation ) );
}

//! Function to convert a matrix to the format used to save dependent variables
void getMatrixInOutputVectorRepresentation(
        const Eigen::MatrixXd& matrix, Eigen::VectorXd& vector )
{
    vector.setZero( matrix.rows( ) * matrix.cols( ) );
    for( int i = 0; i < matrix.rows( ); i++ )
    {
        vector.segment( i * matrix.cols( ), matrix.cols( ) ) =
              matrix.block( i, 0, 1, matrix.cols( ) ).transpose( );
    }
}

//! Function to convert a vector dependent variable output to its original matrix representation
void getOutputVectorInMatrixRepresentation(
        const Eigen::VectorXd& vector, Eigen::MatrixXd& matrix,
        const int rows, const int columns )
{
    if( rows * columns != vector.rows( ) )
    {
        throw std::runtime_error( "Error when getting matrix from output vector: sizes are incompatible" );
    }
    matrix.setZero( rows, columns );
    for( int i = 0; i < rows; i++ )
    {
        matrix.block( i, 0, 1, columns ) = vector.segment( i * columns, columns ).transpose( );
    }
}

//! Function to retrieve matrix block function output in vector representation
Eigen::VectorXd getVectorFunctionFromBlockFunction(
        const std::function< void( Eigen::Block< Eigen::MatrixXd > ) > blockFunction,
                                    const int numberOfRows, const int numberOfColumns )
{
    Eigen::MatrixXd matrixEvaluation = Eigen::MatrixXd::Zero( numberOfRows, numberOfColumns );
    blockFunction( matrixEvaluation.block( 0, 0, numberOfRows, numberOfColumns ) );

    Eigen::VectorXd vectorEvaluation;
    getMatrixInOutputVectorRepresentation( matrixEvaluation, vectorEvaluation );

    return vectorEvaluation;
}

//! Function to compute the Fay-Riddell equilibrium heat flux from body properties
double computeEquilibriumFayRiddellHeatFluxFromProperties(
        const std::shared_ptr< aerodynamics::AtmosphericFlightConditions > flightConditions,
        const std::shared_ptr< system_models::VehicleSystems > vehicleSystems )
{
    return aerodynamics::computeEquilibriumFayRiddellHeatFlux(
                flightConditions->getCurrentDensity( ), flightConditions->getCurrentAirspeed( ),
                flightConditions->getCurrentFreestreamTemperature( ), flightConditions->getCurrentMachNumber( ),
                vehicleSystems->getNoseRadius( ), vehicleSystems->getWallEmissivity( ) );
}


//! Function to return a vector containing only one value given by doubleFunction
Eigen::VectorXd getVectorFromDoubleFunction( const std::function< double( ) >& doubleFunction )
{
    Eigen::VectorXd vector( 1 );
    vector << doubleFunction( );
    return vector;
}

//! Function to evaluate a set of vector-returning functions and concatenate the results.
Eigen::VectorXd evaluateListOfVectorFunctions(
        const std::vector< std::pair< std::function< Eigen::VectorXd( ) >, int > > vectorFunctionList,
        const int totalSize )
{
    Eigen::VectorXd variableList = Eigen::VectorXd::Zero( totalSize );
    int currentIndex = 0;

    for( std::pair< std::function< Eigen::VectorXd( ) >, int > vectorFunction: vectorFunctionList )
    {
        variableList.segment( currentIndex, vectorFunction.second ) = vectorFunction.first( );
        currentIndex += vectorFunction.second;
    }

    // Check consistency with input
    if( currentIndex != totalSize )
    {
        std::string errorMessage = "Error when evaluating lists of functions, sizes are inconsistent: " +
                std::to_string( currentIndex ) + " and " +
                std::to_string( totalSize );
        throw std::runtime_error( errorMessage );
    }

    return variableList;
}

Eigen::VectorXd getNormsOfAccelerationDifferencesFromLists(
                       const std::function< Eigen::VectorXd( ) > firstAccelerationFunction,
                       const std::function< Eigen::VectorXd( ) > secondAccelerationFunction )
{
    Eigen::VectorXd firstAcceleration = firstAccelerationFunction( );
    Eigen::VectorXd secondAcceleration = secondAccelerationFunction( );

    if( firstAcceleration.rows( ) != secondAcceleration.rows( ) )
    {
        throw std::runtime_error( "Error when computing acceleration difference norms, inputs are inconsistent." );
    }

    if( firstAcceleration.rows( ) % 3 != 0 )
    {
        throw std::runtime_error( "Error when computing acceleration difference norms, input size is inconsistent." );
    }

    Eigen::VectorXd accelerationDifference = Eigen::VectorXd::Zero(
                firstAcceleration.rows( ) / 3 );
    for( int i = 0; i < accelerationDifference.rows( ); i++ )
    {
        accelerationDifference( i ) =
                ( firstAcceleration.segment( i * 3, 3 ) - secondAcceleration.segment( i * 3, 3 ) ).norm( );
    }

    return accelerationDifference;
}

//! Funtion to get the size of a dependent variable save settings
int getDependentVariableSaveSize(
        const std::shared_ptr< SingleDependentVariableSaveSettings >& singleDependentVariableSaveSettings )
{
    if ( singleDependentVariableSaveSettings->componentIndex_ >= 0 )
    {
        return 1;
    }
    else
    {
        return getDependentVariableSize(  singleDependentVariableSaveSettings );
    }
}

//! Funtion to get the size of a dependent variable
int getDependentVariableSize(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings )
{
    int variableSize = -1;
    switch( dependentVariableSettings->dependentVariableType_ )
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
    case aerodynamic_force_coefficients_dependent_variable:
        variableSize = 3;
        break;
    case aerodynamic_moment_coefficients_dependent_variable:
        variableSize = 3;
        break;
    case inertial_to_body_fixed_rotation_matrix_variable:
        variableSize = 9;
        break;
    case intermediate_aerodynamic_rotation_matrix_variable:
        variableSize = 9;
        break;
    case relative_body_aerodynamic_orientation_angle_variable:
        variableSize = 1;
        break;
    case body_fixed_airspeed_based_velocity_variable:
        variableSize = 3;
        break;
    case body_fixed_groundspeed_based_velocity_variable:
        variableSize = 3;
        break;
    case total_aerodynamic_g_load_variable:
        variableSize = 1;
        break;
    case stagnation_point_heat_flux_dependent_variable:
        variableSize = 1;
        break;
    case local_temperature_dependent_variable:
        variableSize = 1;
        break;
    case local_dynamic_pressure_dependent_variable:
        variableSize = 1;
        break;
    case local_aerodynamic_heat_rate_dependent_variable:
        variableSize = 1;
        break;
    case geodetic_latitude_dependent_variable:
        variableSize = 1;
        break;
    case control_surface_deflection_dependent_variable:
        variableSize = 1;
        break;
    case total_mass_rate_dependent_variables:
        variableSize = 1;
        break;
    case tnw_to_inertial_frame_rotation_dependent_variable:
        variableSize = 9;
        break;
    case rsw_to_inertial_frame_rotation_dependent_variable:
        variableSize = 9;
        break;
    case periapsis_altitude_dependent_variable:
        variableSize = 1;
        break;
    case total_torque_dependent_variable:
        variableSize = 3;
        break;
    case single_torque_dependent_variable:
        variableSize = 3;
        break;
    case total_torque_norm_dependent_variable:
        variableSize = 1;
        break;
    case single_torque_norm_dependent_variable:
        variableSize = 3;
        break;
    case keplerian_state_dependent_variable:
        variableSize = 6;
        break;
    case spherical_harmonic_acceleration_norm_terms_dependent_variable:
    {
        std::shared_ptr< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >
                sphericalHarmonicAccelerationTermsDependentVariableSaveSettings =
                std::dynamic_pointer_cast< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >(
                                    dependentVariableSettings );
        if( sphericalHarmonicAccelerationTermsDependentVariableSaveSettings == nullptr )
        {
             std::string errorMessage = "Error, input for spherical_harmonic_acceleration_terms_dependent_variable inconsistent when getting parameter size ";
             throw std::runtime_error( errorMessage );
        }
        else
        {
            variableSize = sphericalHarmonicAccelerationTermsDependentVariableSaveSettings->componentIndices_.size( );
        }
        break;
    }
    case spherical_harmonic_acceleration_terms_dependent_variable:
    {
        std::shared_ptr< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >
                sphericalHarmonicAccelerationTermsDependentVariableSaveSettings =
                std::dynamic_pointer_cast< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >(
                                    dependentVariableSettings );
        if( sphericalHarmonicAccelerationTermsDependentVariableSaveSettings == nullptr )
        {
             std::string errorMessage = "Error, input for spherical_harmonic_acceleration_terms_dependent_variable inconsistent when getting parameter size ";
             throw std::runtime_error( errorMessage );
        }
        else
        {
            variableSize = 3 * sphericalHarmonicAccelerationTermsDependentVariableSaveSettings->componentIndices_.size( );
        }
        break;
    }
    case modified_equinocial_state_dependent_variable:
        variableSize = 6;
        break;
    case body_fixed_relative_cartesian_position:
        variableSize = 3;
        break;
    case body_fixed_relative_spherical_position:
        variableSize = 3;
        break;
    case euler_angles_to_body_fixed_313:
        variableSize = 3;
        break;
    case total_gravity_field_variation_acceleration:
        variableSize = 3;
        break;
    case single_gravity_field_variation_acceleration:
        variableSize = 3;
        break;
    case single_gravity_field_variation_acceleration_terms:
    {
        if( std::dynamic_pointer_cast< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
                    dependentVariableSettings ) == nullptr )
        {
             std::string errorMessage = "Error, input for single_gravity_field_variation_acceleration_terms inconsistent when getting parameter size ";
             throw std::runtime_error( errorMessage );
        }
        else
        {
            variableSize = 3 * std::dynamic_pointer_cast< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
                        dependentVariableSettings )->componentIndices_.size( );
        }
        break;
    }
    case acceleration_partial_wrt_body_translational_state:
        variableSize = 18;
        break;
    case current_body_mass_dependent_variable:
        variableSize = 1;
        break;
    case radiation_pressure_coefficient_dependent_variable:
        variableSize = 1;
        break;
    case custom_dependent_variable:
        if( std::dynamic_pointer_cast< CustomDependentVariableSaveSettings >(
                    dependentVariableSettings ) == nullptr )
        {
            std::string errorMessage = "Error, input for custom dependent variable parameter size ";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            variableSize = 3 * std::dynamic_pointer_cast< CustomDependentVariableSaveSettings >(
                        dependentVariableSettings )->dependentVariableSize_;
        }
        break;
    default:
        std::string errorMessage = "Error, did not recognize dependent variable size of type: " +
                std::to_string( dependentVariableSettings->dependentVariableType_ );
        throw std::runtime_error( errorMessage );
    }
    return variableSize;
}

template std::pair< std::function< Eigen::VectorXd( ) >, std::map< int, std::string > > createDependentVariableListFunction< double, double >(
        const std::shared_ptr< DependentVariableSaveSettings > saveSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const std::unordered_map< IntegratedStateType,
        std::vector< std::shared_ptr< SingleStateTypeDerivative< double, double > > > >& stateDerivativeModels );

//template std::pair< std::function< Eigen::VectorXd( ) >, int > getVectorDependentVariableFunction< double, double >(
//        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const std::unordered_map< IntegratedStateType,
//        std::vector< std::shared_ptr< SingleStateTypeDerivative< double, double > > > >& stateDerivativeModels );

//template std::function< double( ) > getDoubleDependentVariableFunction< double, double >(
//        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const std::unordered_map< IntegratedStateType,
//        std::vector< std::shared_ptr< SingleStateTypeDerivative< double, double > > > >& stateDerivativeModels );

} // namespace propagators

} // namespace tudat
