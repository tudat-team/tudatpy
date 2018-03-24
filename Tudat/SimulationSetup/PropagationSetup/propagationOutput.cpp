/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamics.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationOutput.h"

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
        const boost::function< Eigen::Matrix3d( ) > rotationFunction )
{
    return getVectorRepresentationForRotationMatrix( rotationFunction( ) );
}

//! Get the vector representation of a quaternion.
Eigen::VectorXd getVectorRepresentationForRotationQuaternion(
        const boost::function< Eigen::Quaterniond( ) > rotationFunction )
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

//! Function to compute the Fay-Riddell equilibrium heat flux from body properties
double computeEquilibriumFayRiddellHeatFluxFromProperties(
        const boost::shared_ptr< aerodynamics::FlightConditions > flightConditions,
        const boost::shared_ptr< system_models::VehicleSystems > vehicleSystems )
{
    return aerodynamics::computeEquilibriumFayRiddellHeatFlux(
                flightConditions->getCurrentDensity( ), flightConditions->getCurrentAirspeed( ),
                flightConditions->getCurrentFreestreamTemperature( ), flightConditions->getCurrentMachNumber( ),
                vehicleSystems->getNoseRadius( ), vehicleSystems->getWallEmissivity( ) );
}


//! Function to return a vector containing only one value given by doubleFunction
Eigen::VectorXd getVectorFromDoubleFunction( const boost::function< double( ) >& doubleFunction )
{
    Eigen::VectorXd vector( 1 );
    vector << doubleFunction( );
    return vector;
}

//! Function to evaluate a set of vector-returning functions and concatenate the results.
Eigen::VectorXd evaluateListOfVectorFunctions(
        const std::vector< std::pair< boost::function< Eigen::VectorXd( ) >, int > > vectorFunctionList,
        const int totalSize )
{
    Eigen::VectorXd variableList = Eigen::VectorXd::Zero( totalSize );
    int currentIndex = 0;

    for( std::pair< boost::function< Eigen::VectorXd( ) >, int > vectorFunction: vectorFunctionList )
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

//! Funtion to get the size of a dependent variable save settings
int getDependentVariableSaveSize(
        const boost::shared_ptr< SingleDependentVariableSaveSettings >& singleDependentVariableSaveSettings )
{
    if ( singleDependentVariableSaveSettings->componentIndex_ >= 0 )
    {
        return 1;
    }
    else
    {
        return getDependentVariableSize( singleDependentVariableSaveSettings->dependentVariableType_ );
    }
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
    case aerodynamic_force_coefficients_dependent_variable:
        variableSize = 3;
        break;
    case aerodynamic_moment_coefficients_dependent_variable:
        variableSize = 3;
        break;
    case rotation_matrix_to_body_fixed_frame_variable:
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
    case geodetic_latitude_dependent_variable:
        variableSize = 1;
        break;
    case control_surface_deflection_dependent_variable:
        variableSize = 1;
        break;
    case total_mass_rate_dependent_variables:
        variableSize = 1;
        break;
    case lvlh_to_inertial_frame_rotation_dependent_variable:
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
    case modified_equinocial_state_dependent_variable:
        variableSize = 6;
        break;
    default:
        std::string errorMessage = "Error, did not recognize dependent variable size of type: " +
                std::to_string( dependentVariableSettings );
        throw std::runtime_error( errorMessage );
    }
    return variableSize;
}


} // namespace propagators

} // namespace tudat
