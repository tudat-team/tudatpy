/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"

namespace tudat
{

namespace propagators
{

//! Function to get a string representing a 'named identification' of a variable type
std::string getVariableName( const VariableType variableType )
{
    switch ( variableType )
    {
    case independentVariable:
        return "Independent variable ";
    case cpuTimeVariable:
        return "Cumulative computation time variable ";
    case stateVariable:
        return "Integrated state ";
    case dependentVariable:
        return "Dependent variable ";
    case stateTransitionMatrix:
        return "State transition matrix ";
    case sensitivityMatrix:
        return "Sensitivity matrix ";
    default:
        throw std::runtime_error( "Error, variable " +
                                  std::to_string( variableType ) +
                                  "not found when retrieving parameter name " );
    }
}

//! Function to get a string representing a 'named identification' of a variable
std::string getVariableId( const std::shared_ptr< VariableSettings > variableSettings )
{
    std::shared_ptr< SingleDependentVariableSaveSettings > singleDependentVariableSaveSettings =
            std::dynamic_pointer_cast< SingleDependentVariableSaveSettings >( variableSettings );
    if ( singleDependentVariableSaveSettings )
    {
        return getDependentVariableId( singleDependentVariableSaveSettings );
    }
    else
    {
        return getVariableName( variableSettings->variableType_ );
    }
}


//! Function to get a string representing a 'named identification' of a dependent variable type
std::string getDependentVariableName(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings )
{

    PropagationDependentVariables propagationDependentVariables = dependentVariableSettings->dependentVariableType_;
    std::string variableName = "";
    switch( propagationDependentVariables )
    {
    case mach_number_dependent_variable:
        variableName = "Mach number ";
        break;
    case altitude_dependent_variable:
        variableName = "Altitude ";
        break;
    case airspeed_dependent_variable:
        variableName = "Airspeed ";
        break;
    case local_density_dependent_variable:
        variableName = "Density ";
        break;
    case relative_speed_dependent_variable:
        variableName = "Relative speed ";
        break;
    case relative_position_dependent_variable:
        variableName = "Relative position ";
        break;
    case relative_distance_dependent_variable:
        variableName = "Relative distance ";
        break;
    case relative_velocity_dependent_variable:
        variableName = "Relative velocity ";
        break;
    case radiation_pressure_dependent_variable:
        variableName = "Radiation pressure ";
        break;
    case total_acceleration_norm_dependent_variable:
        variableName = "Total acceleration norm ";
        break;
    case single_acceleration_norm_dependent_variable:
        variableName = "Single acceleration norm of type ";
        break;
    case total_acceleration_dependent_variable:
        variableName = "Total acceleration in inertial frame ";
        break;
    case single_acceleration_dependent_variable:
        variableName = "Single acceleration in inertial frame of type ";
        break;
    case aerodynamic_force_coefficients_dependent_variable:
        variableName = "Aerodynamic force coefficients ";
        break;
    case aerodynamic_moment_coefficients_dependent_variable:
        variableName = "Aerodynamic moment coefficients ";
        break;
    case aerodynamic_control_surface_free_force_coefficients_dependent_variable:
        variableName = "Aerodynamic control-surface free force coefficients ";
        break;
    case aerodynamic_control_surface_free_moment_coefficients_dependent_variable:
        variableName = "Aerodynamic control-surface free moment coefficients ";
        break;
    case aerodynamic_control_surface_force_coefficients_increment_dependent_variable:
        variableName = "Aerodynamic force coefficients control surface increment ";
        break;
    case aerodynamic_control_surface_moment_coefficients_increment_dependent_variable:
        variableName = "Aerodynamic moment coefficients control surface increment ";
        break;
    case inertial_to_body_fixed_rotation_matrix_variable:
        variableName = "Rotation matrix to body-fixed frame ";
        break;
    case intermediate_aerodynamic_rotation_matrix_variable:
        variableName = "Rotation matrix from ";
        break;
    case relative_body_aerodynamic_orientation_angle_variable:
    {
        std::shared_ptr< BodyAerodynamicAngleVariableSaveSettings > angleDependentVariableSettings =
                std::dynamic_pointer_cast< BodyAerodynamicAngleVariableSaveSettings >( dependentVariableSettings );
        if( angleDependentVariableSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting dependent variable type string, expected body angle type" );
        }
        else
        {
            reference_frames::AerodynamicsReferenceFrameAngles angleToSave = angleDependentVariableSettings->angle_;
            if( angleToSave == reference_frames::latitude_angle ||
                    angleToSave == reference_frames::longitude_angle )
            {
                variableName = "Spherical position angle ";
            }
            else if( angleToSave == reference_frames::heading_angle ||
                     angleToSave == reference_frames::flight_path_angle )
            {
                variableName = "Velocity orientation angle  ";
            }
            else if( angleToSave == reference_frames::angle_of_attack ||
                     angleToSave == reference_frames::angle_of_sideslip ||
                     angleToSave == reference_frames::bank_angle  )
            {
                variableName = "Aerodynamic body orientation angle ";
            }
        }
        break;
    }
    case body_fixed_airspeed_based_velocity_variable:
        variableName = "Airspeed-based velocity ";
        break;
    case body_fixed_groundspeed_based_velocity_variable:
        variableName = "Groundspeed-based velocity ";
        break;
    case total_aerodynamic_g_load_variable:
        variableName = "Aerodynamic g-load ";
        break;
    case stagnation_point_heat_flux_dependent_variable:
        variableName = "Stagnation-point heat flux ";
        break;
    case local_temperature_dependent_variable:
        variableName = "Local freestream temperature ";
        break;
    case local_dynamic_pressure_dependent_variable:
        variableName = "Local dynamic pressure ";
        break;
    case local_aerodynamic_heat_rate_dependent_variable:
        variableName = "Local aerodynamic heat rate ";
        break;
    case geodetic_latitude_dependent_variable:
        variableName = "Geodetic latitude ";
        break;
    case control_surface_deflection_dependent_variable:
        variableName = "Control Surface Deflection ";
        break;
    case total_mass_rate_dependent_variables:
        variableName = "Body mass rate ";
        break;
    case tnw_to_inertial_frame_rotation_dependent_variable:
        variableName = "TNW to inertial frame rotation matrix ";
        break;
    case rsw_to_inertial_frame_rotation_dependent_variable:
        variableName = "RSW to inertial frame rotation matrix ";
        break;
    case periapsis_altitude_dependent_variable:
        variableName = "Periapsis altitude ";
        break;
    case apoapsis_altitude_dependent_variable:
        variableName = "Apoapsis altitude ";
        break;
    case single_torque_dependent_variable:
        variableName = "Single torque in body-fixed frame of type ";
        break;
    case total_torque_dependent_variable:
        variableName = "Total torque in body-fixed frame ";
        break;
    case single_torque_norm_dependent_variable:
        variableName = "Single torque norm in body-fixed frame of type ";
        break;
    case total_torque_norm_dependent_variable:
        variableName = "Total torque norm in body-fixed frame ";
        break;
    case keplerian_state_dependent_variable:
        variableName = "Kepler elements ";
        break;
    case modified_equinocial_state_dependent_variable:
        variableName = "Modified equinoctial elements ";
        break;
    case spherical_harmonic_acceleration_norm_terms_dependent_variable:
        variableName = "Spherical harmonic acceleration term norms ";
        break;
    case spherical_harmonic_acceleration_terms_dependent_variable:
        variableName = "Spherical harmonic acceleration terms ";
        break;
    case body_fixed_relative_cartesian_position:
        variableName = "Body-fixed relative Cartesian position ";
        break;
    case body_fixed_relative_spherical_position:
        variableName = "Body-fixed relative spherical position ";
        break;
    case euler_angles_to_body_fixed_313:
        variableName = "313 Euler angles to body-fixed frame ";
        break;
    case total_gravity_field_variation_acceleration:
        variableName = "Total time-variable gravity field acceleration correction ";
        break;
    case single_gravity_field_variation_acceleration:
        variableName = "Single-source time-variable gravity field acceleration correction ";
        break;
    case single_gravity_field_variation_acceleration_terms:
        variableName = "Single-source time-variable gravity field per-term acceleration correction ";
        break;
    case total_spherical_harmonic_cosine_coefficient_variation:
        variableName = "Total spherical-harmonic cosine coefficient variation ";
        break;
    case total_spherical_harmonic_sine_coefficient_variation:
        variableName = "Total spherical-harmonic sine coefficient variation ";
        break;
    case acceleration_partial_wrt_body_translational_state:
    {
        std::shared_ptr< AccelerationPartialWrtStateSaveSettings > partialDependentVariableSettings =
                std::dynamic_pointer_cast< AccelerationPartialWrtStateSaveSettings >( dependentVariableSettings );
        if( partialDependentVariableSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting dependent variable type string, expected acceleration partial type" );
        }
        else
        {
            std::string accelerationName = basic_astrodynamics::getAccelerationModelName(
                        partialDependentVariableSettings->accelerationModelType_ );
            accelerationName[ 0 ] = std::toupper( accelerationName[ 0 ] );
            variableName =  accelerationName + "acceleration partial";
        }
        break;
    }
    case current_body_mass_dependent_variable:
        variableName = "Current body mass ";
        break;
    case radiation_pressure_coefficient_dependent_variable:
        variableName = "Radiation pressure coefficient ";
        break;
    case custom_dependent_variable:
        variableName = "Custom dependent variable ";
        break;
    case gravity_field_potential_dependent_variable:
        variableName = "Gravity field potential ";
        break;
    case gravity_field_laplacian_of_potential_dependent_variable:
        variableName = "Gravity field laplacian of potential";
        break;
    case total_acceleration_partial_wrt_body_translational_state:
        variableName = "Total acceleration partial w.r.t. translational state";
        break;
    case minimum_constellation_distance:
        variableName = "Minimum instantaneous constellation distance";
        break;
    case minimum_constellation_ground_station_distance:
        variableName = "Minimum instantaneous ground station visible constellation distance";
        break;
    case body_center_of_mass:
        variableName = "Body center of mass";
        break;
    case body_inertia_tensor:
        variableName = "Body inertia tensor";
        break;
    case received_irradiance:
        variableName = "Received irradiance";
        break;
    case received_fraction:
        variableName = "Received fraction of irradiance";
        break;
    case visible_and_emitting_source_panel_count:
        variableName = "Number of visible and emitting source panels";
        break;
    case visible_source_area:
        variableName = "Visible area";
        break;
    case vehicle_panel_inertial_surface_normals:
        variableName = "Vehicle surface panel inertial surface normals";
        break;
    case vehicle_panel_body_fixed_surface_normals:
        variableName = "Vehicle surface panel body-fixed surface normals";
        break;
    case vehicle_surface_panel_radiation_pressure_force:
        variableName = "Per-target panel vehicle radiation pressure force";
        break;
    case paneled_radiation_source_per_panel_irradiance:
        variableName = "Per-source panel radiation pressure irradiance";
        break;
    case paneled_radiation_source_geometry:
        variableName = "Per-source panel radiation pressure geometry";
        break;
    case nrlmsise_input_data:
        variableName = "NRLMSISE00 input data vector";
        break;
    default:
        std::string errorMessage = "Error, dependent variable " +
                std::to_string( propagationDependentVariables ) +
                "not found when retrieving parameter name ";
        throw std::runtime_error( errorMessage );
    }
    return variableName;
}


//! Function to get a string representing a 'named identification' of a dependent variable
std::string getDependentVariableId(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings )
{
    std::string variableId = getDependentVariableName( dependentVariableSettings );

    if( ( dependentVariableSettings->dependentVariableType_ == single_acceleration_dependent_variable ) ||
            ( dependentVariableSettings->dependentVariableType_ == single_acceleration_norm_dependent_variable ) )
    {
        std::shared_ptr< SingleAccelerationDependentVariableSaveSettings > accelerationDependentVariableSettings =
                std::dynamic_pointer_cast< SingleAccelerationDependentVariableSaveSettings >( dependentVariableSettings );
        if( accelerationDependentVariableSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting dependent variable ID, input is inconsistent (acceleration type )" );
        }
        else
        {
            variableId += basic_astrodynamics::getAccelerationModelName(
                        accelerationDependentVariableSettings->accelerationModelType_ );
        }
    }
    else if( ( dependentVariableSettings->dependentVariableType_ == single_torque_dependent_variable ) ||
             ( dependentVariableSettings->dependentVariableType_ == single_torque_norm_dependent_variable ) )
    {
        std::shared_ptr< SingleTorqueDependentVariableSaveSettings > torqueDependentVariableSettings =
                std::dynamic_pointer_cast< SingleTorqueDependentVariableSaveSettings >( dependentVariableSettings );
        if( torqueDependentVariableSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting dependent variable ID, input is inconsistent (torque type )" );
        }
        else
        {
            variableId += basic_astrodynamics::getTorqueModelName(
                        torqueDependentVariableSettings->torqueModelType_ );
        }
    }
    else if( dependentVariableSettings->dependentVariableType_ == intermediate_aerodynamic_rotation_matrix_variable )
    {
        std::shared_ptr< IntermediateAerodynamicRotationVariableSaveSettings > rotationDependentVariableSettings =
                std::dynamic_pointer_cast< IntermediateAerodynamicRotationVariableSaveSettings >( dependentVariableSettings );
        if( rotationDependentVariableSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting dependent variable ID, input is inconsistent (rotation matrix)" );
        }
        else
        {
            variableId +=
                    reference_frames::getAerodynamicFrameName( rotationDependentVariableSettings->baseFrame_ ) + "to " +
                    reference_frames::getAerodynamicFrameName( rotationDependentVariableSettings->targetFrame_ );
        }
    }
    else if( dependentVariableSettings->dependentVariableType_ == relative_body_aerodynamic_orientation_angle_variable )
    {
        std::shared_ptr< BodyAerodynamicAngleVariableSaveSettings > angleDependentVariableSettings =
                std::dynamic_pointer_cast< BodyAerodynamicAngleVariableSaveSettings >( dependentVariableSettings );
        if( angleDependentVariableSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting dependent variable ID, input is inconsistent (angle)" );
        }
        else
        {
            variableId +=
                    reference_frames::getAerodynamicAngleName( angleDependentVariableSettings->angle_ );
        }
    }

    if( ( dependentVariableSettings->dependentVariableType_ == single_acceleration_dependent_variable ) ||
            ( dependentVariableSettings->dependentVariableType_ == single_acceleration_norm_dependent_variable ) ||
            ( dependentVariableSettings->dependentVariableType_ == spherical_harmonic_acceleration_terms_dependent_variable ) ||
            ( dependentVariableSettings->dependentVariableType_ == spherical_harmonic_acceleration_norm_terms_dependent_variable )  ||
            ( dependentVariableSettings->dependentVariableType_ == acceleration_partial_wrt_body_translational_state ) )
    {
        variableId += ", acting on " + dependentVariableSettings->associatedBody_;
        if( dependentVariableSettings->secondaryBody_ != dependentVariableSettings->associatedBody_ )
        {
            variableId += ", exerted by " + dependentVariableSettings->secondaryBody_;
        }
    }
    else if ((dependentVariableSettings->dependentVariableType_ == visible_and_emitting_source_panel_count) ||
             (dependentVariableSettings->dependentVariableType_ == visible_source_area))
    {
        variableId += " of " + dependentVariableSettings->secondaryBody_
                + " as seen from " + dependentVariableSettings->associatedBody_;
    }
    else if ((dependentVariableSettings->dependentVariableType_ == received_irradiance) ||
            (dependentVariableSettings->dependentVariableType_ == received_fraction))
    {
        variableId += " at " + dependentVariableSettings->associatedBody_
                + " due to " + dependentVariableSettings->secondaryBody_;
    }
    else
    {
        variableId += "of " + dependentVariableSettings->associatedBody_;
        if( dependentVariableSettings->secondaryBody_ != "" )
        {
            variableId += " w.r.t. " + dependentVariableSettings->secondaryBody_;
        }

    }

    if( dependentVariableSettings->dependentVariableType_ == acceleration_partial_wrt_body_translational_state )
    {
        std::shared_ptr< AccelerationPartialWrtStateSaveSettings > partialDependentVariableSettings =
                std::dynamic_pointer_cast< AccelerationPartialWrtStateSaveSettings >( dependentVariableSettings );
        if( partialDependentVariableSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting dependent variable type full string, expected acceleration partial type" );
        }
        else
        {
            variableId += " w.r.t. body translational state of " + partialDependentVariableSettings->derivativeWrtBody_;
        }
    }


    if( dependentVariableSettings->dependentVariableType_ == total_acceleration_partial_wrt_body_translational_state )
    {
        std::shared_ptr< TotalAccelerationPartialWrtStateSaveSettings > partialDependentVariableSettings =
                std::dynamic_pointer_cast< TotalAccelerationPartialWrtStateSaveSettings >( dependentVariableSettings );
        if( partialDependentVariableSettings == nullptr )
        {
            throw std::runtime_error( "Error when getting dependent variable type full string, expected full acceleration partial type" );
        }
        else
        {
            variableId += " w.r.t. body translational state of " + partialDependentVariableSettings->derivativeWrtBody_;
        }
    }

    return variableId;
}

}

}
