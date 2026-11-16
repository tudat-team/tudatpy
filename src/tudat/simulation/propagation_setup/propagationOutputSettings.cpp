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

#include <sstream>

namespace tudat
{

namespace propagators
{

//! Function to get a string representing a 'named identification' of a variable type
std::string getVariableName( const VariableType variableType )
{
    switch( variableType )
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
            throw std::runtime_error( "Error, variable " + std::to_string( variableType ) + "not found when retrieving parameter name " );
    }
}

//! Function to get a string representing a 'named identification' of a variable
std::string getVariableId( const std::shared_ptr< VariableSettings > variableSettings )
{
    std::shared_ptr< SingleDependentVariableSaveSettings > singleDependentVariableSaveSettings =
            std::dynamic_pointer_cast< SingleDependentVariableSaveSettings >( variableSettings );
    if( singleDependentVariableSaveSettings )
    {
        return getDependentVariableId( singleDependentVariableSaveSettings );
    }
    else
    {
        return getVariableName( variableSettings->variableType_ );
    }
}

//! Function to get a string representing a 'named identification' of a dependent variable type
std::string getDependentVariableName( const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings )
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
        case vehicle_part_rotation_matrix_dependent_variable:
            variableName = "Vehicle part rotation matrix ";
            break;
        case intermediate_aerodynamic_rotation_matrix_variable:
            variableName = "Rotation matrix from ";
            break;
        case relative_body_aerodynamic_orientation_angle_variable: {
            std::shared_ptr< BodyAerodynamicAngleVariableSaveSettings > angleDependentVariableSettings =
                    std::dynamic_pointer_cast< BodyAerodynamicAngleVariableSaveSettings >( dependentVariableSettings );
            if( angleDependentVariableSettings == nullptr )
            {
                throw std::runtime_error( "Error when getting dependent variable type string, expected body angle type" );
            }
            else
            {
                reference_frames::AerodynamicsReferenceFrameAngles angleToSave = angleDependentVariableSettings->angle_;
                if( angleToSave == reference_frames::latitude_angle || angleToSave == reference_frames::longitude_angle )
                {
                    variableName = "Spherical position angle ";
                }
                else if( angleToSave == reference_frames::heading_angle || angleToSave == reference_frames::flight_path_angle )
                {
                    variableName = "Velocity orientation angle  ";
                }
                else if( angleToSave == reference_frames::angle_of_attack || angleToSave == reference_frames::angle_of_sideslip ||
                         angleToSave == reference_frames::bank_angle )
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
        case local_wind_velocity_dependent_variable:
            variableName = "Local wind velocity ";
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
        case acceleration_partial_wrt_body_translational_state: {
            std::shared_ptr< AccelerationPartialWrtStateSaveSettings > partialDependentVariableSettings =
                    std::dynamic_pointer_cast< AccelerationPartialWrtStateSaveSettings >( dependentVariableSettings );
            if( partialDependentVariableSettings == nullptr )
            {
                throw std::runtime_error( "Error when getting dependent variable type string, expected acceleration partial type" );
            }
            else
            {
                std::string accelerationName =
                        basic_astrodynamics::getAccelerationModelName( partialDependentVariableSettings->accelerationModelType_ );
                accelerationName[ 0 ] = std::toupper( accelerationName[ 0 ] );
                variableName = accelerationName + "acceleration partial";
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
        case acceleration_derivative_partial_wrt_parameter:
            variableName = "Acceleration derivative partial w.r.t. parameter";
            break;
        case total_acceleration_derivative_partial_wrt_parameter:
            variableName = "Total acceleration derivative partial w.r.t. parameter";
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
        case illuminated_panel_fraction:
            variableName = "Illuminated fractions";
            break;
        case cross_section_change:
            variableName = "Cross section change";
            break;
        case actual_cross_section:
            variableName = "Actual cross section";
            break;
        case full_body_paneled_geometry:
            variableName = "Full body paneled geometry";
            break;
        case aerodynamic_coefficients:
            variableName = "Aerodynamic coefficients";
            break;
        case solar_longitude:
            variableName = "Solar longitude";
            break;
        case number_density:
            variableName = "Number density";
            break;
        case proper_time_rate_kinematic_term:
            variableName = "Proper-time-rate kinematic term -v^2/(2 c^2)";
            break;
        case proper_time_rate_potential_term:
            variableName = "Proper-time-rate potential term -U_ext/c^2";
            break;
        default:
            std::string errorMessage = "Error, dependent variable " + std::to_string( propagationDependentVariables ) +
                    "not found when retrieving parameter name ";
            throw std::runtime_error( errorMessage );
    }
    return variableName;
}

//! Function to cast dependent variable settings to the derived class that is required to create its ID
template< typename SettingsType >
static std::shared_ptr< SettingsType > getDerivedSettingsForId(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings )
{
    std::shared_ptr< SettingsType > derivedSettings = std::dynamic_pointer_cast< SettingsType >( dependentVariableSettings );
    if( derivedSettings == nullptr )
    {
        throw std::runtime_error( "Error when getting dependent variable ID, input is inconsistent for " +
                                  getDependentVariableName( dependentVariableSettings ) );
    }
    return derivedSettings;
}

//! Function to get the dependent variable ID sub-string for a list of spherical harmonic degree/order terms
static std::string getDegreeOrderTermsIdString( const std::vector< std::pair< int, int > >& degreeOrderTerms )
{
    std::string idString = ", degree/order";
    for( const std::pair< int, int >& degreeOrderTerm : degreeOrderTerms )
    {
        idString += " (" + std::to_string( degreeOrderTerm.first ) + "," + std::to_string( degreeOrderTerm.second ) + ")";
    }
    return idString;
}

//! Function to get the dependent variable ID sub-string for a list of bodies
static std::string getBodyListIdString( const std::vector< std::string >& bodies )
{
    std::string idString = ", bodies";
    for( const std::string& body : bodies )
    {
        idString += " " + body;
    }
    return idString;
}

//! Function to get the dependent variable ID sub-string for a single gravity field variation
static std::string getGravityFieldVariationIdString( const gravitation::BodyDeformationTypes deformationType,
                                                     const std::string& identifier )
{
    std::string idString = ( identifier != "" ) ? ( ", tide raised by " + identifier ) : "";
    return idString + ", deformation type " + std::to_string( deformationType );
}

//! Function to get the dependent variable ID sub-string for a partial w.r.t. a body's translational state
static std::string getStatePartialIdString( const std::string& derivativeWrtBody )
{
    return " w.r.t. body translational state of " + derivativeWrtBody;
}

//! Function to get the dependent variable ID sub-string for a partial w.r.t. an estimatable parameter
static std::string getParameterPartialIdString(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > parameterSettings )
{
    if( parameterSettings == nullptr )
    {
        throw std::runtime_error( "Error when getting dependent variable ID, no parameter settings found for parameter partial" );
    }
    std::string idString = " w.r.t. parameter " +
            estimatable_parameters::getParameterTypeString( parameterSettings->parameterType_.first ) + " of " +
            parameterSettings->parameterType_.second.first;
    if( parameterSettings->parameterType_.second.second != "" )
    {
        idString += ", " + parameterSettings->parameterType_.second.second;
    }
    return idString;
}

//! Function to get the dependent variable ID sub-string for a body and the radiation source w.r.t. which it is defined
static std::string getBodyAndSourceIdString( const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings )
{
    return " of body " + dependentVariableSettings->associatedBody_ + " w.r.t. source " + dependentVariableSettings->secondaryBody_;
}

//! Function to get the dependent variable ID sub-string for the model type (acceleration, torque, frame, angle)
static std::string getModelTypeIdString( const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings )
{
    switch( dependentVariableSettings->dependentVariableType_ )
    {
        case single_acceleration_dependent_variable:
        case single_acceleration_norm_dependent_variable:
            return basic_astrodynamics::getAccelerationModelName(
                    getDerivedSettingsForId< SingleAccelerationDependentVariableSaveSettings >( dependentVariableSettings )
                            ->accelerationModelType_ );
        case single_torque_dependent_variable:
        case single_torque_norm_dependent_variable:
            return basic_astrodynamics::getTorqueModelName(
                    getDerivedSettingsForId< SingleTorqueDependentVariableSaveSettings >( dependentVariableSettings )->torqueModelType_ );
        case intermediate_aerodynamic_rotation_matrix_variable: {
            std::shared_ptr< IntermediateAerodynamicRotationVariableSaveSettings > rotationSettings =
                    getDerivedSettingsForId< IntermediateAerodynamicRotationVariableSaveSettings >( dependentVariableSettings );
            return reference_frames::getAerodynamicFrameName( rotationSettings->baseFrame_ ) + "to " +
                    reference_frames::getAerodynamicFrameName( rotationSettings->targetFrame_ );
        }
        case relative_body_aerodynamic_orientation_angle_variable:
            return reference_frames::getAerodynamicAngleName(
                    getDerivedSettingsForId< BodyAerodynamicAngleVariableSaveSettings >( dependentVariableSettings )->angle_ );
        default:
            return "";
    }
}

//! Function to get the dependent variable ID sub-string for the associated and secondary body
static std::string getBodiesIdString( const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings )
{
    const std::string& associatedBody = dependentVariableSettings->associatedBody_;
    const std::string& secondaryBody = dependentVariableSettings->secondaryBody_;
    switch( dependentVariableSettings->dependentVariableType_ )
    {
        case single_acceleration_dependent_variable:
        case single_acceleration_norm_dependent_variable:
        case spherical_harmonic_acceleration_terms_dependent_variable:
        case spherical_harmonic_acceleration_norm_terms_dependent_variable:
        case acceleration_partial_wrt_body_translational_state:
            return ", acting on " + associatedBody + ( ( secondaryBody != associatedBody ) ? ( ", exerted by " + secondaryBody ) : "" );
        case visible_and_emitting_source_panel_count:
        case visible_source_area:
            return " of " + secondaryBody + " as seen from " + associatedBody;
        case received_irradiance:
        case received_fraction:
            return " at " + associatedBody + " due to " + secondaryBody;
        default:
            return "of " + associatedBody + ( ( secondaryBody != "" ) ? ( " w.r.t. " + secondaryBody ) : "" );
    }
}

//! Function to get the dependent variable ID sub-string for the settings that are specific to a derived settings class. Every
//! member that distinguishes two dependent variables must be included, since createDependentVariableListFunction throws an error
//! for different settings with identical IDs.
static std::string getDerivedSettingsIdString( const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings )
{
    switch( dependentVariableSettings->dependentVariableType_ )
    {
        case spherical_harmonic_acceleration_terms_dependent_variable:
        case spherical_harmonic_acceleration_norm_terms_dependent_variable:
            return getDegreeOrderTermsIdString(
                    getDerivedSettingsForId< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >( dependentVariableSettings )
                            ->componentIndices_ );
        case total_spherical_harmonic_cosine_coefficient_variation:
        case total_spherical_harmonic_sine_coefficient_variation:
            return getDegreeOrderTermsIdString(
                    getDerivedSettingsForId< TotalGravityFieldVariationSettings >( dependentVariableSettings )->componentIndices_ );
        case single_gravity_field_variation_acceleration: {
            std::shared_ptr< SingleVariationSphericalHarmonicAccelerationSaveSettings > variationSettings =
                    getDerivedSettingsForId< SingleVariationSphericalHarmonicAccelerationSaveSettings >( dependentVariableSettings );
            return getGravityFieldVariationIdString( variationSettings->deformationType_, variationSettings->identifier_ );
        }
        case single_gravity_field_variation_acceleration_terms: {
            std::shared_ptr< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings > termsSettings =
                    getDerivedSettingsForId< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
                            dependentVariableSettings );
            return getDegreeOrderTermsIdString( termsSettings->componentIndices_ ) +
                    getGravityFieldVariationIdString( termsSettings->deformationType_, termsSettings->identifier_ );
        }
        case acceleration_partial_wrt_body_translational_state:
            return getStatePartialIdString(
                    getDerivedSettingsForId< AccelerationPartialWrtStateSaveSettings >( dependentVariableSettings )->derivativeWrtBody_ );
        case total_acceleration_partial_wrt_body_translational_state:
            return getStatePartialIdString(
                    getDerivedSettingsForId< TotalAccelerationPartialWrtStateSaveSettings >( dependentVariableSettings )
                            ->derivativeWrtBody_ );
        case acceleration_derivative_partial_wrt_parameter: {
            std::shared_ptr< AccelerationDerivativePartialWrtParameterSaveSettings > partialSettings =
                    getDerivedSettingsForId< AccelerationDerivativePartialWrtParameterSaveSettings >( dependentVariableSettings );
            return getParameterPartialIdString( partialSettings->parameterSettings_ ) + ", acceleration type " +
                    basic_astrodynamics::getAccelerationModelName( partialSettings->accelerationModelType_ );
        }
        case total_acceleration_derivative_partial_wrt_parameter:
            return getParameterPartialIdString(
                    getDerivedSettingsForId< TotalAccelerationDerivativePartialWrtParameterSaveSettings >( dependentVariableSettings )
                            ->parameterSettings_ );
        case local_wind_velocity_dependent_variable:
            return ", in " +
                    reference_frames::getAerodynamicFrameName(
                            getDerivedSettingsForId< LocalWindVelocityDependentVariableSaveSettings >( dependentVariableSettings )
                                    ->targetFrame_ );
        case aerodynamic_control_surface_force_coefficients_increment_dependent_variable:
        case aerodynamic_control_surface_moment_coefficients_increment_dependent_variable:
            return ", control surface " +
                    getDerivedSettingsForId< ControlSurfaceCoefficientDependentVariableSettings >( dependentVariableSettings )
                            ->controlSurfaceName_;
        case minimum_constellation_distance:
            return getBodyListIdString(
                    getDerivedSettingsForId< MinimumConstellationDistanceDependentVariableSaveSettings >( dependentVariableSettings )
                            ->bodiesToCheck_ );
        case minimum_constellation_ground_station_distance: {
            std::shared_ptr< MinimumConstellationStationDistanceDependentVariableSaveSettings > stationDistanceSettings =
                    getDerivedSettingsForId< MinimumConstellationStationDistanceDependentVariableSaveSettings >(
                            dependentVariableSettings );
            return getBodyListIdString( stationDistanceSettings->bodiesToCheck_ ) + ", elevation limit " +
                    std::to_string( stationDistanceSettings->elevationAngleLimit_ );
        }
        case custom_dependent_variable: {
            // Custom functions cannot be compared, so the settings object itself identifies the variable
            std::ostringstream settingsAddress;
            settingsAddress << dependentVariableSettings.get( );
            return "(settings object " + settingsAddress.str( ) + ")";
        }
        case illuminated_panel_fraction:
            return getBodyAndSourceIdString( dependentVariableSettings ) + " for panel " +
                    getDerivedSettingsForId< IlluminatedPanelFractionDependentVariableSaveSettings >( dependentVariableSettings )
                            ->panelTypeId_;
        case cross_section_change:
        case actual_cross_section:
            return getBodyAndSourceIdString( dependentVariableSettings ) + ", acceleration type " +
                    getDerivedSettingsForId< CrossSectionDependentVariableSaveSettings >( dependentVariableSettings )->accelerationType_;
        case full_body_paneled_geometry:
        case aerodynamic_coefficients:
            return " of body " + dependentVariableSettings->associatedBody_;
        default:
            return "";
    }
}

//! Function to get a string representing a 'named identification' of a dependent variable
std::string getDependentVariableId( const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings )
{
    std::string variableId = getDependentVariableName( dependentVariableSettings ) + getModelTypeIdString( dependentVariableSettings ) +
            getBodiesIdString( dependentVariableSettings ) + getDerivedSettingsIdString( dependentVariableSettings );
    if( dependentVariableSettings->componentIndex_ >= 0 )
    {
        variableId += ", component " + std::to_string( dependentVariableSettings->componentIndex_ );
    }
    return variableId;
}

}  // namespace propagators

}  // namespace tudat
