/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

std::string getParameterTypeString( const EstimatebleParametersEnum parameterType )
{
    std::string parameterDescription;
    switch( parameterType )
    {
    case arc_wise_initial_body_state:
        parameterDescription = "arc-wise translational state ";
        break;
    case initial_body_state:
        parameterDescription = "translational state ";
        break;
    case initial_rotational_body_state:
        parameterDescription = "rotational state ";
        break;
    case initial_mass_state:
        parameterDescription = "rotational state ";
        break;
    case gravitational_parameter:
        parameterDescription = "gravitational parameter ";
        break;
    case constant_drag_coefficient:
        parameterDescription = "constant drag coefficient ";
        break;
    case radiation_pressure_coefficient:
        parameterDescription = "radiation pressure coefficient ";
        break;
    case arc_wise_radiation_pressure_coefficient:
        parameterDescription = "arc-wise radiation pressure coefficient ";
        break;
    case arc_wise_constant_drag_coefficient:
        parameterDescription = "arc-wise drag coefficient ";
        break;
    case spherical_harmonics_cosine_coefficient_block:
        parameterDescription = "cosine spherical harmonic coefficient block ";
        break;
    case spherical_harmonics_sine_coefficient_block:
        parameterDescription = "sine spherical harmonic coefficient block ";
        break;
    case constant_rotation_rate:
        parameterDescription = "constant rotation rate ";
        break;
    case rotation_pole_position:
        parameterDescription = "pole position ";
        break;
    case constant_additive_observation_bias:
        parameterDescription = "absolute observation bias ";
        break;
    case arcwise_constant_additive_observation_bias:
        parameterDescription = "absolute arc-wise observation bias ";
        break;
    case constant_relative_observation_bias:
        parameterDescription = "relative observation bias ";
        break;
    case arcwise_constant_relative_observation_bias:
        parameterDescription = "relative arc-wise observation bias ";
        break;
    case ground_station_position:
        parameterDescription = "ground station position ";
        break;
    case reference_point_position:
        parameterDescription = "reference point position ";
        break;
    case equivalence_principle_lpi_violation_parameter:
        parameterDescription = " equivalence principle violation parameter ";
        break;
    case empirical_acceleration_coefficients:
        parameterDescription = " empirical acceleration coefficients ";
        break;
    case arc_wise_empirical_acceleration_coefficients:
        parameterDescription = " arc-wise empirical acceleration coefficients ";
        break;
    case full_degree_tidal_love_number:
        parameterDescription = " tidal Love number at full degree ";
        break;
    case single_degree_variable_tidal_love_number:
        parameterDescription = " tidal Love number at separate orders of single degree ";
        break;
    case ppn_parameter_gamma:
        parameterDescription = "PPN parameter gamma ";
        break;
    case ppn_parameter_beta:
        parameterDescription = "PPN parameter beta ";
        break;
    case direct_dissipation_tidal_time_lag:
        parameterDescription = " direct tidal time-lag ";
        break;
    case mean_moment_of_inertia:
        parameterDescription = " mean moment of inertia ";
        break;
    case periodic_spin_variation:
        parameterDescription = " periodic spin variation for full planetary rotational model ";
        break;
    case polar_motion_amplitude:
        parameterDescription = " polar motion amplitude for full planetary rotational model";
        break;
    case core_factor:
        parameterDescription = " core factor of the celestial body ";
        break;
    case free_core_nutation_rate:
        parameterDescription = " free core nutation rate of the celestial body";
        break;
    case desaturation_delta_v_values:
        parameterDescription = " momentum wheel desaturation Delta V ";
        break;
    case scaled_longitude_libration_amplitude:
        parameterDescription = " scaled longitude libration amplitude ";
        break;
    case constant_thrust_magnitude_parameter:
        parameterDescription = " constant thrust magnitude ";
        break;
    case constant_specific_impulse:
        parameterDescription = " constant specific impulse ";
        break;
    case constant_time_drift_observation_bias:
        parameterDescription = "time drift observation bias ";
        break;
    case arc_wise_time_drift_observation_bias:
        parameterDescription = "arc-wise time drift observation bias ";
        break;
    case constant_time_observation_bias:
        parameterDescription = "time observation bias ";
        break;
    case arc_wise_time_observation_bias:
        parameterDescription = "arc-wise time observation bias ";
        break;
    case inverse_tidal_quality_factor:
        parameterDescription = " inverse of tidal quality factor ";
        break;
    case yarkovsky_parameter:
        parameterDescription = " Yarkovsky parameter A2 ";
        break;
    case custom_estimated_parameter:
        parameterDescription = " Custom parameter ";
        break;
    case polynomial_gravity_field_variation_amplitudes:
        parameterDescription = " Polynomial gravity field variations ";
        break;
    case periodic_gravity_field_variation_amplitudes:
        parameterDescription = " Periodic gravity field variations ";
        break;
    default:
        std::string errorMessage = "Error when getting parameter string, did not recognize parameter " +
                std::to_string( parameterType );
        throw std::runtime_error( errorMessage );
    }
    return parameterDescription;
}

//! Function to determine whether the given parameter represents an initial dynamical state, or a static parameter.
bool isParameterDynamicalPropertyInitialState( const EstimatebleParametersEnum parameterType )
{
    bool flag;
    switch( parameterType )
    {
    case arc_wise_initial_body_state:
        flag = true;
        break;
    case initial_body_state:
        flag = true;
        break;
    case initial_rotational_body_state:
        flag = true;
        break;
    case initial_mass_state:
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
    case arc_wise_radiation_pressure_coefficient:
        isDoubleParameter = false;
        break;
    case arc_wise_constant_drag_coefficient:
        isDoubleParameter = false;
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
    case arcwise_constant_additive_observation_bias:
        isDoubleParameter = false;
        break;
    case constant_relative_observation_bias:
        isDoubleParameter = false;
        break;
    case arcwise_constant_relative_observation_bias:
        isDoubleParameter = false;
        break;
    case ppn_parameter_gamma:
        isDoubleParameter = true;
        break;
    case ppn_parameter_beta:
        isDoubleParameter = true;
        break;
    case ground_station_position:
        isDoubleParameter = false;
        break;
    case reference_point_position:
        isDoubleParameter = false;
        break;
    case equivalence_principle_lpi_violation_parameter:
        isDoubleParameter = true;
        break;
   case full_degree_tidal_love_number:
        isDoubleParameter = false;
        break;
    case single_degree_variable_tidal_love_number:
         isDoubleParameter = false;
         break;
    case empirical_acceleration_coefficients:
         isDoubleParameter = false;
         break;
     case arc_wise_empirical_acceleration_coefficients:
          isDoubleParameter = false;
          break;
    case direct_dissipation_tidal_time_lag:
         isDoubleParameter = true;
        break;
    case mean_moment_of_inertia:
         isDoubleParameter = true;
        break;
    case desaturation_delta_v_values:
        isDoubleParameter = false;
       break;
    case periodic_spin_variation:
        isDoubleParameter = false;
        break;
    case polar_motion_amplitude:
        isDoubleParameter = false;
        break;
    case core_factor:
        isDoubleParameter = true;
        break;
    case free_core_nutation_rate:
        isDoubleParameter = true;
        break;
    case scaled_longitude_libration_amplitude:
        isDoubleParameter = true;
        break;
    case constant_thrust_magnitude_parameter:
        isDoubleParameter = true;
        break;
    case constant_specific_impulse:
        isDoubleParameter = true;
        break;
    case constant_time_drift_observation_bias:
        isDoubleParameter = false;
        break;
    case arc_wise_time_drift_observation_bias:
        isDoubleParameter = false;
        break;
    case constant_time_observation_bias:
        isDoubleParameter = false;
        break;
    case arc_wise_time_observation_bias:
        isDoubleParameter = false;
        break;
    case inverse_tidal_quality_factor:
        isDoubleParameter = true;
        break;
    case yarkovsky_parameter:
        isDoubleParameter = true;
        break;
    case custom_estimated_parameter:
        isDoubleParameter = false;
        break;
    case polynomial_gravity_field_variation_amplitudes:
        isDoubleParameter = false;
        break;
    case periodic_gravity_field_variation_amplitudes:
        isDoubleParameter = false;
        break;
    default:
        throw std::runtime_error( "Error, parameter type " + std::to_string( parameterType ) +
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
    case initial_rotational_body_state:
        flag = true;
        break;
    case periodic_spin_variation:
        flag = true;
        break;
    case polar_motion_amplitude:
        flag = true;
        break;
    case core_factor:
        flag = true;
        break;
    case free_core_nutation_rate:
        flag = true;
        break;
    case scaled_longitude_libration_amplitude:
        flag = true;
        break;
    default:
        flag = false;
        break;
    }
    return flag;
}

//! Function to determine whether the given parameter influences an observation link directly
bool isParameterObservationLinkProperty( const EstimatebleParametersEnum parameterType )
{
    bool flag;
    switch( parameterType )
    {
    case constant_additive_observation_bias:
        flag = true;
        break;
    case arcwise_constant_additive_observation_bias:
        flag = true;
        break;
    case constant_relative_observation_bias:
        flag = true;
        break;
    case arcwise_constant_relative_observation_bias:
        flag = true;
        break;
    case constant_time_drift_observation_bias:
        flag = true;
        break;
    case arc_wise_time_drift_observation_bias:
        flag = true;
        break;
    case constant_time_observation_bias:
        flag = true;
        break;
    case arc_wise_time_observation_bias:
        flag = true;
        break;
    default:
        flag = false;
        break;
    }
    return flag;
}

//! Function to determine whether the given parameter influences an observation link directly
bool isParameterObservationLinkTimeProperty( const EstimatebleParametersEnum parameterType )
{
    bool flag;
    switch( parameterType )
    {
        case constant_time_observation_bias:
            flag = true;
            break;
        case arc_wise_time_observation_bias:
            flag = true;
            break;
        default:
            flag = false;
            break;
    }
    return flag;
}

//! Function to determine whether the given parameter influences a body's tidal gravity field variations.
bool isParameterTidalProperty( const EstimatebleParametersEnum parameterType )
{
    bool flag;
    switch( parameterType )
    {
    case full_degree_tidal_love_number:
        flag = true;
        break;
    case single_degree_variable_tidal_love_number:
        flag = true;
        break;
    default:
        flag = false;
        break;
    }
    return flag;
}

//! Function to determine whether the given parameter influences a body's tidal gravity field variations.
bool isParameterNonTidalGravityFieldVariationProperty( const EstimatebleParametersEnum parameterType )
{
    bool flag;
    switch( parameterType )
    {
    case polynomial_gravity_field_variation_amplitudes:
        flag = true;
        break;
    case periodic_gravity_field_variation_amplitudes:
        flag = true;
        break;
    default:
        flag = false;
        break;
    }
    return flag;
}

//! Function to determine whether the given parameter represents an arc-wise initial dynamical state.
bool isParameterArcWiseInitialStateProperty( const EstimatebleParametersEnum parameterType )
{
    bool flag;
    switch( parameterType )
    {
        case arc_wise_initial_body_state:
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

