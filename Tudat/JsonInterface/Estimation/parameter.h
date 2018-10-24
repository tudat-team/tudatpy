/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_PARAMETER_H
#define TUDAT_JSONINTERFACE_PARAMETER_H

#include "Tudat/SimulationSetup/EstimationSetup/estimatableParameterSettings.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace estimatable_parameters
{


// PropagationDependentVariables

//! Map of `PropagationDependentVariables` string representations.
static std::map< EstimatebleParametersEnum, std::string > estimatedParameterTypes =
{
    { arc_wise_initial_body_state, "arcWiseInitialBodyState" },
    { initial_body_state, "initialBodyState" },
    { gravitational_parameter, "gravitationalParameter" },
    { constant_drag_coefficient, "constantDragCoefficient" },
    { radiation_pressure_coefficient, "radiationPressureCoefficient" },
    { arc_wise_radiation_pressure_coefficient, "arcWiseRadiationPressureCoefficient" },
    { spherical_harmonics_cosine_coefficient_block, "sphericalHarmonicsCosineCoefficientBlock" },
    { spherical_harmonics_sine_coefficient_block, "sphericalHarmonicsSineCoefficientBlock" },
    { constant_rotation_rate, "constantRotationRate" },
    { rotation_pole_position, "rotationPolePosition" },
    { constant_additive_observation_bias, "constantAdditiveObservationBias" },
    { arcwise_constant_additive_observation_bias, "arcwiseConstantAdditiveObservationBias" },
    { constant_relative_observation_bias, "constantRelativeObservationBias" },
    { arcwise_constant_relative_observation_bias, "arcwiseConstantRelativeObservationBias" },
    { ppn_parameter_gamma, "ppnParameterGamma" },
    { ppn_parameter_beta, "ppnParameterBeta" },
    { ground_station_position, "groundStationPosition" },
    { equivalence_principle_lpi_violation_parameter, "equivalencePrincipleLpiViolationParameter" },
    { empirical_acceleration_coefficients, "empiricalAccelerationCoefficients" },
    { arc_wise_empirical_acceleration_coefficients, "arcWiseEmpiricalAccelerationCoefficients" },
    { full_degree_tidal_love_number, "fullDegreeTidalLoveNumber" },
    { single_degree_variable_tidal_love_number, "singleDegreeVariableTidalLoveNumber" },
    { direct_dissipation_tidal_time_lag, "directDissipationTidalTimeLag" }
};


//! Convert `PropagationDependentVariables` to `json`.
inline void to_json( nlohmann::json& jsonObject, const EstimatebleParametersEnum& parameterType )
{
    jsonObject = json_interface::stringFromEnum( parameterType, estimatedParameterTypes );
}

//! Convert `json` to `PropagationDependentVariables`.
inline void from_json( const nlohmann::json& jsonObject, EstimatebleParametersEnum& parameterType )
{
    parameterType = json_interface::enumFromString( jsonObject, estimatedParameterTypes );
}


//! Create a `json` object from a shared pointer to a `EstimatableParameterSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< EstimatableParameterSettings >& parameterSettings );

//! Create a shared pointer to a `EstimatableParameterSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< EstimatableParameterSettings >& parameterSettings );


} // namespace propagators

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_PARAMETER_H
