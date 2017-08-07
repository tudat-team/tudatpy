/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_VARIABLE_H
#define TUDAT_JSONINTERFACE_VARIABLE_H

#include <Tudat/SimulationSetup/PropagationSetup/propagationOutputSettings.h>

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace propagators
{

// VariableType

//! Map of `VariableType`s string representations.
static std::map< VariableType, std::string > variableTypes =
{
    { epochVariable, "epoch" },
    { stateVariable, "state" },
    { dependentVariable, "dependent" }
};

//! `VariableType`s not supported by `json_interface`.
static std::vector< VariableType > unsupportedVariableTypes = { };

//! Convert `VariableType` to `json`.
inline void to_json( json& jsonObject, const VariableType& variableType )
{
    jsonObject = json_interface::stringFromEnum( variableType, variableTypes );
}

//! Convert `json` to `VariableType`.
inline void from_json( const json& jsonObject, VariableType& variableType )
{
    variableType = json_interface::enumFromString( jsonObject.get< std::string >( ), variableTypes );
}


// PropagationDependentVariables

//! Map of `PropagationDependentVariables` string representations.
static std::map< PropagationDependentVariables, std::string > dependentVariableTypes =
{
    { mach_number_dependent_variable, "machNumber" },
    { altitude_dependent_variable, "altitude" },
    { airspeed_dependent_variable, "airspeed" },
    { local_density_dependent_variable, "localDensity" },
    { relative_speed_dependent_variable, "relativeSpeed" },
    { relative_position_dependent_variable, "relativePosition" },
    { relative_distance_dependent_variable, "relativeDistance" },
    { relative_velocity_dependent_variable, "relativeVelocity" },
    { radiation_pressure_dependent_variable, "radiationPressure" },
    { total_acceleration_norm_dependent_variable, "totalAccelerationNorm" },
    { single_acceleration_norm_dependent_variable, "accelerationNorm" },
    { total_acceleration_dependent_variable, "totalAcceleration" },
    { single_acceleration_dependent_variable, "acceleration" },
    { aerodynamic_force_coefficients_dependent_variable, "aerodynamicForceCoefficients" },
    { aerodynamic_moment_coefficients_dependent_variable, "aerodynamicMomentCoefficients" },
    { rotation_matrix_to_body_fixed_frame_variable, "rotationMatrixToBodyFixedFrame" },
    { intermediate_aerodynamic_rotation_matrix_variable, "intermediateAerodynamicRotationMatrix" },
    { relative_body_aerodynamic_orientation_angle_variable, "relativeBodyAerodynamicOrientationAngle" },
    { body_fixed_airspeed_based_velocity_variable, "bodyFixedAirspeedBasedVelocity" },
    { total_aerodynamic_g_load_variable, "totalAerodynamicGLoad" },
    { stagnation_point_heat_flux_dependent_variable, "stagnationPointHeatFlux" },
    { local_temperature_dependent_variable, "localTemperature" },
    { geodetic_latitude_dependent_variable, "geodeticLatitude" },
    { control_surface_deflection_dependent_variable, "controlSurfaceDeflection" },
    { total_mass_rate_dependent_variables, "totalMassRates" },
    { lvlh_to_inertial_frame_rotation_dependent_variable, "lvlhToInertialFrameRotation" },
    { periapsis_altitude_dependent_variable, "periapsisAltitude" },
    { total_torque_norm_dependent_variable, "totalTorqueNorm" },
    { single_torque_norm_dependent_variable, "torqueNorm" },
    { total_torque_dependent_variable, "totalTorque" },
    { single_torque_dependent_variable, "torque" },
    { body_fixed_groundspeed_based_velocity_variable, "bodyFixedGroundspeedBasedVelocity" }
};

// //! `PropagationDependentVariables` not supported by `json_interface`.
// static std::vector< PropagationDependentVariables > unsupportedDependentVariables = { };

//! Convert `PropagationDependentVariables` to `json`.
inline void to_json( json& jsonObject, const PropagationDependentVariables& dependentVariable )
{
    jsonObject = json_interface::stringFromEnum( dependentVariable, dependentVariableTypes );
}

//! Convert `json` to `PropagationDependentVariables`.
inline void from_json( const json& jsonObject, PropagationDependentVariables& dependentVariable )
{
    dependentVariable = json_interface::enumFromString( jsonObject.get< std::string >( ), dependentVariableTypes );
}


// VariableSettings

//! Create a `json` object from a shared pointer to a `VariableSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< VariableSettings >& variableSettings );

//! Create a shared pointer to a `VariableSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< VariableSettings >& variableSettings );


// SingleDependentVariableSaveSettings

//! Create a `json` object from a shared pointer to a `SingleDependentVariableSaveSettings` object.
void to_json( json& jsonObject,
              const boost::shared_ptr< SingleDependentVariableSaveSettings >& dependentVariableSettings );

//! Create a shared pointer to a `SingleDependentVariableSaveSettings` object from a `json` object.
void from_json( const json& jsonObject,
                boost::shared_ptr< SingleDependentVariableSaveSettings >& dependentVariableSettings );

} // namespace propagators


namespace json_interface
{

//! Get a shared pointer to a `VariableSettings` object created from a `jsonObject[ keyPath ]`.
/*!
 * Get a shared pointer to a `VariableSettings` object created from `jsonObject[ keyPath ]`.
 *
 * If `jsonObject[ keyPath ]` is an object, the variable is created directly using the `from_json` function.
 * If `jsonObject[ keyPath ]` is a string, depending on its value the following can happen:
 * - "epoch": the variable is created automatically as a `VariableSettings` with type `epochVariable`.
 * - "bodyName.state": the variable is created automatically as a `BodyVariableSettings` with type `stateVariable`
 * and associated to body named "bodyName".
 * - "customName": the variable will be retrieved (if defined) from the root `json` object at "variables.customName".
 *
 * \param jsonObject `json` object containing, at `keyPath`, the variable settings or the variable name and a root
 * `json` object.
 * \param keyPath The key path at which the variable settings are specified in `jsonObject`.
 * \return Shared pointer to a `VariableSettings` object created from `jsonObject[ keyPath ]`.
 */
boost::shared_ptr< propagators::VariableSettings > getVariable( const json& jsonObject, const KeyPath& keyPath );

//! Get a vector of shared pointers to `VariableSettings` objects created from `jsonObject[ keyPath ]`.
/*!
 * Get a vector of shared pointers to `VariableSettings` objects created from `jsonObject[ keyPath ]`.
 * Each of the variable settings of the vector are created using the `getVariable` function.
 * \param jsonObject `json` object containing, at `keyPath`, the variables settings or the variable names and a root
 * `json` object.
 * \param keyPath The key path at which the variables settings are specified in `jsonObject`.
 * \return Vector of shared pointers to `VariableSettings` objects created from `jsonObject[ keyPath ]`.
 */
std::vector< boost::shared_ptr< propagators::VariableSettings > > getVariables(
        const json& jsonObject, const KeyPath& keyPath );

//! Get a vector of shared pointers to `SingleDependentVariableSaveSettings` objects created from
//! `jsonObject[ keyPath ]`.
/*!
 * Get a vector of shared pointers to `SingleDependentVariableSaveSettings` objects created from
 * `jsonObject[ keyPath ]`.
 * Each of the variable settings of the vector are created using the `getVariable` function. Only the variables
 * defined in `jsonObject[ keyPath ]` that are of type `dependentVariable` are included in the vector.
 * \param jsonObject `json` object containing, at `keyPath`, the dependent variables settings or the dependent variable
 * names and a root `json` object.
 * \param keyPath The key path at which the dependent variables settings are specified in `jsonObject`.
 * \return Vector of shared pointers to `SingleDependentVariableSaveSettings` objects created from
 * `jsonObject[ keyPath ]`.
 */
std::vector< boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > > getDependentVariables(
        const json& jsonObject, const KeyPath& keyPath );

} // json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_VARIABLE_H
