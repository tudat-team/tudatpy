/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_DEPENDENTVARIABLE_H
#define TUDAT_JSONINTERFACE_DEPENDENTVARIABLE_H

#include <Tudat/SimulationSetup/PropagationSetup/propagationOutputSettings.h>

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace propagators
{

//! Map of `PropagationDependentVariables` string representations.
static std::map< PropagationDependentVariables, std::string > dependentVariables =
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
    { single_acceleration_norm_dependent_variable, "singleAccelerationNorm" },
    { total_acceleration_dependent_variable, "totalAcceleration" },
    { single_acceleration_dependent_variable, "singleAcceleration" },
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
    { single_torque_norm_dependent_variable, "singleTorqueNorm" },
    { total_torque_dependent_variable, "totalTorque" },
    { single_torque_dependent_variable, "singleTorque" },
    { body_fixed_groundspeed_based_velocity_variable, "bodyFixedGroundspeedBasedVelocity" }
};

// //! `PropagationDependentVariables` not supported by `json_interface`.
// static std::vector< PropagationDependentVariables > unsupportedDependentVariables = { };

//! Convert `PropagationDependentVariables` to `json`.
inline void to_json( json& jsonObject, const PropagationDependentVariables& dependentVariable )
{
    jsonObject = json_interface::stringFromEnum( dependentVariable, dependentVariables );
}

//! Convert `json` to `PropagationDependentVariables`.
inline void from_json( const json& jsonObject, PropagationDependentVariables& dependentVariable )
{
    dependentVariable = json_interface::enumFromString( jsonObject.get< std::string >( ), dependentVariables );
}

//! Create a `json` object from a shared pointer to a `SingleDependentVariableSaveSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< SingleDependentVariableSaveSettings >& variableSettings );

//! Create a shared pointer to a `SingleDependentVariableSaveSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< SingleDependentVariableSaveSettings >& variableSettings );

} // namespace propagators

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_DEPENDENTVARIABLE_H
