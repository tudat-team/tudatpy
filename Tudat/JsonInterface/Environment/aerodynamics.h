/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_AERODYNAMICS_H
#define TUDAT_JSONINTERFACE_AERODYNAMICS_H

#include "Tudat/SimulationSetup/EnvironmentSetup/createAerodynamicCoefficientInterface.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace aerodynamics
{

// AerodynamicCoefficientsIndependentVariables

//! Map of `AerodynamicCoefficientsIndependentVariables` string representations.
static std::map< AerodynamicCoefficientsIndependentVariables, std::string > aerodynamicVariables =
{
    { mach_number_dependent, "machNumber" },
    { angle_of_attack_dependent, "angleOfAttack" },
    { angle_of_sideslip_dependent, "angleOfSideslip" },
    { altitude_dependent, "altitude" },
    { control_surface_deflection_dependent, "controlSurfaceDeflection" },
    { undefined_independent_variable, "undefined" }
};

//! `AerodynamicCoefficientsIndependentVariables` not supported by `json_interface`.
static std::vector< AerodynamicCoefficientsIndependentVariables > unsupportedAerodynamicVariables = { };

//! Convert `AerodynamicCoefficientsIndependentVariables` to `json`.
inline void to_json( nlohmann::json& jsonObject, const AerodynamicCoefficientsIndependentVariables& aerodynamicVariable )
{
    jsonObject = json_interface::stringFromEnum( aerodynamicVariable, aerodynamicVariables );
}

//! Convert `json` to `AerodynamicCoefficientsIndependentVariables`.
inline void from_json( const nlohmann::json& jsonObject, AerodynamicCoefficientsIndependentVariables& aerodynamicVariable )
{
    aerodynamicVariable = json_interface::enumFromString( jsonObject, aerodynamicVariables );
}

}  // namespace aerodynamics


namespace simulation_setup
{

// AerodynamicCoefficientTypes

//! Map of `AerodynamicCoefficientTypes` string representations.
static std::map< AerodynamicCoefficientTypes, std::string > aerodynamicCoefficientTypes =
{
    { constant_aerodynamic_coefficients, "constant" },
    { hypersonic_local_inclincation_coefficients, "hypersonicLocalInclincation" },
    { tabulated_coefficients, "tabulated" }
};

//! `AerodynamicCoefficientTypes` not supported by `json_interface`.
static std::vector< AerodynamicCoefficientTypes > unsupportedAerodynamicCoefficientTypes =
{
    hypersonic_local_inclincation_coefficients
};

//! Convert `AerodynamicCoefficientTypes` to `json`.
inline void to_json( nlohmann::json& jsonObject, const AerodynamicCoefficientTypes& aerodynamicCoefficientType )
{
    jsonObject = json_interface::stringFromEnum( aerodynamicCoefficientType, aerodynamicCoefficientTypes );
}

//! Convert `json` to `AerodynamicCoefficientTypes`.
inline void from_json( const nlohmann::json& jsonObject, AerodynamicCoefficientTypes& aerodynamicCoefficientType )
{
    aerodynamicCoefficientType = json_interface::enumFromString( jsonObject, aerodynamicCoefficientTypes );
}


// AerodynamicCoefficientSettings

//! Create a `json` object from a shared pointer to a `AerodynamicCoefficientSettings` object.
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< AerodynamicCoefficientSettings >& aerodynamicSettings );

//! Create a `json` object from a shared pointer to a `AerodynamicCoefficientSettings` object.
void from_json(const nlohmann::json& jsonObject, std::shared_ptr< AerodynamicCoefficientSettings >& aerodynamicSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_AERODYNAMICS_H
