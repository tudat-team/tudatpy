/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <Tudat/SimulationSetup/EnvironmentSetup/createAerodynamicCoefficientInterface.h>

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace simulation_setup
{

//! Map of `AerodynamicCoefficientTypes` supported by `json_interface`.
static std::map< std::string, AerodynamicCoefficientTypes > aerodynamicCoefficientTypes =
{
    { "constantAerodynamicCoefficients",	     constant_aerodynamic_coefficients },
    { "hypersonicLocalInclincationCoefficients", hypersonic_local_inclincation_coefficients },
    { "tabulatedCoefficients",			         tabulated_coefficients }
};

//! Convert `AerodynamicCoefficientTypes` to `json`.
void to_json( json& jsonObject, const AerodynamicCoefficientTypes& aerodynamicCoefficientType );

//! Convert `json` to `AerodynamicCoefficientTypes`.
void from_json( const json& jsonObject, AerodynamicCoefficientTypes& aerodynamicCoefficientType );

//! Create a `json` object from a shared pointer to a `AerodynamicCoefficientSettings` object.
//! Called automatically by `nlohmann::json` when using a constructor such as
//! `json( aerodynamicCoefficientSettings )`.
void to_json( json& jsonObject,
              const boost::shared_ptr< AerodynamicCoefficientSettings >& aerodynamicCoefficientSettings );

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to an `AerodynamicCoefficientSettings` object from a `json` object.
/*!
 * Create a shared pointer to an `AerodynamicCoefficientSettings` object from a `json` object.
 * \param settings `json` object containing the settings for one aerodynamic model.
 * \param keyTree Key tree at which the object containing the aerodynamic coefficients settings can be accessed.
 * Empty if `settings` contains ONLY the aerodynamic coefficients settings.
 * \param fallbackArea Fallback reference area to be used when no reference area is speciefied in `settings`.
 * \return Shared pointer to an `AerodynamicCoefficientSettings` object.
 */
boost::shared_ptr< simulation_setup::AerodynamicCoefficientSettings > createAerodynamicCoefficientSettings(
        const json& settings, const KeyTree& keyTree = { }, const double& fallbackArea = TUDAT_NAN );

} // namespace json_interface


} // namespace tudat

#endif // TUDAT_JSONINTERFACE_AERODYNAMICS_H
