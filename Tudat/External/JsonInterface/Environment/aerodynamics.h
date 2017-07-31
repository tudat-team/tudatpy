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

//! Map of `AerodynamicCoefficientTypes` string representations.
static std::map< AerodynamicCoefficientTypes, std::string > aerodynamicCoefficientTypes =
{
    { constant_aerodynamic_coefficients, "constantAerodynamicCoefficients" },
    { hypersonic_local_inclincation_coefficients, "hypersonicLocalInclincationCoefficients" },
    { tabulated_coefficients, "tabulatedCoefficients" }
};

//! `AerodynamicCoefficientTypes` not supported by `json_interface`.
static std::vector< AerodynamicCoefficientTypes > unsupportedAerodynamicCoefficientTypes = { };

//! Convert `AerodynamicCoefficientTypes` to `json`.
inline void to_json( json& jsonObject, const AerodynamicCoefficientTypes& aerodynamicCoefficientType )
{
    jsonObject = json_interface::stringFromEnum( aerodynamicCoefficientType, aerodynamicCoefficientTypes );
}

//! Convert `json` to `AerodynamicCoefficientTypes`.
inline void from_json( const json& jsonObject, AerodynamicCoefficientTypes& aerodynamicCoefficientType )
{
    aerodynamicCoefficientType =
            json_interface::enumFromString( jsonObject.get< std::string >( ), aerodynamicCoefficientTypes );
}

//! Create a `json` object from a shared pointer to a `AerodynamicCoefficientSettings` object.
void to_json( json& jsonObject,
              const boost::shared_ptr< AerodynamicCoefficientSettings >& aerodynamicCoefficientSettings );

//! Create a `json` object from a shared pointer to a `AerodynamicCoefficientSettings` object.
void from_json( const json& jsonObject,
                boost::shared_ptr< AerodynamicCoefficientSettings >& aerodynamicCoefficientSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_AERODYNAMICS_H
