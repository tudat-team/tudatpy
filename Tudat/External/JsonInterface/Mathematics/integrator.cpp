/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "integrator.h"

namespace tudat
{

namespace numerical_integrators
{

//! Convert `AvailableIntegrators` to `json`.
void to_json( json& jsonObject, const AvailableIntegrators& availableIntegrator )
{
    jsonObject = json_interface::stringFromEnum( availableIntegrator, integratorTypes );
}

//! Convert `json` to `AvailableIntegrators`.
void from_json( const json& jsonObject, AvailableIntegrators& availableIntegrator )
{
    availableIntegrator = json_interface::enumFromString( jsonObject.get< std::string >( ), integratorTypes );
}

//! Convert `RungeKuttaCoefficients::CoefficientSets` to `json`.
void to_json( json& jsonObject, const RungeKuttaCoefficients::CoefficientSets& rungeKuttaCoefficientSet )
{
    jsonObject = json_interface::stringFromEnum( rungeKuttaCoefficientSet, rungeKuttaCoefficientSets );
}

//! Convert `json` to `RungeKuttaCoefficients::CoefficientSets`.
void from_json( const json& jsonObject, RungeKuttaCoefficients::CoefficientSets& rungeKuttaCoefficientSet )
{
    rungeKuttaCoefficientSet = json_interface::enumFromString(
		jsonObject.get< std::string >( ), rungeKuttaCoefficientSets );
}

} // namespace numerical_integrators

} // namespace tudat
