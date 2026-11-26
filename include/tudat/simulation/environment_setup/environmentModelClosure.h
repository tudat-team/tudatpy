/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ENVIRONMENTMODELCLOSURE_H
#define TUDAT_ENVIRONMENTMODELCLOSURE_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to check whether a model is of a given type (or of a type derived from it)
template< typename ModelType, typename BaseType >
bool isModelOfType( const std::shared_ptr< BaseType > model )
{
    return std::dynamic_pointer_cast< ModelType >( model ) != nullptr;
}

/*
 *  Closure of environment models: links of an environment model to environment models of other bodies (or of the same body). Models
 *  are created without closure, after which the closure is applied by closeSystemOfBodies, or, for models added
 *  later, by the function adding them. Which environment models require closure is defined (in a single place) in
 *  environmentModelClosure.cpp; the closure itself is applied by the applyEnvironmentModelClosure functions of each type of
 *  environment model (in the setup layer), which use the applyEnvironmentModelClosure member function of the model.
 */

//! Types of environment models of a body that may require closure
enum EnvironmentModelClosureType {
    rotation_model_closure,
    gravity_field_variations_closure,
    flight_conditions_closure,
    aerodynamic_coefficients_closure,
    vehicle_systems_closure,
    radiation_source_closure,
    radiation_pressure_target_closure,
    atmosphere_closure,
    ground_stations_closure,
    body_deformation_closure
};

//! Function to apply the closure of the environment models of a body (of the given types, or of all types if none are given) that
//! require closure. Closures that require other closures are applied after these. Closures that require flight conditions that do
//! not yet exist are not applied; they are applied when the flight conditions are added (see addFlightConditions)
void closeEnvironmentModels(
        const std::string& bodyName,
        const SystemOfBodies& bodies,
        const std::vector< EnvironmentModelClosureType >& closureTypes = std::vector< EnvironmentModelClosureType >( ) );

//! Function to apply the closure of the environment models of all bodies (see closeEnvironmentModels)
void closeSystemOfBodies( const SystemOfBodies& bodies );

//! Function to apply the closure of (newly added) flight conditions of a body, and of the models of the body that require these
//! flight conditions and of which the closure has not yet been applied. Models of which the closure has been applied are not
//! modified.
void closeFlightConditionsAndDependentModels( const std::string& bodyName, const SystemOfBodies& bodies );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_ENVIRONMENTMODELCLOSURE_H
