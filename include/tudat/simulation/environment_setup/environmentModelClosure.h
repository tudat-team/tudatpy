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

/*
 *  Closure of environment models: links of an environment model to environment models of other bodies (or of the same body). Models
 *  are created (or copied) without closure, after which the closure is applied by closeSystemOfBodies. The functions below define
 *  (in a single place) which environment models require closure; the closure itself is applied by the applyEnvironmentModelClosure
 *  functions of each type of environment model (in the setup layer), which use the applyEnvironmentModelClosure member function of
 *  the model.
 */

//! Function to determine whether a rotation model requires closure (does not check whether it has been applied)
bool doesEnvironmentModelRequireClosure( const std::shared_ptr< ephemerides::RotationalEphemeris > rotationModel );

//! Function to determine whether a set of gravity field variations requires closure (does not check whether it has been applied)
bool doesEnvironmentModelRequireClosure( const std::shared_ptr< gravitation::GravityFieldVariationsSet > variationSet );

//! Function to determine whether flight conditions require closure (does not check whether it has been applied)
bool doesEnvironmentModelRequireClosure( const std::shared_ptr< aerodynamics::FlightConditions > flightConditions );

//! Function to determine whether an aerodynamic coefficient interface requires closure (does not check whether it has been applied)
bool doesEnvironmentModelRequireClosure( const std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface );

//! Function to determine whether vehicle systems require closure (does not check whether it has been applied)
bool doesEnvironmentModelRequireClosure( const std::shared_ptr< system_models::VehicleSystems > vehicleSystems );

//! Function to determine whether a radiation source model requires closure (does not check whether it has been applied)
bool doesEnvironmentModelRequireClosure( const std::shared_ptr< electromagnetism::RadiationSourceModel > radiationSourceModel );

//! Function to determine whether a radiation pressure target model requires closure (does not check whether it has been applied)
bool doesEnvironmentModelRequireClosure( const std::shared_ptr< electromagnetism::RadiationPressureTargetModel > targetModel );

//! Function to determine whether an atmosphere model (including its wind model) requires closure (does not check whether it has
//! been applied)
bool doesEnvironmentModelRequireClosure( const std::shared_ptr< aerodynamics::AtmosphereModel > atmosphereModel );

//! Function to determine whether a ground station requires closure (does not check whether it has been applied)
bool doesEnvironmentModelRequireClosure( const std::shared_ptr< ground_stations::GroundStation > groundStation );

//! Function to determine whether a body deformation model requires closure (does not check whether it has been applied)
bool doesEnvironmentModelRequireClosure( const std::shared_ptr< basic_astrodynamics::BodyDeformationModel > deformationModel );

//! Function to register (in the body) which of the environment models of a body require closure, to be called after the models
//! have been created or copied without closure
void registerEnvironmentModelClosures( const std::shared_ptr< Body > body );

//! Function to apply the closure of a single type of environment model of a body. Returns false if the closure cannot yet be
//! applied (environment models it depends on are not yet closed).
bool applyEnvironmentModelClosure( const EnvironmentModelClosureType closureType,
                                   const std::string& bodyName,
                                   const SystemOfBodies& bodies );

//! Function to apply all closures of environment models (registered in the bodies) that have not yet been applied. Closures are
//! applied repeatedly, until no more closures can be applied (since some closures require other closures to be applied first).
//! Returns the closures that could not be applied, per body.
std::map< std::string, std::vector< EnvironmentModelClosureType > > closeSystemOfBodies( const SystemOfBodies& bodies );

//! Function to retrieve the closures of environment models (registered in the bodies) that have not been applied, per body
std::map< std::string, std::vector< EnvironmentModelClosureType > > getEnvironmentModelsWithoutClosure( const SystemOfBodies& bodies );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_ENVIRONMENTMODELCLOSURE_H
