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
 *  are created (or copied) without closure, after which the closure is applied by closeSystemOfBodies. Which environment models
 *  require closure is defined (in a single place) in environmentModelClosure.cpp; the closure itself is applied by the
 *  applyEnvironmentModelClosure functions of each type of environment model (in the setup layer), which use the
 *  applyEnvironmentModelClosure member function of the model.
 */

//! Function to register (in the body) which of the environment models of a body require closure, to be called after the models
//! have been created or copied without closure
void registerEnvironmentModelClosures( const std::shared_ptr< Body > body );

//! Function to apply all closures of environment models (registered in the bodies) that have not yet been applied. Closures are
//! applied repeatedly, until no more closures can be applied (since some closures require other closures to be applied first).
//! Returns the closures that could not be applied, per body (for which a warning is printed if printWarnings is true).
std::map< std::string, std::vector< EnvironmentModelClosureType > > closeSystemOfBodies( const SystemOfBodies& bodies,
                                                                                         const bool printWarnings = true );

//! Function to retrieve the closures of environment models (registered in the bodies) that have not been applied, per body
std::map< std::string, std::vector< EnvironmentModelClosureType > > getEnvironmentModelsWithoutClosure( const SystemOfBodies& bodies );

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_ENVIRONMENTMODELCLOSURE_H
