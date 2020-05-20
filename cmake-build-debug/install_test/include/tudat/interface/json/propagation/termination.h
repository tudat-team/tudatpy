/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_TERMINATION_H
#define TUDAT_JSONINTERFACE_TERMINATION_H

#include "tudat/simulation/propagation/propagationTerminationSettings.h"
#include "tudat/interface/json/support/valueAccess.h"
#include "tudat/interface/json/support/valueConversions.h"

namespace tudat
{

namespace propagators
{

// PropagationHybridTerminationSettings

//! Create a `json` object from a shared pointer to a `PropagationHybridTerminationSettings` object.
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< PropagationHybridTerminationSettings >& hybridTerminationSettings );

//! Create a shared pointer to a `PropagationHybridTerminationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject,
                std::shared_ptr< PropagationHybridTerminationSettings >& hybridTerminationSettings );


// PropagationTerminationSettings

//! Create a `json` object from a shared pointer to a `PropagationTerminationSettings` object.
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< PropagationTerminationSettings >& terminationSettings );

//! Create a shared pointer to a `PropagationTerminationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject,
                std::shared_ptr< PropagationTerminationSettings >& terminationSettings );

} // namespace propagators

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_TERMINATION_H
