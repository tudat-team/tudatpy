/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_SAVE_H
#define TUDAT_JSONINTERFACE_SAVE_H

#include <Tudat/SimulationSetup/PropagationSetup/propagationOutputSettings.h>

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace propagators
{

//! Create a `json` object from a shared pointer to a `DependentVariableSaveSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< DependentVariableSaveSettings >& saveSettings );

//! Create a shared pointer to a `DependentVariableSaveSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< DependentVariableSaveSettings >& saveSettings );

} // namespace propagators

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_SAVE_H
