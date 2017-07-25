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

#include "propagator.h"

namespace tudat
{

namespace propagators
{

//! Convert `IntegratedStateType` to `json`.
void to_json( json& jsonObject, const IntegratedStateType& integratedStateType )
{
    jsonObject = json_interface::stringFromEnum( integratedStateType, integratedStateTypes );
}

//! Convert `json` to `IntegratedStateType`.
void from_json( const json& jsonObject, IntegratedStateType& integratedStateType )
{
    integratedStateType = json_interface::enumFromString( jsonObject.get< std::string >( ), integratedStateTypes );
}

} // namespace propagators

} // namespace tudat
