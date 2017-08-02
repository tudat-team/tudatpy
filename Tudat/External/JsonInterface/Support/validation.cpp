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

#include "validation.h"

namespace tudat
{

namespace json_interface
{

//! Create a `json` object from a shared pointer to a `ValidationSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< ValidationSettings >& validationSettings )
{
    if ( ! validationSettings )
    {
    return;
    }
    using K = Keys::Validation;

    jsonObject[ K::usingDefaultValueForMissingKey ] = validationSettings->usingDefaultValueForMissingKey;
    jsonObject[ K::unusedKey ] = validationSettings->unusedKey;
}

//! Create a shared pointer to a `ValidationSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< ValidationSettings >& validationSettings )
{
    using K = Keys::Validation;

    validationSettings = boost::make_shared< ValidationSettings >( );
    updateFromJSONIfDefined( validationSettings->usingDefaultValueForMissingKey,
                             jsonObject, K::usingDefaultValueForMissingKey );
    updateFromJSONIfDefined( validationSettings->unusedKey, jsonObject, K::unusedKey );
}

} // namespace json_interface

} // namespace tudat
