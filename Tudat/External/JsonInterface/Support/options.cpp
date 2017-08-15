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

#include "options.h"

namespace tudat
{

namespace json_interface
{

//! Create a `json` object from a shared pointer to a `ApplicationOptions` object.
void to_json( json& jsonObject, const boost::shared_ptr< ApplicationOptions >& applicationOptions )
{
    if ( ! applicationOptions )
    {
        return;
    }
    using K = Keys::Options;

    jsonObject[ K::defaultValueUsedForMissingKey ] = applicationOptions->defaultValueUsedForMissingKey;
    jsonObject[ K::unusedKey ] = applicationOptions->unusedKey;
    jsonObject[ K::unidimensionalArrayInference ] = applicationOptions->unidimensionalArrayInference;
    assignIfNotEmpty( jsonObject, K::populatedFile, applicationOptions->populatedFile );
}

//! Create a shared pointer to a `ApplicationOptions` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< ApplicationOptions >& applicationOptions )
{
    using K = Keys::Options;

    applicationOptions = boost::make_shared< ApplicationOptions >( );
    updateFromJSONIfDefined( applicationOptions->defaultValueUsedForMissingKey,
                             jsonObject, K::defaultValueUsedForMissingKey );
    updateFromJSONIfDefined( applicationOptions->unusedKey, jsonObject, K::unusedKey );
    updateFromJSONIfDefined( applicationOptions->unidimensionalArrayInference,
                             jsonObject, K::unidimensionalArrayInference );
    updateFromJSONIfDefined( applicationOptions->populatedFile, jsonObject, K::populatedFile );
}

} // namespace json_interface

} // namespace tudat
