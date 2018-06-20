/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/JsonInterface/Support/options.h"

namespace tudat
{

namespace json_interface
{

//! Create a `json` object from a shared pointer to a `ApplicationOptions` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< ApplicationOptions >& applicationOptions )
{
    if ( ! applicationOptions )
    {
        return;
    }
    using K = Keys::Options;

    jsonObject[ K::notifyOnPropagationStart ] = applicationOptions->notifyOnPropagationStart_;
    jsonObject[ K::notifyOnPropagationTermination ] = applicationOptions->notifyOnPropagationTermination_;
    jsonObject[ K::defaultValueUsedForMissingKey ] = applicationOptions->defaultValueUsedForMissingKey_;
    jsonObject[ K::unusedKey ] = applicationOptions->unusedKey_;
    assignIfNotEmpty( jsonObject, K::fullSettingsFile, applicationOptions->fullSettingsFile_ );
    jsonObject[ K::tagOutputFilesIfPropagationFails ] = applicationOptions->tagOutputFilesIfPropagationFails_;
}

//! Create a shared pointer to a `ApplicationOptions` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< ApplicationOptions >& applicationOptions )
{
    using K = Keys::Options;

    applicationOptions = std::make_shared< ApplicationOptions >( );

    updateFromJSONIfDefined( applicationOptions->notifyOnPropagationStart_, jsonObject, K::notifyOnPropagationStart );

    updateFromJSONIfDefined( applicationOptions->notifyOnPropagationTermination_,
                             jsonObject, K::notifyOnPropagationTermination );

    updateFromJSONIfDefined( applicationOptions->defaultValueUsedForMissingKey_,
                             jsonObject, K::defaultValueUsedForMissingKey );

    updateFromJSONIfDefined( applicationOptions->unusedKey_, jsonObject, K::unusedKey );

    updateFromJSONIfDefined( applicationOptions->fullSettingsFile_, jsonObject, K::fullSettingsFile );

    updateFromJSONIfDefined( applicationOptions->tagOutputFilesIfPropagationFails_,
                             jsonObject, K::tagOutputFilesIfPropagationFails );
}

} // namespace json_interface

} // namespace tudat
