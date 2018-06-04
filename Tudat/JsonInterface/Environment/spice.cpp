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

#include "Tudat/JsonInterface/Environment/spice.h"

namespace tudat
{

namespace json_interface
{

//! Create a `json` object from a shared pointer to a `SpiceSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< SpiceSettings >& spiceSettings )
{
    if ( ! spiceSettings )
    {
        return;
    }
    using K = Keys::Spice;

    jsonObject[ K::useStandardKernels ] = spiceSettings->useStandardKernels_;
    if ( spiceSettings->useStandardKernels_ )
    {
        assignIfNotEmpty( jsonObject, K::alternativeKernels, spiceSettings->alternativeKernels_ );
    }
    else
    {
        jsonObject[ K::kernels ] = spiceSettings->kernels_;
    }

    jsonObject[ K::preloadEphemeris ] = spiceSettings->preloadEphemeris_;
    if ( spiceSettings->preloadEphemeris_ )
    {
        jsonObject[ K::interpolationStep ] = spiceSettings->interpolationStep_;
        jsonObject[ K::interpolationOffsets ] = { spiceSettings->getInitialOffset( ),
                                                  spiceSettings->getFinalOffset( ) };
    }
}

//! Create a shared pointer to a `SpiceSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< SpiceSettings >& spiceSettings )
{
    using K = Keys::Spice;

    spiceSettings = std::make_shared< SpiceSettings >( );

    updateFromJSON( spiceSettings->useStandardKernels_, jsonObject, K::useStandardKernels );
    if ( spiceSettings->useStandardKernels_ )
    {
        updateFromJSONIfDefined( spiceSettings->alternativeKernels_, jsonObject, K::alternativeKernels );
    }
    else
    {
        updateFromJSON( spiceSettings->kernels_, jsonObject, K::kernels );
    }

    updateFromJSONIfDefined( spiceSettings->preloadEphemeris_, jsonObject, K::preloadEphemeris );
    if ( spiceSettings->preloadEphemeris_ )
    {
        spiceSettings->interpolationOffsets_ =
                getValue( jsonObject, K::interpolationOffsets, spiceSettings->interpolationOffsets_ );
        spiceSettings->interpolationStep_ =
                getValue( jsonObject, K::interpolationStep, spiceSettings->interpolationStep_ );
    }
}


//! Load in Tudat the Spice kernels specified in \p spiceSettings.
void loadSpiceKernels( const std::shared_ptr< SpiceSettings >& spiceSettings )
{
    using namespace spice_interface;

    clearSpiceKernels( );

    if ( spiceSettings )
    {
        if ( spiceSettings->useStandardKernels_ )
        {
            std::vector< std::string > alternativeKernelsFiles;
            for ( boost::filesystem::path kernel : spiceSettings->alternativeKernels_ )
            {
                alternativeKernelsFiles.push_back( kernel.string( ) );
            }
            loadStandardSpiceKernels( alternativeKernelsFiles );
        }
        else
        {
            for ( const boost::filesystem::path kernel : spiceSettings->kernels_ )
            {
                spice_interface::loadSpiceKernelInTudat( kernel.string( ) );
            }
        }
    }
}

} // namespace json_interface

} // namespace tudat
