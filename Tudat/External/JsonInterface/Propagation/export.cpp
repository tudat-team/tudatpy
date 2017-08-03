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

#include "export.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `ExportSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< ExportSettings >& exportSettings )
{
    if ( ! exportSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Export;

    jsonObject[ K::file ] = exportSettings->outputFile;
    jsonObject[ K::variables ] = exportSettings->variables;
    jsonObject[ K::epochsInFirstColumn ] = exportSettings->epochsInFirstColumn;
    jsonObject[ K::onlyFinalStep ] = exportSettings->onlyFinalStep;
    jsonObject[ K::numericalPrecision ] = exportSettings->numericalPrecision;
}

//! Create a shared pointer to a `ExportSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< ExportSettings >& exportSettings )
{
    using namespace json_interface;
    using K = Keys::Export;

    exportSettings = boost::make_shared< ExportSettings >( getValue< path >( jsonObject, K::file ),
                                                           getVariables( jsonObject, K::variables ) );
    updateFromJSONIfDefined( exportSettings->epochsInFirstColumn, jsonObject, K::epochsInFirstColumn );
    updateFromJSONIfDefined( exportSettings->onlyFinalStep, jsonObject, K::onlyFinalStep );
    updateFromJSONIfDefined( exportSettings->numericalPrecision, jsonObject, K::numericalPrecision );
}

} // namespace simulation_setup

} // namespace tudat
