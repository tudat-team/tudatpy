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

#include "Tudat/JsonInterface/Propagation/export.h"

namespace tudat
{

namespace json_interface
{

//! Create a `json` object from a shared pointer to a `ExportSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< ExportSettings >& exportSettings )
{
    if ( ! exportSettings )
    {
        return;
    }
    using K = Keys::Export;

    jsonObject[ K::file ] = exportSettings->outputFile_;
    jsonObject[ K::variables ] = exportSettings->variables_;
    assignIfNotEmpty( jsonObject, K::header, exportSettings->header_ );
    jsonObject[ K::epochsInFirstColumn ] = exportSettings->epochsInFirstColumn_;
    jsonObject[ K::onlyInitialStep ] = exportSettings->onlyInitialStep_;
    jsonObject[ K::onlyFinalStep ] = exportSettings->onlyFinalStep_;
    jsonObject[ K::numericalPrecision ] = exportSettings->numericalPrecision_;
}

//! Create a shared pointer to a `ExportSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< ExportSettings >& exportSettings )
{
    using namespace propagators;
    using K = Keys::Export;

    exportSettings = boost::make_shared< ExportSettings >(
                getValue< boost::filesystem::path >( jsonObject, K::file ),
                getValue< std::vector< boost::shared_ptr< VariableSettings > > >( jsonObject, K::variables ) );
    updateFromJSONIfDefined( exportSettings->header_, jsonObject, K::header );
    updateFromJSONIfDefined( exportSettings->epochsInFirstColumn_, jsonObject, K::epochsInFirstColumn );
    updateFromJSONIfDefined( exportSettings->onlyInitialStep_, jsonObject, K::onlyInitialStep );
    updateFromJSONIfDefined( exportSettings->onlyFinalStep_, jsonObject, K::onlyFinalStep );
    updateFromJSONIfDefined( exportSettings->numericalPrecision_, jsonObject, K::numericalPrecision );
}

} // namespace simulation_setup

} // namespace tudat
