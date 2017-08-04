/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_PATH_H
#define TUDAT_JSONINTERFACE_PATH_H

#include <boost/filesystem.hpp>

#include "json/src/json.hpp"
using json = nlohmann::json;

#include <Tudat/InputOutput/basicInputOutput.h>

namespace tudat
{

namespace json_interface
{

typedef boost::filesystem::path path;

//! -DOC
static std::vector< std::pair< std::string, std::string > > pathPlaceholders =
{
    std::make_pair( "TUDAT_ROOT_PATH", input_output::getTudatRootPath( ) ),
    std::make_pair( "ATMOSPHERE_TABLES_PATH", input_output::getAtmosphereTablesPath( ) ),
    std::make_pair( "GRAVITY_MODELS_PATH", input_output::getGravityModelsPath( ) ),
    std::make_pair( "SPICE_KERNELS_PATH", input_output::getSpiceKernelPath( ) )
};

//! Replace recognized paths with placeholders such as ${TUDAT_ROOT_PATH}.
std::string addPathPlaceholders( std::string path );

//! Replace recognized placeholders such as ${TUDAT_ROOT_PATH} with the actual paths.
std::string removePathPlaceholders( std::string path );

} // namespace json_interface

} // namespace tudat


namespace boost
{

namespace filesystem
{

//! Create a `json` object from a `path`.
void to_json( json& j, const path& p );

//! Create a path from a `json` object.
void from_json( const json& j, path& p );

}  // namespace filesystem

}  // namespace boost

#endif // TUDAT_JSONINTERFACE_PATH_H
