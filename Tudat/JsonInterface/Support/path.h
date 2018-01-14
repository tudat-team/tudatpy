/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <json/src/json.hpp>

#include "Tudat/InputOutput/basicInputOutput.h"
namespace tudat
{

namespace json_interface
{

//! Names for replacable paths in JSON files.
//! E.g., the text "${TUDAT_ROOT_PATH}" will be replaced with the actual path when reading JSON files.
static std::vector< std::pair< std::string, std::string > > pathPlaceholders =
{
    std::make_pair( "TUDAT_ROOT_PATH", input_output::getTudatRootPath( ) ),
    std::make_pair( "ATMOSPHERE_TABLES_PATH", input_output::getAtmosphereTablesPath( ) ),
    std::make_pair( "GRAVITY_MODELS_PATH", input_output::getGravityModelsPath( ) ),
    std::make_pair( "SPICE_KERNELS_PATH", input_output::getSpiceKernelPath( ) )
};

//! Return \p path with the recognized paths replaced by placeholders (such as ${TUDAT_ROOT_PATH}).
/*!
 * @copybrief pathAddingPlaceholders
 * \param path Path without placeholders.
 * \return Path with placeholders.
 */
std::string pathAddingPlaceholders( std::string path );

//! Return \p path with the recognized path placeholders (such as ${TUDAT_ROOT_PATH}) replaced by the actual paths.
/*!
 * @copybrief pathRemovingPlaceholders
 * If any of the placeholders is not recognized, it will not be replaced and no warning or error will be generated.
 * \param path Path with placeholders.
 * \return Path without placeholders.
 */
std::string pathRemovingPlaceholders( std::string path );

} // namespace json_interface

} // namespace tudat


namespace boost
{

namespace filesystem
{

//! Create a `json` object from a `path`.
void to_json( nlohmann::json& j, const path& p );

//! Create a path from a `json` object.
void from_json( const nlohmann::json& j, path& p );

}  // namespace filesystem

}  // namespace boost

#endif // TUDAT_JSONINTERFACE_PATH_H
