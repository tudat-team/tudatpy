/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_MODULAR_H
#define TUDAT_JSONINTERFACE_MODULAR_H

#include <boost/filesystem.hpp>

#include "json/src/json.hpp"
using json = nlohmann::json;

namespace tudat
{

namespace json_interface
{

typedef boost::filesystem::path path;


//! Get the path for a JSON file.
/*!
 * Get the path for a JSON file.
 * This function will first try to locate the file "basePath/file".
 * If this fails, and `file` contains no extension, the function will try to locate the file "basePath/file/main.json".
 * If this also fails, the function will then try to locate the file "basePath/file.json".
 * If this also fails, an error will be thrown.
 * \param file Name of the JSON file (with extension) or directory at which a "main.json" file is located.
 * \param basePath Parent directory in which to look for the requested file.
 * \return Path for the JSON file.
 * \throw Error if the file does not exist.
 */
path getPathForJSONFile( const std::string& file, const path& basePath = boost::filesystem::current_path( ) );


//! Parse a modular `json` object containing #import commands.
/*!
 * Parses a modular `json` object, by recursively iterating over all the values of the keys defined in `jsonObject`
 * that are either an object or an array of objects, looking inside those objects for keys with string values following
 * the expression "#import relativePathToJsonFile", and replacing them with the contents of "relativePathToJsonFile".
 *
 * The number of spaces between #import and relativePathToJsonFile is irrelevant, but must be at least 1.
 * relativePathToJsonFile can contain spaces. Some characters must be escaped (e.g. \").
 *
 * Flattening will be applied during the replacement if (and only if) in relativePathToJsonFile only one key is defined,
 * with the same name as the key that it is replacing.
 *
 * Example WITHOUT FLATTENING:
 * jsonObject:
 * {
 *   "age": 25,
 *   "person": "#import someone.json"
 * }
 * someone.json:
 * {
 *   "name": "Aleix",
 *   "surname": "Pinardell"
 * }
 * returned json object:
 * {
 *   "age": 25,
 *   "person":
 *   {
 *     "name": "Aleix",
 *     "surname": "Pinardell
 *   }
 * }
 *
 * Example WITH FLATTENING:
 * jsonObject:
 * {
 *   "age": 25,
 *   "name": "#import someone.json"
 * }
 * someone.json:
 * {
 *   "name": "Aleix Pinardell"
 * }
 * returned json object:
 * {
 *   "age": 25,
 *   "name": "Aleix Pinardell
 * }
 *
 * \param jsonObject
 * \param parentDirectoryPath
 */
void parseModularJSON( json& jsonObject, const path& parentDirectoryPath );


//! -DOC
json getParsedModularJSON( const path& inputFilePath );

}  // namespace json_interface

}  // namespace tudat

#endif // TUDAT_JSONINTERFACE_MODULAR_H
