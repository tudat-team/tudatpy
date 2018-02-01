/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <json/src/json.hpp>

namespace tudat
{

namespace json_interface
{

//! Read a JSON file into a `json` object.
/*!
 * @copybrief readJSON If \p inclusionFile is not empty, the strings matching the expression "#path(relativePath)"
 * are replaced by the corresponding relative path relative to the file in which it is being included.
 * \param filePath Path to the JSON file.
 * \param parentFilePath File from which the contents of \p filePath are being read
 * (empty if \p filePath is equal to \p parentFilePath).
 * \param rootFilePath Path to the root file that was passed as command-line argument to the application
 * (empty if \p filePath is equal to \p rootFilePath).
 * \return Read `json` object.
 */
nlohmann::json readJSON( const boost::filesystem::path& filePath, const boost::filesystem::path& parentFilePath = boost::filesystem::path( ), const boost::filesystem::path& rootFilePath = boost::filesystem::path( ) );

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
 * \throws Error if the file does not exist.
 */
boost::filesystem::path getPathForJSONFile( const std::string& file, const boost::filesystem::path& basePath = boost::filesystem::current_path( ) );

//! Read and parse a (normal) `json` object from a file, and then parse its imported modular files.
/*!
 * Read and parse a (normal) `json` object from a file, and then parse its imported modular files using
 * parseModularJSON().
 * \param inputFilePath Path to the root JSON file.
 * \return Object containing the keys defined in the original file and all the imported files.
 */
nlohmann::json getDeserializedJSON( const boost::filesystem::path& inputFilePath );

}  // namespace json_interface

}  // namespace tudat

#endif // TUDAT_JSONINTERFACE_MODULAR_H
