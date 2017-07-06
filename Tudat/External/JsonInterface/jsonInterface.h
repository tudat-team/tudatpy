/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_H
#define TUDAT_JSONINTERFACE_H

#include <boost/filesystem.hpp>

#include "json/src/json.hpp"

#include "integratorSettings.h"


namespace tudat
{

namespace json_interface
{

using json = nlohmann::json;
typedef boost::filesystem::path path;

json parseInputFile( std::string inputFile )
{
    // Add .json extension if necessary
    std::string inputFileWithExtension = inputFile;
    if ( inputFileWithExtension.find(".json") == std::string::npos )
    {
        inputFileWithExtension += ".json";
    }

    // Determine the absolute path of the directory where the input file is located
    path inputDirectory = boost::filesystem::canonical( inputFileWithExtension ).parent_path();

    // Determine the filename (without file extension) of the input file
    // path inputFilename = path( inputFile ).stem();

    // Build absolute path
    path inputPath = inputDirectory / path( inputFileWithExtension ).filename();

    // Check that input file exists
    if ( ! boost::filesystem::exists( inputPath ) )
    {
        std::cerr << "The input file \"" << inputPath << "\" does not exist." << std::endl;
        std::terminate();
    }

    // Read file
    std::ifstream stream( inputPath.string() );
    // std::string contents( ( std::istreambuf_iterator<char>(ifs) ), ( std::istreambuf_iterator<char>() ) );

    return json::parse( stream );
}

} // namespace json_interfaces

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_H
