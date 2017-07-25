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

#include "jsonInterface.h"

namespace tudat
{

namespace json_interface
{

//! -DOC
path pathForJSONFile( const std::string& file, const path& basePath )
{
    try
    {
        // Get absolute path to input file
        path filePath = boost::filesystem::canonical( file, basePath );
        if ( boost::filesystem::is_directory( filePath ) )
        {
            // If it's a directory, try to find a JSON file with the same name in that directory
            // E.g. if asked to load "bodies", and "bodies" is a directory, try to load "bodies/bodies.json"
            return pathForJSONFile( file, filePath );
        }
        else
        {
            return filePath;
        }
    }
    catch ( ... )
    {
        try
        {
            // Add .json extension if necessary
            if ( file.find( ".json" ) == std::string::npos )
            {
                return boost::filesystem::canonical( file + ".json", basePath );
            }
            else
            {
                throw;
            }
        }
        catch ( ... )
        {
            throw std::runtime_error( "The file \"" + file +  "\" does not exist." );
        }
    }
}

//! -DOC
json parseModularJSONFile( const path& filePath )
{
    json jsonObject = json::parse( std::ifstream( filePath.string( ) ) );
    json jsonObjectWithImportedFiles = jsonObject;

    for ( json::iterator it = jsonObject.begin( ); it != jsonObject.end( ); ++it )
    {
        std::string importCommand;
        try
        {
            importCommand = it.value( ).get< std::string >( );
        }
        catch ( ... )
        {
            continue;
        }
        boost::regex expression { "#import\\s+(.+)" };
        boost::regex_token_iterator< std::string::iterator > bit{
            importCommand.begin( ), importCommand.end( ), expression, 1 };
        boost::regex_token_iterator< std::string::iterator > bend;
        while ( bit != bend )
        {
            const std::string fileToImport = *bit++;
            jsonObjectWithImportedFiles[ it.key( ) ] =
                    parseModularJSONFile( pathForJSONFile( fileToImport, filePath.parent_path( ) ) );
        }
    }

    return jsonObjectWithImportedFiles;
}


//! -DOC
bool defined( const json& jsonObject, const KeyTree& keyTree )
{
    if ( getValuePointer< json >( jsonObject, keyTree ) )
    {
        return true;
    }
    else
    {
        return false;
    }
}


}  // namespace json_interface

}  // namespace tudat
