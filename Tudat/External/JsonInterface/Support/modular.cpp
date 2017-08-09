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

#include <boost/regex.hpp>

#include "modular.h"
#include "keys.h"
#include "utilities.h"

namespace tudat
{

namespace json_interface
{

//! Get the path for a JSON file.
path getPathForJSONFile( const std::string& file, const path& basePath )
{
    try
    {
        // Get absolute path to input file
        path filePath = boost::filesystem::canonical( file, basePath );
        if ( boost::filesystem::is_directory( filePath ) )
        {
            // If it's a directory, try to find a "main.json" file in that directory
            return getPathForJSONFile( "main.json", filePath );
        }
        else
        {
            return filePath;
        }
    }
    catch ( ... )
    {
        std::string fileWithExtension = file;
        try
        {
            // Add .json extension if no extension found
            if ( path( file ).extension( ).empty( ) )
            {
                fileWithExtension += ".json";
                return boost::filesystem::canonical( fileWithExtension, basePath );
            }
            else
            {
                throw;
            }
        }
        catch ( ... )
        {
            std::cerr << "The file " << basePath / fileWithExtension << " does not exist." << std::endl;
            throw;
        }
    }
}

//! Parse a modular `json` object.
void parseModularJSON( json& jsonObject, const path& filePath, json parentObject = json( ) )
{
    if ( parentObject.is_null( ) )
    {
        parentObject = jsonObject;
    }
    if ( jsonObject.is_structured( ) )
    {
        for ( json::iterator it = jsonObject.begin( ); it != jsonObject.end( ); ++it )
        {
            json& subjson = it.value( );
            parseModularJSON( subjson, filePath, jsonObject );
        }
    }
    else if ( jsonObject.is_string( ) )
    {
        boost::regex expression( R"(\$(?:\((.+?)\))?(?:\{(.+?)\})?)" );
        boost::cmatch groups;
        boost::regex_match( jsonObject.get< std::string >( ).c_str( ), groups, expression );
        const bool fileMatch = groups[ 1 ].matched;
        const bool varsMatch = groups[ 2 ].matched;
        if ( fileMatch || varsMatch )
        {
            const std::string file( groups[ 1 ] );
            const std::string vars( groups[ 2 ] );
            const path importPath = fileMatch ? getPathForJSONFile( file, filePath.parent_path( ) ) : filePath;
            const json importedJsonObject =
                    fileMatch ? json::parse( std::ifstream( importPath.string( ) ) ) : parentObject;
            std::vector< std::string > keys;
            std::vector< KeyPath > keyPaths;
            if ( varsMatch )
            {
                for ( const std::string variable : split( vars, ',' ) )
                {
                    const std::vector< std::string > keyVar = split( variable, ':' );
                    if ( keyVar.size( ) == 2 )
                    {
                        keys.push_back( keyVar.front( ) );
                    }
                    keyPaths.push_back( split( keyVar.back( ), '.' ) );
                }
            }
            else
            {
                keyPaths = { KeyPath( { } ) };
            }
            json parsedJsonObject;
            for ( unsigned int i = 0; i < keyPaths.size( ); ++i )
            {
                const KeyPath keyPath = keyPaths.at( i );
                json subJsonObject = importedJsonObject;
                if ( keyPath.empty( ) )
                {
                    parseModularJSON( subJsonObject, importPath, parentObject );
                }
                else
                {
                    try
                    {
                        // Recursively update jsonObject for every key in keyPath
                        for ( const std::string key : keyPath )
                        {
                            try
                            {
                                // Try to access element at key
                                subJsonObject = subJsonObject.at( key );
                            }
                            catch ( ... )
                            {
                                // Key may be convertible to int.
                                subJsonObject = subJsonObject.at( std::stoi( key ) );
                            }
                            parseModularJSON( subJsonObject, importPath, parentObject );
                        }
                    }
                    catch ( ... )
                    {
                        std::cerr << "Could not load ";
                        if ( ! keyPath.empty( ) )
                        {
                            std::cerr << "value for " << keyPath << " from ";
                        }
                        std::cerr << "file " << importPath << " referenced from file " << filePath << std::endl;
                        throw;
                    }
                }
                if ( keys.size( ) > 0 )
                {
                    parsedJsonObject[ keys.at( i ) ] = subJsonObject;
                }
                else
                {
                    parsedJsonObject[ i ] = subJsonObject;
                }
            }
            if ( keyPaths.size( ) == 1 && keys.size( ) == 0 && vars.find( ',' ) == std::string::npos )
            {
                jsonObject = parsedJsonObject.front( );
            }
            else
            {
                jsonObject = parsedJsonObject;
            }
        }
    }
}

//! Get the corresponding line and column for a certain (byte) position in a `std::ifstream`.
/*!
 * Get the corresponding line and column for a certain (byte) position in a `std::ifstream`.
 * \param stream The stream (passed by reference, its current byte position will be set to zero and then reset to the
 * original one before returning, so effectively this function does not modify the stream).
 * \param position The byte position (or character index) of interest.
 * \return Line and column for the character at `position` in `stream`.
 */
std::pair< unsigned int, unsigned int > getLineAndCol( std::ifstream& stream, const std::streampos position )
{
    std::streampos originalPos = stream.tellg( );
    stream.seekg( 0, std::ifstream::beg );
    
    std::string lineContents;
    unsigned int line = 1;
    std::streampos col;
    std::streampos lastPos = stream.tellg( );
    while ( std::getline( stream, lineContents ) )
    {
        const std::streampos currentPos = stream.tellg( );
        if ( currentPos > position )
        {
            col = position - lastPos;
            break;
        }
        lastPos = currentPos;
        line++;
    }
    
    stream.seekg( originalPos );
    return { line, col };
}

//! Read and parse a (normal) `json` object from a file, and then parse its imported modular files.
json getParsedModularJSON( const path& filePath )
{
    std::ifstream stream( filePath.string( ) );
    json jsonObject;
    try
    {
        jsonObject = json::parse( stream );
    }
    catch ( const nlohmann::detail::parse_error& error )
    {
        std::pair< unsigned int, unsigned int > errorLineCol = getLineAndCol( stream, error.byte );
        std::cerr << "Parse error in file " << filePath
                  << " at line " << errorLineCol.first << ", col " << errorLineCol.second << "." << std::endl;
        throw error;
    }
    
    parseModularJSON( jsonObject, filePath );
    return jsonObject;
}

}  // namespace json_interface

}  // namespace tudat
