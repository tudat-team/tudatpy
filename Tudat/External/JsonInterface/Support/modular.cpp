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
#include "valueAccess.h"

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

//! Update paths for a JSON object.
/*!
 * @copybrief updatePaths
 * <br/>
 * Looks for expressions ${FILE_STEM}, ${FILE_NAME}, ${PARENT_FILE_STEM}, ${PARENT_FILE_NAME},
 * ${ROOT_FILE_STEM} and ${ROOT_FILE_NAME} in all the strings of \p jsonObject, replacing them, respetively, bu the
 * stem of filename of the file in which it is defined, in which it was included, or in the root input file (the one
 * provided to the command line app as argument).
 * <br/>
 * Then, it fixes paths in included files by searching string matching the expressions "@path(myPath)", only when
 * myPath is a relative path. This is necessary because myPath must be relative to the root file, but in input files
 * the user provides paths relative to the directory in which the definition file is located.
 * <br/>
 * This feature can be "turned off" by providing paths directly without the "@path" keyword (e.g. "myPath"). This
 * is safe when no modular files are included, but can result in problems otherwise. If paths are provided directly
 * without the "@path" keyword, then the user is responsible for defining them always relative to the root file in
 * which they are going to be included (and not relative to the file in which they are defined).
 * \param jsonObject The `json` object containing / being a string.
 * \param filePath Path to the file in which the path is defined.
 * \param parentFilePath Path to the file in which the file in which the path is defined was included
 * (equal to \p filePath if \p jsonObject was defined at \p parentFilePath).
 * \param rootFilePath Path to the root file that was passed as command-line argument to the application.
 * (equal to \p filePath if \p jsonObject was defined at \p rootFilePath).
 */
void updatePaths( json& jsonObject, const path& filePath, const path& parentFilePath, const path& rootFilePath )
{
    if ( jsonObject.is_structured( ) )
    {
        for ( json::iterator it = jsonObject.begin( ); it != jsonObject.end( ); ++it )
        {
            json& subjson = it.value( );
            updatePaths( subjson, filePath, parentFilePath, rootFilePath );
        }
    }
    else if ( jsonObject.is_string( ) )
    {
        std::string str = jsonObject.get< std::string >( );

        // Replace: ${FILE_DIR}, ${FILE_STEM}, ${FILE_NAME}, etc.
        str = regex_replace( str, boost::regex( R"(\$\{FILE_DIR\})" ), filePath.parent_path( ).string( ) );
        str = regex_replace( str, boost::regex( R"(\$\{FILE_STEM\})" ), filePath.stem( ).string( ) );
        str = regex_replace( str, boost::regex( R"(\$\{FILE_NAME\})" ), filePath.filename( ).string( ) );
        str = regex_replace( str, boost::regex( R"(\$\{PARENT_FILE_DIR\})" ), parentFilePath.parent_path( ).string( ) );
        str = regex_replace( str, boost::regex( R"(\$\{PARENT_FILE_STEM\})" ), parentFilePath.stem( ).string( ) );
        str = regex_replace( str, boost::regex( R"(\$\{PARENT_FILE_NAME\})" ), parentFilePath.filename( ).string( ) );
        str = regex_replace( str, boost::regex( R"(\$\{ROOT_FILE_DIR\})" ), rootFilePath.parent_path( ).string( ) );
        str = regex_replace( str, boost::regex( R"(\$\{ROOT_FILE_STEM\})" ), rootFilePath.stem( ).string( ) );
        str = regex_replace( str, boost::regex( R"(\$\{ROOT_FILE_NAME\})" ), rootFilePath.filename( ).string( ) );

        // Fix relative paths
        boost::cmatch groups;
        boost::regex_match( str.c_str( ), groups, boost::regex( R"(\@path\((.+)\))" ) );
        if ( groups[ 1 ].matched )
        {
            const path providedPath = std::string( groups[ 1 ] );
            if ( providedPath.is_relative( ) )
            {
                const path rel = boost::filesystem::relative( filePath.parent_path( ), rootFilePath.parent_path( ) );
                str = ( rel / providedPath ).string( );
            }
            else
            {
                const path rel = boost::filesystem::relative( providedPath, rootFilePath.parent_path( ) );
                str = ( ( ! rel.empty( ) && rel.size( ) < providedPath.size( ) ) ? rel : providedPath ).string( );
            }
        }

        jsonObject = str;
    }
}

//! Read a JSON file into a `json` object.
json readJSON( const path& filePath, const path& parentFilePath, const path& rootFilePath )
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
    updatePaths( jsonObject, filePath, parentFilePath.empty( ) ? filePath : parentFilePath,
                 rootFilePath.empty( ) ? filePath : rootFilePath );
    return jsonObject;
}

//! Merge an array of JSON input objects.
/*!
 * @copybrief mergeJSON
 * <br/>
 * If \p jsonObject is not an array, this function will do nothing.
 * <br/>
 * If \p jsonObject is an array of objects, the first object will be used to create the initial `json` object.
 * Then, the subsequent objects will be used to create this initial object, by updating (or settings) the keys
 * defined in these subsequent objects with the provided objects.
 * \param jsonObject The `json` containing an object (no merging will be applied) or an array of objects
 * (the elements will be merged).
 * \param filePath The file from which \p jsonObject was retrieved.
 */
void mergeJSON( json& jsonObject, const path& filePath )
{
    if ( jsonObject.is_array( ) )
    {
        json jsonArray = jsonObject;
        for ( json::iterator it = jsonArray.begin( ); it != jsonArray.end( ); ++it )
        {
            json subjson = it.value( );
            mergeJSON( subjson, filePath );
            if ( it == jsonArray.begin( ) )
            {
                jsonObject = subjson;
            }
            else
            {
                for ( json::iterator subit = subjson.begin( ); subit != subjson.end( ); ++subit )
                {
                    const KeyPath keyPath = split( subit.key( ), '.' );
                    try
                    {
                        json newValue = subit.value( );
                        for ( unsigned int j = 0; j < keyPath.size( ); ++j )
                        {
                            json updatedsonObject = jsonObject;
                            for ( unsigned int i = 0; i < keyPath.size( ) - j; ++i )
                            {
                                const std::string key = keyPath.at( i );
                                if ( i < keyPath.size( ) - 1 - j )
                                {
                                    updatedsonObject = valueAt( updatedsonObject, key );
                                }
                                else
                                {
                                    valueAt( updatedsonObject, key ) = newValue;
                                }
                            }
                            newValue = updatedsonObject;
                        }
                        jsonObject = newValue;
                    }
                    catch ( ... )
                    {
                        std::cerr << "Could not update value for " << keyPath << " referenced from file "
                                  << filePath << std::endl;
                        throw;
                    }
                }
            }
        }
    }
}

//! Parse a modular `json` object.
/*!
 * Parse a modular `json` object containing strings (or being equal to a string) of matching the expression
 * `"$(path){key}"`, which will be replaced with the value for the key "key" in the file at "path".
 * <br/>
 * "path" can be absolute or relative (to the parent directory of \p filePath).
 * <br/>
 * If the part "(path)" is not provided, the key "key" will be searched in \p parentObject
 * (or in \p jsonObject if \p parentObject is not defined).
 * <br/>
 * If the part "{key}" is not provided, the whole `json` object contained in the file at "path" will be used.
 * <br/>
 * "key" can be a single key, or a key path separated by dots (e.g. "key.subkey.subsubkey").
 * Elements of arrays can be accessed by indicating their index, starting from 0 (e.g. "key.2.subkey.0").
 * <br/>
 * "key" can contain a single key (path) or many, separated by commas (e.g. "key1,key2.subkey"). If the character ","
 * is found in "key", \p jsonObject will be converted to an array. The key "key1," creates an array with one element.
 * <br/>
 * "key" can contain colons, in which case \p jsonObject will be converted to an object/map (e.g. "a:keyA,b:keyB.0").
 * <br/>
 * Illegal characters in referenced file paths: (){}
 * <br/>
 * Illegal characters in referenced keys: (){}.:,
 * \remark This function results in a recursive loop when a key value references itself, making the app terminate
 * without providing an informative error message.
 * \param jsonObject The `json` object to be parsed (passed by reference).
 * \param filePath Path of the file from which \p jsonObject was retrieved.
 * \param parentObject The `json` object form which \p jsonObject was retrieved, if any (default is null `json`).
 * \param rootFilePath Path to the root file that was passed as command-line argument to the application
 * (empty if \p filePath is equal to \p rootFilePath).
 */
void parseModularJSON( json& jsonObject, const path& filePath,
                       json parentObject = json( ), path rootFilePath = path( ) )
{
    if ( parentObject.is_null( ) )
    {
        parentObject = jsonObject;
    }
    if ( rootFilePath.empty( ) )
    {
        rootFilePath = filePath;
    }
    if ( jsonObject.is_structured( ) )
    {
        for ( json::iterator it = jsonObject.begin( ); it != jsonObject.end( ); ++it )
        {
            json& subjson = it.value( );
            parseModularJSON( subjson, filePath, jsonObject, rootFilePath );
        }
    }
    else if ( jsonObject.is_string( ) )
    {
        boost::cmatch groups;
        boost::regex_match( jsonObject.get< std::string >( ).c_str( ), groups,
                            boost::regex( R"(\$(?:\((.+?)\))?(?:\{(.+?)\})?)" ) );
        const bool fileMatch = groups[ 1 ].matched;
        const bool varsMatch = groups[ 2 ].matched;
        if ( fileMatch || varsMatch )
        {
            const std::string file( groups[ 1 ] );
            const std::string vars( groups[ 2 ] );
            const path importPath = fileMatch ? getPathForJSONFile( file, filePath.parent_path( ) ) : filePath;
            const json importedJsonObject = fileMatch ? readJSON( importPath, filePath, rootFilePath ) : parentObject;
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
                    parseModularJSON( subJsonObject, importPath, parentObject, rootFilePath );
                }
                else
                {
                    try
                    {
                        // Recursively update jsonObject for every key in keyPath
                        for ( const std::string key : keyPath )
                        {
                            subJsonObject = valueAt( subJsonObject, key );
                            parseModularJSON( subJsonObject, importPath, parentObject, rootFilePath );
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

//! Read and parse a (normal) `json` object from a file, and then parse its imported modular files.
json getParsedModularJSON( const path& filePath )
{
    json jsonObject = readJSON( filePath );
    parseModularJSON( jsonObject, filePath );
    mergeJSON( jsonObject, filePath );
    return jsonObject;
}

}  // namespace json_interface

}  // namespace tudat
