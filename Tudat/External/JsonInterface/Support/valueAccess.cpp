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

#include "valueAccess.h"

#include "options.h"

namespace tudat
{

namespace json_interface
{

// KEY ACCESS

//! Get an index for a `json` array from a int-convertible key
/*!
 * @copybrief indexFromKey
 * \remark If \p key is a negative integer, reverse access is used (e.g. -1 = last element).
 * \param key The int-convertible key.
 * \param jsonArray The `json` array.
 * \return The array index.
 * \throws std::exception If \p key is not int-convertible.
 */
int indexFromKey( const std::string& key, const json& jsonArray )
{
    int arrayIndex = std::stoi( key );
    if ( arrayIndex < 0 )
    {
        arrayIndex += jsonArray.size( );
    }
    return arrayIndex;
}

//! Access/modify a key of a `json` object or array.
json& valueAt( json& jsonObject, const std::string& key )
{
    try
    {
        return jsonObject[ key ];
    }
    catch ( ... )
    {
        return jsonObject[ indexFromKey( key, jsonObject ) ];
    }
}

//! Access a key of a `json` object or array.
const json& valueAt( const json& jsonObject, const std::string& key )
{
    try
    {
        return jsonObject.at( key );
    }
    catch ( ... )
    {
        return jsonObject.at( indexFromKey( key, jsonObject ) );
    }
}

//! Access a key path of a `json` object or array.
json valueAt( json jsonObject, const KeyPath& keyPath )
{
    for ( const std::string key : keyPath )
    {
        const json constJsonObject = jsonObject;
        jsonObject = valueAt( constJsonObject, key );
    }
    return jsonObject;
}

//! Whether the key at \p keyPath is defined for \p jsonObject.
bool defined( const json& jsonObject, const KeyPath& keyPath )
{
    try
    {
        valueAt( jsonObject, keyPath );
        return true;
    }
    catch ( ... )
    {
        return false;
    }
}


// SPECIAL KEYS ACCESS

//! Get the a shared pointer to \p jsonObject at key SpecialKeys::rootObject.
json getRootObject( const json& jsonObject )
{
    return getValue< json >( jsonObject, SpecialKeys::rootObject );
}

//! Get the absolute key path from which \p jsonObject was retrieved.
KeyPath getKeyPath( const json& jsonObject )
{
    return getValue< KeyPath >( jsonObject, SpecialKeys::keyPath, SpecialKeys::root );
}

//! Get the key at which \p jsonObject was obtained.
std::string getParentKey( const json& jsonObject, const std::string& errorMessage )
{
    try
    {
        return getKeyPath( jsonObject ).back( );
    }
    catch ( ... )
    {
        std::cerr << errorMessage << std::endl;
        throw;
    }
}

//! Convert \p j to object if \p j is array.
void convertToObjectIfArray( json& j, const bool onlyIfElementsAreStructured )
{
    if ( ! j.is_array( ) )
    {
        return;
    }
    if ( onlyIfElementsAreStructured )
    {
        if ( j.empty( ) )
        {
            return;
        }
        if ( ! j.front( ).is_structured( ) )
        {
            return;
        }
    }
    json jsonArray = j;
    j = json( );
    for ( unsigned int i = 0; i < jsonArray.size( ); ++i )
    {
        j[ std::to_string( i ) ] = jsonArray.at( i );
    }
}

//! Get the response type to an event for a `json` object.
ExceptionResponseType getResponseToEventNamed( const json& jsonObject, const std::string& eventName,
                                               const ExceptionResponseType defaultResponse )
{
    ExceptionResponseType response = defaultResponse;
    try
    {
        json rootObject = jsonObject;
        try
        {
            rootObject = valueAt( jsonObject, SpecialKeys::rootObject );
        }
        catch ( ... ) { }
        response = valueAt( rootObject, Keys::options / eventName ).get< ExceptionResponseType >( );
    }
    catch ( ... ) { }
    return response;
}


// ACCESS HISTORY

//! Global variable containing all the key paths that were accessed since clearAccessHistory() was called for the
//! last time (or since this variable was initialized).
std::set< KeyPath > accessedKeyPaths = { };

//! Get all the key paths defined for \p jsonObject.
/*!
 * @copybrief getKeyPaths
 * \param jsonObject The `json` object.
 * \param baseKeyPath The key path in which to look for sub-key paths. Default is SpecialKeys::root, i.e. a key path
 * will be added to the returned `std::set` for each key defined in \p jsonObject.
 * \return A set containing a key path for each key defined in \p jsonObject at \p baseKeyPath.
 */
std::set< KeyPath > getKeyPaths( const json& jsonObject, const KeyPath& baseKeyPath = SpecialKeys::root )
{
    std::set< KeyPath > keyPaths;
    for ( json::const_iterator it = jsonObject.begin( ); it != jsonObject.end( ); ++it )
    {
        const std::string key = it.key( );
        json subObject = it.value( );
        KeyPath keyPath = baseKeyPath / key;
        keyPaths.insert( keyPath );
        convertToObjectIfArray( subObject, true );
        if ( subObject.is_object( ) )
        {
            const std::set< KeyPath > subkeyPaths = getKeyPaths( subObject, keyPath );
            for ( const KeyPath subkey : subkeyPaths )
            {
                keyPaths.insert( subkey );
            }
        }
    }
    return keyPaths;
}

//! Check for key paths that are defined in \p jsonObject but not contained by the global variable accessedKeyPaths.
void checkUnusedKeys( const json& jsonObject, const ExceptionResponseType response )
{
    if ( response == printWarning || response == throwError )
    {
        bool unusedKeys = false;
        const std::set< KeyPath > definedKeyPaths = getKeyPaths( jsonObject );
        for ( const KeyPath keyPath : definedKeyPaths )
        {
            if ( accessedKeyPaths.count( keyPath ) == 0 )
            {
                unusedKeys = true;
                std::cerr << "Unused key: " << keyPath << std::endl;
            }
        }
        if ( unusedKeys && response == throwError )
        {
            throw std::runtime_error( "Validation failed because there are unsued keys." );
        }
    }
}

}  // namespace json_interface

}  // namespace tudat
