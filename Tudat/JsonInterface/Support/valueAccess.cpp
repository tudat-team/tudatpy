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

#include "Tudat/JsonInterface/Support/valueAccess.h"

#include "Tudat/JsonInterface/Support/options.h"

namespace tudat
{

namespace json_interface
{

// KEY ACCESS

//! Access/mutate a key of a `json` object or array.
nlohmann::json& valueAt( nlohmann::json& jsonObject, const std::string& key, const bool mutator )
{
    const int intKey = indexFromKey( key );
    if ( mutator )
    {
        if ( intKey >= 0 && jsonObject.is_array( ) )
        {
            return jsonObject[ intKey ];
        }
        else
        {
            return jsonObject[ key ];
        }
    }
    else
    {
        if ( intKey >= 0 && jsonObject.is_array( )  )
        {
            return jsonObject.at( intKey );
        }
        else
        {
            return jsonObject.at( key );
        }
    }
}

//! Access a key path of a `json` object or array.
nlohmann::json valueAt( nlohmann::json jsonObject, const KeyPath& keyPath )
{
    for ( const std::string key : keyPath )
    {
        jsonObject = valueAt( jsonObject, key );
    }
    return jsonObject;
}

//! Whether the key at \p keyPath is defined for \p jsonObject.
bool isDefined( const nlohmann::json& jsonObject, const KeyPath& keyPath )
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
nlohmann::json getRootObject( const nlohmann::json& jsonObject )
{
    return getValue< nlohmann::json >( jsonObject, SpecialKeys::rootObject );
}

//! Get the absolute key path from which \p jsonObject was retrieved.
KeyPath getKeyPath( const nlohmann::json& jsonObject )
{
    return getValue( jsonObject, SpecialKeys::keyPath, KeyPath( SpecialKeys::root ) );
}

//! Get the key at which \p jsonObject was obtained.
std::string getParentKey( const nlohmann::json& jsonObject, const std::string& errorMessage )
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

//! Get the response type to an event for a `json` object.
ExceptionResponseType getResponseToEvent( const nlohmann::json& jsonObject, const std::string& eventName,
                                          const ExceptionResponseType defaultResponse )
{
    ExceptionResponseType response = defaultResponse;
    try
    {
        nlohmann::json rootObject = jsonObject;
        try
        {
            rootObject = valueAt( jsonObject, SpecialKeys::rootObject );
        }
        catch ( ... ) { }
        response = valueAt( rootObject, Keys::options / eventName );
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
std::set< KeyPath > getKeyPaths( const nlohmann::json& jsonObject, const KeyPath& baseKeyPath = SpecialKeys::root )
{
    std::set< KeyPath > keyPaths;
    for ( nlohmann::json::const_iterator it = jsonObject.begin( ); it != jsonObject.end( ); ++it )
    {
        const std::string key = it.key( );
        nlohmann::json subObject = it.value( );
        KeyPath keyPath = baseKeyPath / key;
        keyPaths.insert( keyPath );
        convertToObjectIfContainsObjects( subObject );
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
void checkUnusedKeys( const nlohmann::json& jsonObject, const ExceptionResponseType response )
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


// JSON ARRAY

//! Returns whether any of the elements in \p jsonArray (and their subelements) are of type object.
bool arrayContainsObjects( nlohmann::json& jsonArray )
{
    for ( unsigned int i = 0; i < jsonArray.size( ); ++i )
    {
        nlohmann::json element = jsonArray.at( i );
        if ( element.is_object( ) )
        {
            return true;
        }
        else if ( element.is_array( ) )
        {
            if ( arrayContainsObjects( element ) )
            {
                return true;
            }
        }
    }
    return false;
}

//! Convert \p j to object if \p j is an array containing objects.
void convertToObjectIfContainsObjects( nlohmann::json& j )
{
    if ( ! j.is_array( ) || j.empty( ) )
    {
        return;
    }
    if ( ! arrayContainsObjects( j ) )
    {
        return;
    }

    nlohmann::json jsonArray = j;
    j = nlohmann::json( );
    for ( unsigned int i = 0; i < jsonArray.size( ); ++i )
    {
        j[ "@" + std::to_string( i ) ] = jsonArray.at( i );
    }
}

//! Convert \p jsonObject to a json array.
nlohmann::json getAsArray( const nlohmann::json& jsonObject )
{
    nlohmann::json jsonArray = jsonObject;

    bool isObjectWithIntConvertibleKeys = jsonObject.is_object( );
    std::vector< unsigned int > indices;
    std::vector< nlohmann::json > values;
    if ( isObjectWithIntConvertibleKeys )
    {
        for ( nlohmann::json::const_iterator it = jsonObject.begin( ); it != jsonObject.end( ); ++it )
        {
            const std::string key = it.key( );
            if ( ! contains( SpecialKeys::all, key ) )
            {
                const int intKey = indexFromKey( key );
                if ( intKey >= 0 )
                {
                    indices.push_back( intKey );
                }
                else
                {
                    isObjectWithIntConvertibleKeys = false;
                    break;
                }
                values.push_back( getValue< nlohmann::json >( jsonObject, key ) );
            }
        }
    }

    if ( isObjectWithIntConvertibleKeys )
    {
        if ( values.empty( ) )
        {
            jsonArray = nlohmann::json( );
        }
        else
        {
            std::vector< nlohmann::json > vector( *std::max_element( indices.begin( ), indices.end( ) ) + 1 );
            for ( unsigned int i = 0; i < indices.size( ); ++i )
            {
                vector[ indices.at( i ) ] = values.at( i );
            }
            jsonArray = vector;
        }
    }

    return jsonArray;
}

//! Whether \p j is convertible to a json array.
bool isConvertibleToArray( const nlohmann::json& j )
{
    return getAsArray( j ).is_array( );
}

}  // namespace json_interface

}  // namespace tudat
