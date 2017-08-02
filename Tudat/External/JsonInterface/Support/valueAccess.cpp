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

#include "validation.h"

namespace tudat
{

namespace json_interface
{

//! -DOC
bool defined( const json& jsonObject, const KeyPath& keyPath )
{
    if ( getOptional< json >( jsonObject, keyPath ) )
    {
        return true;
    }
    return false;
}

//! -DOC
boost::shared_ptr< json > getRootObject( const json& jsonObject )
{
    return getOptional< json >( jsonObject, SpecialKeys::rootObject );
}

//! -DOC
KeyPath getKeyPath( const json& jsonObject )
{
    return getValue< KeyPath >( jsonObject, SpecialKeys::keyPath, SpecialKeys::root );
}

//! -DOC
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

//! -DOC
void convertToObjectIfArray( json& jsonObject, const bool onlyIfElementsAreStructured )
{
    if ( ! jsonObject.is_array( ) )
    {
        return;
    }
    if ( onlyIfElementsAreStructured )
    {
        if ( jsonObject.empty( ) )
        {
            return;
        }
        if ( ! jsonObject.front( ).is_structured( ) )
        {
            return;
        }
    }
    json jsonArray = jsonObject;
    jsonObject = json( );
    for ( unsigned int i = 0; i < jsonArray.size( ); ++i )
    {
        jsonObject[ std::to_string( i ) ] = jsonArray.at( i );
    }
}


//! -DOC
std::set< KeyPath > accessedKeyPaths = { };

//! -DOC
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

//! -DOC
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
