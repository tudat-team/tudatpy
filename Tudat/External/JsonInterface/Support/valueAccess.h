/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_VALUEACCESS_H
#define TUDAT_JSONINTERFACE_VALUEACCESS_H

#include <unordered_map>
#include <set>
#include <functional>

#include <boost/core/demangle.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "json/src/json.hpp"
using json = nlohmann::json;

#include <Tudat/InputOutput/basicInputOutput.h>

#include "errorHandling.h"
#include "keys.h"
#include "units.h"
#include "utilities.h"

namespace tudat
{

namespace json_interface
{

//! -DOC
boost::shared_ptr< json > getRootObject( const json& jsonObject );

//! -DOC
KeyPath getKeyPath( const json& jsonObject );

//! -DOC
std::string getParentKey( const json& jsonObject,
                          const std::string& errorMessage = "Could not determine parent key: context is missing." );


//! Get the value of a parameter defined by a `KeyPath` from a `json` object.
/*!
 * Get the value of a parameter defined by a `KeyPath` from a `json` object.
 * An error will be thrown if the requested key does not exist, or the provided value is not of the expected type `T`.
 * \param jsonObject JSON object from which to get the value.
 * \param keyPath Vector of keys defining the value to be accessed (key.subkey.subsubkey ...).
 * \return Value of the requested key.
 * \throw UndefinedKeyError If the requested key is not defined.
 * \throw IllegalValueError If the provided value for the requested key is not of type `T`.
 */
template< typename ValueType >
ValueType getValue( const json& jsonObject, const KeyPath& keyPath )
{
    boost::shared_ptr< json > rootObject;
    json currentObject = jsonObject;
    KeyPath currentKeyPath = keyPath;
    KeyPath canonicalKeyPath = keyPath;

    if ( ! contains( currentKeyPath, SpecialKeys::keyPath ) )
    {
        canonicalKeyPath = currentKeyPath.canonical( getKeyPath( currentObject ) );
        if ( ! contains( currentKeyPath, SpecialKeys::rootObject ) )
        {
            rootObject = getRootObject( currentObject );
            if ( currentKeyPath.size( ) > 0 )
            {
                if ( *currentKeyPath.begin( ) == SpecialKeys::root || contains( currentKeyPath, SpecialKeys::up ) )
                {
                    // Path is absolute or contains ..

                    // Use absolute key path
                    currentKeyPath = canonicalKeyPath;

                    // Use jsonObject's root object instead
                    if ( ! rootObject )
                    {
                        std::cerr << "Could not use absolute key path because "
                                  << "the JSON object does not contain a root object." << std::endl;
                        throw;
                    }
                    currentObject = *rootObject;
                }
            }
        }
    }

    if ( currentKeyPath.isAbsolute( ) )
    {
        // Remove "~"
        currentKeyPath.erase( currentKeyPath.begin( ) );
    }

    try
    {
        // Recursively update jsonObject for every key in keyPath
        for ( const std::string key : currentKeyPath )
        {
            try
            {
                // Try to access element at key
                currentObject = currentObject.at( key );
            }
            catch ( ... )
            {
                // Key may be convertible to int.
                currentObject = currentObject.at( std::stoi( key ) );
            }
        }
    }
    catch ( ... )
    {
        // At least one of the keys was missing
        throw UndefinedKeyError( canonicalKeyPath );
    }

    if ( ! contains( currentKeyPath, SpecialKeys::all ) )
    {
        // If jsonObject is an array, convert to object
        if ( currentObject.is_array( ) )
        {
            json jsonArray = currentObject;
            currentObject = json( );
            for ( unsigned int i = 0; i < jsonArray.size( ); ++i )
            {
                currentObject[ std::to_string( i ) ] = jsonArray.at( i );
            }
        }

        // Define keys rootObject and keyPath of jsonObject to be returned
        if ( currentObject.is_object( ) )
        {
            currentObject[ SpecialKeys::rootObject ] = rootObject ? *rootObject : jsonObject;
            currentObject[ SpecialKeys::keyPath ] = canonicalKeyPath;
        }
    }

    try
    {
        // Try to convert to the requested type
        return currentObject.get< ValueType >( );
    }
    catch ( const nlohmann::detail::type_error& )
    {
        // Could not convert to the requested type
        throw IllegalValueError< ValueType >( canonicalKeyPath, currentObject );
    }
    catch ( const UnknownEnumError& )
    {
        // Could not convert string to enum
        throw IllegalValueError< ValueType >( canonicalKeyPath, currentObject );
    }
    catch ( const UndefinedKeyError& error )
    {
        // Some of the keys that had to be defined for jsonObject are missing.
        if ( currentObject.is_object( ) )
        {
            // If the provided jsonObject is an object, the user wants to know which key is missing.
            // Thus, we re-throw the undefined key error.
            throw error;
        }
        else
        {
            // If the provided jsonObject is not an object (e.g. is a number or a string)
            // the user wants to know that an illegal value was provided for the expected object.
            // Thus, we throw an illegal value error.
            throw IllegalValueError< ValueType >( canonicalKeyPath, currentObject );
        }
    }
}


//! -DOC
template< typename NumberType >
NumberType getNumeric( const json& jsonObject, const KeyPath& keyPath )
{
    try
    {
        return getValue< NumberType >( jsonObject, keyPath );
    }
    catch ( const IllegalValueError< NumberType >& error )
    {
        // Could not convert to the requested type
        try
        {
            // Convert to string (with units) and then parse the number and convert to SI
            return parseMagnitudeWithUnits< NumberType >( error.value.template get< std::string >( ) );
        }
        catch ( ... )
        {
            // Could not convert to the requested type nor parse the string as a number
            throw error;
        }
    }
}


//! -DOC
template< typename NumberType >
NumberType getEpoch( const json& jsonObject, const KeyPath& keyPath )
{
    try
    {
        return getValue< NumberType >( jsonObject, keyPath );
    }
    catch ( const IllegalValueError< NumberType >& error )
    {
        // Could not convert to the requested type
        try
        {
            // Convert to string and then parse as a formatted date
            return convertToSecondsSinceJ2000< NumberType >( error.value.template get< std::string >( ) );
        }
        catch ( ... )
        {
            // Could not convert to the requested type nor parse the string as a date
            throw error;
        }
    }
}


//! -DOC
template< typename ValueType >
boost::shared_ptr< ValueType > getOptional(
        const json& jsonObject, const KeyPath& keyPath,
        const std::function< ValueType( const json&, const KeyPath& ) > getFunction = getValue< ValueType > )
{
    try
    {
        return boost::make_shared< ValueType >( getFunction( jsonObject, keyPath ) );
    }
    catch ( const UndefinedKeyError& error )
    {
        return NULL;
    }
}


//! -DOC
template< typename NumberType >
boost::shared_ptr< NumberType > getOptionalNumeric( const json& jsonObject, const KeyPath& keyPath )
{
    return getOptional( jsonObject, keyPath, getNumeric< NumberType > );
}


//! -DOC
template< typename NumberType >
boost::shared_ptr< NumberType > getOptionalEpoch( const json& jsonObject, const KeyPath& keyPath )
{
    return getOptional( jsonObject, keyPath, getEpoch< NumberType > );
}


//! -DOC
template< typename ValueType >
ValueType getValue( const json& jsonObject, const KeyPath& keyPath, const ValueType& optionalValue )
{
    try
    {
        return getValue< ValueType >( jsonObject, keyPath );
    }
    catch ( const UndefinedKeyError& error )
    {
        return optionalValue;
    }
}


//! -DOC
template< typename NumberType >
NumberType getNumeric( const json& jsonObject, const KeyPath& keyPath,
                       const NumberType& optionalValue, bool allowNaN = false )
{
    try
    {
        return getNumeric< NumberType >( jsonObject, keyPath );
    }
    catch ( const UndefinedKeyError& error )
    {
        if ( ! allowNaN && isnan( optionalValue ) )
        {
            throw error;
        }
        else
        {
            return optionalValue;
        }
    }
}


//! -DOC
template< typename NumberType >
NumberType getEpoch( const json& jsonObject, const KeyPath& keyPath,
                     const NumberType& optionalValue, bool allowNaN = false )
{
    try
    {
        return getEpoch< NumberType >( jsonObject, keyPath );
    }
    catch ( const UndefinedKeyError& error )
    {
        if ( ! allowNaN && isnan( optionalValue ) )
        {
            throw error;
        }
        else
        {
            return optionalValue;
        }
    }
}


//! -DOC
template< typename ValueType >
void updateFromJSON( ValueType& value, const json& jsonObject, const KeyPath& keyPath, bool throwIfNotDefined = true )
{
    try
    {
        value = getValue< ValueType >( jsonObject, keyPath );
    }
    catch ( const UndefinedKeyError& error )
    {
        if ( throwIfNotDefined || ( error.keyPath.back( ) != keyPath.back( ) ) )
        {
            throw error;
        }
    }
}

///! -DOC
bool defined( const json& jsonObject, const KeyPath& keyPath );


//! Support for single-body and body-to-body maps

template < typename T >
using SingleBodyMap = std::unordered_map< std::string, std::vector< boost::shared_ptr< T > > >;

template < typename T >
using BodyToBodyMap = std::unordered_map< std::string, SingleBodyMap< T > >;

template < typename T >
using NoPointerSingleBodyMap = std::unordered_map< std::string, std::vector< T > >;

template < typename T >
using NoPointerBodyToBodyMap = std::unordered_map< std::string, NoPointerSingleBodyMap< T > >;

/*
//! Get a `BodyToBodyMap` of `T` objects from a `json` object.
template < typename T >
BodyToBodyMap< T > getBodyToBodyMap(
        const json& jsonObject, const KeyPath& keyPath,
        std::function< boost::shared_ptr< T >( const json&, const KeyPath& ) > createFunction )
{
    BodyToBodyMap< T > bodyToBodyMap;

    NoPointerBodyToBodyMap< json > jsonBodyToBodyMap =
            getValue< NoPointerBodyToBodyMap< json > >( jsonObject, keyPath, { } );
    for ( auto entry : jsonBodyToBodyMap )
    {
        const std::string bodyUndergoing = entry.first;
        const NoPointerSingleBodyMap< json > jsonSingleBodyMap = entry.second;
        for ( auto subentry : jsonSingleBodyMap )
        {
            const std::string bodyExerting = subentry.first;
            const std::vector< json > jsonVector = subentry.second;
            for ( unsigned int i = 0; i < jsonVector.size( ); ++i )
            {
                bodyToBodyMap[ bodyUndergoing ][ bodyExerting ].push_back(
                            createFunction( jsonObject, keyPath / bodyUndergoing / bodyExerting / i ) );
            }
        }
    }

    return bodyToBodyMap;
}
*/

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_VALUEACCESS_H
