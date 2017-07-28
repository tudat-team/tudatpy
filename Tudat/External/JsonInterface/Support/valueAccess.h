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

#include "keys.h"
#include "units.h"
#include "utilities.h"

namespace tudat
{

namespace json_interface
{

//! -DOC
class UnknownEnumError : public std::runtime_error
{
public:
    UnknownEnumError( ) : runtime_error( "Unknown enum." ) { }
};

//! -DOC
class JSONError : public std::runtime_error
{
public:
    //! Constructor.
    JSONError( const std::string& errorMessage, const KeyPath& keyPath )
        : runtime_error( errorMessage.c_str( ) ), keyPath( keyPath ) { }

    //! Error message.
    virtual const char* what( ) const throw( )
    {
        std::ostringstream stream;
        stream << runtime_error::what( ) << ": " << keyPath;
        std::cerr << stream.str( ).c_str( ) << std::endl;  // FIXME
        return stream.str( ).c_str( );
    }

    const KeyPath keyPath;
};

//! -DOC
class UndefinedKeyError : public JSONError
{
public:
    //! Constructor.
    UndefinedKeyError( const KeyPath& keyPath ) : JSONError( "Undefined key", keyPath ) { }
};

//! -DOC
class IllegalValueError : public JSONError
{
public:
    //! Associated (illegal) value.
    json value;

    //! Constructor.
    IllegalValueError( const KeyPath& keyPath, const json& value )
        : JSONError( "Illegal value for key", keyPath ), value( value ) { }

    //! Error message.
    virtual const char* what( ) const throw( )
    {
        std::ostringstream stream;
        stream << JSONError::what( ) << " = " << value;
        // std::cerr << stream.str( ).c_str( ) << std::endl;  // FIXME
        return stream.str( ).c_str( );
    }
};


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
ValueType getValue( json jsonObject, KeyPath keyPath )
{
    boost::shared_ptr< json > rootObject;
    json parentObject = jsonObject;

    KeyPath canonicalKeyPath = keyPath;
    if ( ! contains( keyPath, SpecialKeys::keyPath ) )
    {
        canonicalKeyPath = keyPath.canonical( getKeyPath( jsonObject ) );
        if ( ! contains( keyPath, SpecialKeys::rootObject ) )
        {
            rootObject = getRootObject( jsonObject );
            if ( keyPath.size( ) > 0 )
            {
                if ( *keyPath.begin( ) == SpecialKeys::root || contains( keyPath, SpecialKeys::up ) )
                {
                    // Path is absolute or contains ..

                    // Use absolute key path
                    keyPath = canonicalKeyPath;

                    // Use jsonObject's root object instead
                    if ( ! rootObject )
                    {
                        std::cerr << "Could not use absolute key path because "
                                  << "the JSON object does not contain a root object." << std::endl;
                        throw;
                    }
                    jsonObject = *rootObject;
                }
            }
        }
    }

    if ( keyPath.isAbsolute( ) )
    {
        // Remove "~"
        keyPath.erase( keyPath.begin( ) );
    }

    try
    {
        // Recursively update jsonObject for every key in keyPath
        for ( const std::string key : keyPath )
        {
            try
            {
                // Try to access element at key
                jsonObject = jsonObject.at( key );
            }
            catch ( ... )
            {
                // Key may be convertible to int.
                jsonObject = jsonObject.at( std::stoi( key ) );
            }
        }
    }
    catch ( ... )
    {
        // At least one of the keys was missing
        throw UndefinedKeyError( canonicalKeyPath );
    }

    if ( ! contains( keyPath, SpecialKeys::all ) )
    {
        // If jsonObject is an array, convert to object
        if ( jsonObject.is_array( ) )
        {
            json jsonArray = jsonObject;
            jsonObject = json( );
            for ( unsigned int i = 0; i < jsonArray.size( ); ++i )
            {
                jsonObject[ std::to_string( i ) ] = jsonArray.at( i );
            }
        }

        // Define keys rootObject and keyPath of jsonObject to be returned
        if ( jsonObject.is_object( ) )
        {
            jsonObject[ SpecialKeys::rootObject ] = rootObject ? *rootObject : parentObject;
            jsonObject[ SpecialKeys::keyPath ] = canonicalKeyPath;
        }
    }

    try
    {
        // Try to convert to the requested type
        return jsonObject.get< ValueType >( );
    }
    catch ( const nlohmann::detail::type_error& )
    {
        // Could not convert to the requested type
        throw IllegalValueError( canonicalKeyPath, jsonObject );
    }
    catch ( const UnknownEnumError& )
    {
        // Could not convert string to enum
        throw IllegalValueError( canonicalKeyPath, jsonObject );
    }
}


//! -DOC
template< typename NumberType >
NumberType getNumeric( json jsonObject, const KeyPath& keyPath )
{
    try
    {
        return getValue< NumberType >( jsonObject, keyPath );
    }
    catch ( const IllegalValueError& error )
    {
        // Could not convert to the requested type
        try
        {
            // Convert to string (with units) and then parse the number and convert to SI
            return parseMagnitudeWithUnits< NumberType >( error.value.get< std::string >( ) );
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
NumberType getEpoch( json jsonObject, const KeyPath& keyPath )
{
    try
    {
        return getValue< NumberType >( jsonObject, keyPath );
    }
    catch ( const IllegalValueError& error )
    {
        // Could not convert to the requested type
        try
        {
            // Convert to string and then parse as a formatted date
            return convertToSecondsSinceJ2000< NumberType >( error.value.get< std::string >( ) );
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
boost::shared_ptr< ValueType > getValuePointer( const json& jsonObject, const KeyPath& keyPath )
{
    try
    {
        return boost::make_shared< ValueType >( getValue< ValueType >( jsonObject, keyPath ) );
    }
    catch ( const UndefinedKeyError& error )
    {
        return NULL;
    }
}


//! -DOC
template< typename NumberType >
boost::shared_ptr< NumberType > getNumericPointer( const json& jsonObject, const KeyPath& keyPath )
{
    try
    {
        return boost::make_shared< NumberType >( getNumeric< NumberType >( jsonObject, keyPath ) );
    }
    catch ( const UndefinedKeyError& error )
    {
        return NULL;
    }
}


//! -DOC
template< typename NumberType >
boost::shared_ptr< NumberType > getEpochPointer( const json& jsonObject, const KeyPath& keyPath )
{
    try
    {
        return boost::make_shared< NumberType >( getEpoch< NumberType >( jsonObject, keyPath ) );
    }
    catch ( const UndefinedKeyError& error )
    {
        return NULL;
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


//! -DOC
template< typename EnumType >
EnumType enumFromString( const std::string& stringValue,
                         const std::map< std::string, EnumType >& possibleValues )
{
    try
    {
        return possibleValues.at( stringValue );
    }
    catch ( ... )
    {
        std::cerr << "Unsupported string \"" << stringValue << "\" for enum " <<
                     boost::core::demangled_name( typeid( EnumType ) ) << std::endl;
        throw UnknownEnumError( );
    }
}

//! -DOC
template< typename EnumType >
std::string stringFromEnum( const EnumType enumValue, const std::map< std::string, EnumType >& possibleValues )
{
    for ( auto ent : possibleValues )
    {
        if ( ent.second == enumValue )
        {
            return ent.first;
        }
    }
    std::cerr << "Unknown string representation for enum value "
              << boost::core::demangled_name( typeid( EnumType ) ) << "::" << enumValue << std::endl;
    throw UnknownEnumError( );
}


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
