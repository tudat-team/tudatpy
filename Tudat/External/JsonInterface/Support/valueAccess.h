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

// KEY ACCESS

//! Whether the key at \p keyPath is defined for \p jsonObject.
/*!
 * @copybrief defined
 * \param jsonObject The `json` object.
 * \param keyPath The key path specifying to key to be checked.
 * \return @copybrief defined
 */
bool defined( const json& jsonObject, const KeyPath& keyPath );

//! Get the a shared pointer to \p jsonObject at key SpecialKeys::rootObject.
/*!
 * copybrief getRootObject
 * \param jsonObject The `json` object.
 * \return Sshared pointer to \p jsonObject at key SpecialKeys::rootObject.
 * `NULL` if key SpecialKeys::rootObject is not defined for \p jsonObject.
 */
boost::shared_ptr< json > getRootObject( const json& jsonObject );

//! Get the absolute key path from which \p jsonObject was retrieved.
/*!
 * copybrief getKeyPath
 * \param jsonObject The `json` object.
 * \return Absolute key path from which \p jsonObject was retrieved.
 * \throws UndefinedKeyError If the key SpecialKeys::keyPath is not defined for \p jsonObject.
 */
KeyPath getKeyPath( const json& jsonObject );

//! Get the key at which \p jsonObject was obtained.
/*!
 * copybrief getParentKey
 * \param jsonObject The `json` object.
 * \param errorMessage Error message to be printed if \p jsonObject has no key SpecialKeys::keyPath.
 * \return Key at which \p jsonObject was obtained.
 * \throws UndefinedKeyError If the key SpecialKeys::keyPath is not defined for \p jsonObject.
 */
std::string getParentKey( const json& jsonObject,
                          const std::string& errorMessage = "Could not determine parent key: context is missing." );

//! Convert \p j to object if \p j is array.
/*!
 * @copybrief convertToObjectIfArray This method does nothing if \p j is not an array, if \p is empty, or if \p is
 * an array of primitive elements and \p onlyIfElementsAreStructured is set to `true`.
 * \param j `json` object to be converted.
 * \param onlyIfElementsAreStructured Only perform the conversion if the elements of \p j are either an object or an
 * array, i.e. if they are not primitive. Defualt is `false` (all arrays will be converted regardless of their elements
 * type).
 */
void convertToObjectIfArray( json& j, const bool onlyIfElementsAreStructured = false );


// ACCESS HISTORY

//! Global variable containing all the key paths that were accessed since clearAccessHistory() was called for the
//! last time (or since this variable was initialized).
extern std::set< KeyPath > accessedKeyPaths;

//! Clear the global variable accessedKeyPaths.
/*!
 * @copybrief clearAccessHistory
 */
inline void clearAccessHistory( )
{
    accessedKeyPaths.clear( );
}

//! Check for key paths that are defined in \p jsonObject but not contained by the global variable accessedKeyPaths.
/*!
 * @copybrief checkUnusedKeys
 * \param jsonObject The `json` object.
 * \param response Response type when finding unused keys in \p jsonObject.
 * \throws std::runtime_error If \p response is set to ExceptionResponseType::throwError and at least one unsued key
 * was found in \p jsonObject.
 */
void checkUnusedKeys( const json& jsonObject, const ExceptionResponseType response );


// GET FROM JSON

//! Get the value of \p jsonObject at \p keyPath.
/*!
 * Get the value of \p jsonObject at \p keyPath.
 * \param jsonObject The `json` object.
 * \param keyPath Key path from which to retrieve the value.
 * \return Value at the requested key path.
 * \throws UndefinedKeyError If the requested key path is not defined.
 * \throws UndefinedKeyError If some of the subkeys needed to create the shared pointer of `ValueType` are missing
 * (only when \p jsonObject at \p keyPath is of type object).
 * \throws IllegalValueError<ValueType> If the obtained value for the requested key path is not convertible to
 * `ValueType`.
 * \throws IllegalValueError<SubvalueType> If some of the subkeys needed to create the shared pointer of `ValueType`
 * are not convertible to `SubvalueType` (only when \p jsonObject at \p keyPath is of type object).
 */
template< typename ValueType >
ValueType getValue( const json& jsonObject, const KeyPath& keyPath )
{
    boost::shared_ptr< json > rootObjectPointer;
    json currentObject = jsonObject;
    KeyPath currentKeyPath = keyPath;
    KeyPath canonicalKeyPath = keyPath;

    if ( ! contains( currentKeyPath, SpecialKeys::keyPath ) )
    {
        canonicalKeyPath = currentKeyPath.canonical( getKeyPath( currentObject ) );
        if ( ! contains( currentKeyPath, SpecialKeys::rootObject ) )
        {
            rootObjectPointer = getRootObject( currentObject );
            if ( currentKeyPath.size( ) > 0 )
            {
                if ( *currentKeyPath.begin( ) == SpecialKeys::root || contains( currentKeyPath, SpecialKeys::up ) )
                {
                    // Path is absolute or contains ..

                    // Use absolute key path
                    currentKeyPath = canonicalKeyPath;

                    // Use jsonObject's root object instead
                    if ( ! rootObjectPointer )
                    {
                        std::cerr << "Could not use absolute key path because "
                                  << "the JSON object does not contain a root object." << std::endl;
                        throw;
                    }
                    currentObject = *rootObjectPointer;
                }
            }
        }
    }

    if ( currentKeyPath.isAbsolute( ) )
    {
        // Remove "~"
        currentKeyPath.erase( currentKeyPath.begin( ) );
    }

    accessedKeyPaths.insert( canonicalKeyPath );

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

    if ( ! containsAnyOf( currentKeyPath, SpecialKeys::all ) )
    {
        // If jsonObject is an array, convert to object
        convertToObjectIfArray( currentObject );

        if ( currentObject.is_object( ) )
        {
            // Define keys rootObject and keyPath of jsonObject to be returned
            currentObject[ SpecialKeys::rootObject ] = rootObjectPointer ? *rootObjectPointer : jsonObject;
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
}

//! Convert \p jsonObject to `ValueType`.
/*!
 * Convert \p jsonObject to `ValueType`.
 * \param jsonObject The `json` object.
 * \return \p jsonObject as `ValueType`.
 * \throws IllegalValueError<ValueType> If \p jsonObject is not convertible to `ValueType`.
 * \throws IllegalValueError<SubvalueType> If some of the subkeys needed to create the shared pointer of `ValueType`
 * are not convertible to `SubvalueType` (only when \p jsonObject at \p keyPath is of type object).
 */
template< typename ValueType >
ValueType getAs( const json& jsonObject )
{
    return getValue< ValueType >( jsonObject, { } );
}

//! Get the numeric value of \p jsonObject at \p keyPath.
/*!
 * Get the numeric value of \p jsonObject at \p keyPath, trying to parse it as a physical magnitude with units and
 * convert to SI units if it is as string.
 * \param jsonObject The `json` object.
 * \param keyPath Key path from which to retrieve the value.
 * \return Numeric value at the requested key path.
 * \throws UndefinedKeyError If the requested key path is not defined.
 * \throws IllegalValueError<NumberType> If the obtained value for the requested key path is not convertible to
 * `NumberType`.
 */
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

//! Get the numeric value of \p jsonObject at \p keyPath.
/*!
 * Get the numeric value of \p jsonObject at \p keyPath, trying to parse it as a date and convert to seconds since
 * J2000 if it is as string.
 * \param jsonObject The `json` object.
 * \param keyPath Key path from which to retrieve the value.
 * \return Numeric value at the requested key path.
 * \throws UndefinedKeyError If the requested key path is not defined.
 * \throws IllegalValueError<NumberType> If the obtained value for the requested key path is not convertible to
 * `NumberType`.
 */
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

//! Get a pointer to the value of \p jsonObject at \p keyPath.
/*!
 * Get a pointer to the value of \p jsonObject at \p keyPath, or `NULL` if the key path does not exist.
 * \param jsonObject The `json` object.
 * \param keyPath Key path from which to retrieve the value.
 * \param getFunction Function used to retrieve the value. Default function is getValue. Can be set to
 * getNumeric or getEpoch to enable parsing of the obtained value when provided as a string.
 * \return Value at the requested key path.
 * \throws UndefinedKeyError If some of the subkeys needed to create the shared pointer of `ValueType` are missing
 * (only when \p jsonObject at \p keyPath is of type object).
 * \throws IllegalValueError<ValueType> If the obtained value for the requested key path is not convertible to
 * `ValueType`.
 */
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
        error.rethrowIfNotTriggeredByMissingValueAt( keyPath );
        return NULL;
    }
}

//! Get a pointer to the numeric value of \p jsonObject at \p keyPath.
/*!
 * Get a pointer to the numeric value of \p jsonObject at \p keyPath, trying to parse it as a physical magnitude with
 * units and convert to SI units if it is as string, or `NULL` if the key path does not exist.
 * \param jsonObject The `json` object.
 * \param keyPath Key path from which to retrieve the value.
 * \return Value at the requested key path.
 * \throws IllegalValueError<NumberType> If the obtained value for the requested key path is not convertible to
 * `NumberType`.
 */
template< typename NumberType >
boost::shared_ptr< NumberType > getOptionalNumeric( const json& jsonObject, const KeyPath& keyPath )
{
    return getOptional( jsonObject, keyPath, getNumeric< NumberType > );
}

//! Get a pointer to the numeric value of \p jsonObject at \p keyPath.
/*!
 * Get a pointer to the numeric value of \p jsonObject at \p keyPath, trying to parse it as a date and convert to
 * seconds since J2000 if it is as string, or `NULL` if the key path does not exist.
 * \param jsonObject The `json` object.
 * \param keyPath Key path from which to retrieve the value.
 * \return Value at the requested key path.
 * \throws IllegalValueError<NumberType> If the obtained value for the requested key path is not convertible to
 * `NumberType`.
 */
template< typename NumberType >
boost::shared_ptr< NumberType > getOptionalEpoch( const json& jsonObject, const KeyPath& keyPath )
{
    return getOptional( jsonObject, keyPath, getEpoch< NumberType > );
}

//! Get the value of \p jsonObject at \p keyPath, or return \p optionalValue if not defined.
/*!
 * Get the value of \p jsonObject at \p keyPath, or return \p optionalValue if not defined.
 * \param jsonObject The `json` object.
 * \param keyPath Key path from which to retrieve the value.
 * \param optionalValue The default value to be returned if the key is not defined.
 * \return Value at the requested key path, or \p optionalValue if not defined.
 * \throws UndefinedKeyError If some of the subkeys needed to create the shared pointer of `ValueType` are missing
 * (only when \p jsonObject at \p keyPath is of type object).
 * \throws IllegalValueError<ValueType> If the obtained value for the requested key path is not convertible to
 * `ValueType`.
 * \throws IllegalValueError<SubvalueType> If some of the subkeys needed to create the shared pointer of `ValueType`
 * are not convertible to `SubvalueType` (only when \p jsonObject at \p keyPath is of type object).
 */
template< typename ValueType >
ValueType getValue( const json& jsonObject, const KeyPath& keyPath, const ValueType& optionalValue )
{
    try
    {
        return getValue< ValueType >( jsonObject, keyPath );
    }
    catch ( const UndefinedKeyError& error )
    {
        error.rethrowIfNotTriggeredByMissingValueAt( keyPath );
        return optionalValue;
    }
}

//! Get the numeric value of \p jsonObject at \p keyPath, or return \p optionalValue if not defined.
/*!
 * Get the numeric value of \p jsonObject at \p keyPath, trying to parse it as a physical magnitude with units and
 * convert to SI units if it is as string, or return \p optionalValue if not defined.
 * \param jsonObject The `json` object.
 * \param keyPath Key path from which to retrieve the value.
 * \param optionalValue The default value to be returned if the key is not defined.
 * \param allowNaN Whether the returned value is allowed to be NaN. Default is `false` (cannot be NaN).
 * \return Value at the requested key path, or \p optionalValue if not defined.
 * \throws UndefinedKeyError If the requested key path is not defined AND \p optionalValue is NaN AND \p allowNaN
 * is set to `false`.
 * \throws IllegalValueError<ValueType> If the obtained value for the requested key path is not convertible to
 * `ValueType`.
 */
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
        error.rethrowIfNotTriggeredByMissingValueAt( keyPath );  // not needed? NumberType cannot have subkeys
        if ( ! allowNaN && isNaN( optionalValue ) )
        {
            throw error;
        }
        else
        {
            return optionalValue;
        }
    }
}

//! Get the numeric value of \p jsonObject at \p keyPath, or return \p optionalValue if not defined.
/*!
 * Get the numeric value of \p jsonObject at \p keyPath, trying to parse it as a date and convert to seconds since
 * J2000 if it is as string, or return \p optionalValue if not defined.
 * \param jsonObject The `json` object.
 * \param keyPath Key path from which to retrieve the value.
 * \param optionalValue The default value to be returned if the key is not defined.
 * \param allowNaN Whether the returned value is allowed to be NaN. Default is `false` (cannot be NaN).
 * \return Value at the requested key path, or \p optionalValue if not defined.
 * \throws UndefinedKeyError If the requested key path is not defined AND \p optionalValue is NaN AND \p allowNaN
 * is set to `false`.
 * \throws IllegalValueError<ValueType> If the obtained value for the requested key path is not convertible to
 * `ValueType`.
 */
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
        error.rethrowIfNotTriggeredByMissingValueAt( keyPath );  // not needed? NumberType cannot have subkeys
        if ( ! allowNaN && isNaN( optionalValue ) )
        {
            throw error;
        }
        else
        {
            return optionalValue;
        }
    }
}


// UPDATE FROM JSON

//! Set \p value to be the value of \p jsonObject at \p keyPath.
/*!
 * Set \p value to be the value of \p jsonObject at \p keyPath.
 * \param value The value to be updated (i.e. re-initialized).
 * \param jsonObject The `json` object.
 * \param keyPath Key path from which to retrieve the value.
 * \throws UndefinedKeyError If the requested key path is not defined.
 * \throws UndefinedKeyError If some of the subkeys needed to create the shared pointer of `ValueType` are missing
 * (only when \p jsonObject at \p keyPath is of type object).
 * \throws IllegalValueError<ValueType> If the obtained value for the requested key path is not convertible to
 * `ValueType`.
 * \throws IllegalValueError<SubvalueType> If some of the subkeys needed to create the shared pointer of `ValueType`
 * are not convertible to `SubvalueType` (only when \p jsonObject at \p keyPath is of type object).
 */
template< typename T >
void updateFromJSON( T& value, const json& jsonObject, const KeyPath& keyPath )
{
    value = getValue< T >( jsonObject, keyPath );
}

//! Set \p value to be the value of \p jsonObject at \p keyPath if defined.
/*!
 * Set \p value to be the value of \p jsonObject at \p keyPath if defined.
 * \remark This function does nothing if \p jsonObject is not defined at \p keyPath.
 * \param value The value to be updated (i.e. re-initialized).
 * \param jsonObject The `json` object.
 * \param keyPath Key path from which to retrieve the value.
 * \throws UndefinedKeyError If some of the subkeys needed to create the shared pointer of `ValueType` are missing
 * (only when \p jsonObject at \p keyPath is of type object).
 * \throws IllegalValueError<ValueType> If the obtained value for the requested key path is not convertible to
 * `ValueType`.
 * \throws IllegalValueError<SubvalueType> If some of the subkeys needed to create the shared pointer of `ValueType`
 * are not convertible to `SubvalueType` (only when \p jsonObject at \p keyPath is of type object).
 */
template< typename T >
void updateFromJSONIfDefined( T& value, const json& jsonObject, const KeyPath& keyPath )
{
    try
    {
        updateFromJSON( value, jsonObject, keyPath );
    }
    catch ( const UndefinedKeyError& error )
    {
        error.rethrowIfNotTriggeredByMissingValueAt( keyPath );
    }
}


// UPDATE JSON

//! Assign \p value to `jsonObject[ key ]` if \p value is not NaN.
/*!
 * @copybrief assignIfNotNaN
 * \remark This function does nothing if \p value is NaN.
 * \remark The comparison operator must be defined for `EquatableType`.
 * \param jsonObject The `json` object being updated.
 * \param key The key of \p jsonObject being updated.
 * \param value The value to be used for updating `jsonObject[ key ]`.
 */
template< typename EquatableType >
void assignIfNotNaN( json& jsonObject, const std::string& key, const EquatableType& value )
{
    if ( ! isNaN( value ) )
    {
        jsonObject[ key ] = value;
    }
}

//! Assign \p object to `jsonObject[ key ]` if \p object is not NULL.
/*!
 * @copybrief assignIfNotNull
 * \remark This function does nothing if \p object is `NULL`.
 * \param jsonObject The `json` object being updated.
 * \param key The key of \p jsonObject being updated.
 * \param object Shared pointer to the object that is being to be used to update `jsonObject[ key ]`.
 */
template< typename T >
void assignIfNotNull( json& jsonObject, const std::string& key, const boost::shared_ptr< T >& object )
{
    if ( object )
    {
        jsonObject[ key ] = object;
    }
}

//! Assign \p object to `jsonObject[ key ]` if \p container is not empty.
/*!
 * @copybrief assignIfNotEmpty
 * \remark This function does nothing if \p container is empty.
 * \remark The method `empty()` must be defined for `ContainerType`.
 * \param jsonObject The `json` object being updated.
 * \param key The key of \p jsonObject being updated.
 * \param container Container (e.g. `std::vector`, `std::set`, `std::map` or `std::unordered_map`) that is being to be
 * used to update `jsonObject[ key ]`.
 */
template< typename ContainerType >
void assignIfNotEmpty( json& jsonObject, const std::string& key, const ContainerType& container )
{
    if ( ! container.empty( ) )
    {
        jsonObject[ key ] = container;
    }
}


// Typedefs for single-body and body-to-body maps.

template < typename T >
using SingleBodyMap = std::unordered_map< std::string, std::vector< boost::shared_ptr< T > > >;

template < typename T >
using BodyToBodyMap = std::unordered_map< std::string, SingleBodyMap< T > >;

template < typename T >
using NoPointerSingleBodyMap = std::unordered_map< std::string, std::vector< T > >;

template < typename T >
using NoPointerBodyToBodyMap = std::unordered_map< std::string, NoPointerSingleBodyMap< T > >;

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_VALUEACCESS_H
