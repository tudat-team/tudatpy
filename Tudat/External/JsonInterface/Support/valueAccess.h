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

//! Access/modify a key of a `json` object or array.
/*!
 * Access/modify a key of a `json` object or array. \p jsonObject is passed by reference and the returned value is a
 * reference, so this method can be used to modify a \p jsonObject (e.g. `valueAt( jsonObject, myKey ) = newValue`).
 * Supports json arrays. If the field "key" does not exist, this function will try to convert it to integer and
 * access \p jsonObject at that index.
 * \param jsonObject The `json` object.
 * \param key The key to access.
 * \return A reference to the value of the accessed key.
 */
json& valueAt( json& jsonObject, const std::string& key );

//! Access a key of a `json` object or array.
/*!
 * Access a key of a `json` object or array. \p jsonObject is constant, so the returned value cannot be modified.
 * Supports json arrays. If the field "key" does not exist, this function will try to convert it to integer and
 * access \p jsonObject at that index.
 * \param jsonObject The constant `json` object.
 * \param key The key to access.
 * \return A constant reference to the value of the accessed key.
 */
const json& valueAt( const json& jsonObject, const std::string& key );

//! Access a key path of a `json` object or array.
/*!
 * Access a key path of a `json` object or array. \p jsonObject is constant, so the returned value cannot be modified.
 * Supports json arrays. If the any key in \p keyPath does not exist, this function will try to convert it to integer
 * and access \p jsonObject at that index.
 * \param jsonObject The constant `json` object.
 * \param keyPath The key path to access.
 * \return The value at the accessed key path.
 */
json valueAt( json jsonObject, const KeyPath& keyPath );

//! Whether the key at \p keyPath is defined for \p jsonObject.
/*!
 * @copybrief defined
 * \param jsonObject The `json` object.
 * \param keyPath The key path specifying to key to be checked.
 * \return @copybrief defined
 */
bool defined( const json& jsonObject, const KeyPath& keyPath );


// SPECIAL KEYS ACCESS

//! Get the \p jsonObject at key SpecialKeys::rootObject.
/*!
 * copybrief getRootObject
 * \param jsonObject The `json` object.
 * \return JSON representation of the root object of \p jsonObject.
 * \throws UndefinedKeyError If the key SpecialKeys::rootObject is not defined for \p jsonObject.
 */
json getRootObject( const json& jsonObject );

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

//! Get the response type to an event for a `json` object.
/*!
 * @copybrief getResponseToEventNamed
 * \param jsonObject The `json` object in which the type of response to the event is defined, or a `json` object with
 * a root object containing this information.
 * \param eventName The name of the event.
 * \param defaultResponse Default response if the response to the event is not defined. By default, continueSilently.
 * \return Response type to the event.
 */
ExceptionResponseType getResponseToEventNamed( const json& jsonObject, const std::string& eventName,
                                               const ExceptionResponseType defaultResponse = continueSilently );


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
    json rootObject;
    json currentObject = jsonObject;
    KeyPath currentKeyPath = keyPath;
    KeyPath canonicalKeyPath = keyPath;

    if ( ! contains( currentKeyPath, SpecialKeys::keyPath ) )
    {
        canonicalKeyPath = currentKeyPath.canonical( getKeyPath( currentObject ) );
        if ( ! contains( currentKeyPath, SpecialKeys::rootObject ) )
        {
            if ( defined( currentObject, SpecialKeys::rootObject ) )
            {
                rootObject = getValue< json >( currentObject, SpecialKeys::rootObject );
            }

            if ( currentKeyPath.size( ) > 0 )
            {
                if ( currentKeyPath.front( ) == SpecialKeys::root || contains( currentKeyPath, SpecialKeys::up ) )
                {
                    // Path is absolute or contains ..

                    // Use absolute key path
                    currentKeyPath = canonicalKeyPath;

                    // Use jsonObject's root object instead
                    if ( rootObject.is_null( ) )
                    {
                        throw UndefinedKeyError( canonicalKeyPath );
                    }
                    currentObject = rootObject;
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
        currentObject = valueAt( currentObject, currentKeyPath );
    }
    catch ( ... )
    {
        // At least one of the keys was missing
        throw UndefinedKeyError( canonicalKeyPath );
    }

    const bool chunkObjectWasArray = currentObject.is_array( );
    if ( ! containsAnyOf( currentKeyPath, SpecialKeys::all ) )
    {
        // If jsonObject is an array, convert to object
        convertToObjectIfArray( currentObject );

        if ( currentObject.is_object( ) )
        {
            // Define keys rootObject and keyPath of jsonObject to be returned
            currentObject[ SpecialKeys::rootObject ] = rootObject.is_null( ) ? jsonObject : rootObject;
            currentObject[ SpecialKeys::keyPath ] = canonicalKeyPath;
        }
    }

    const IllegalValueError illegalValueError( canonicalKeyPath, currentObject, typeid( ValueType ) );
    try
    {
        // Try to convert to the requested type
        const ValueType convertedObject = currentObject.get< ValueType >( );

        // Check if unidimensional array inference was applied
        json convertedJsonObject = json( convertedObject );
        if ( convertedJsonObject.is_array( ) && ! chunkObjectWasArray )
        {
            const ExceptionResponseType response =
                    getResponseToEventNamed( jsonObject, Keys::Options::unidimensionalArrayInference );
            if ( response == throwError )
            {
                throw illegalValueError;
            }
            else
            {
                if ( ! containsAnyOf( currentKeyPath, SpecialKeys::all ) )
                {
                    convertToObjectIfArray( convertedJsonObject );
                    convertedJsonObject[ SpecialKeys::rootObject ] = rootObject.is_null( ) ? jsonObject : rootObject;
                    convertedJsonObject[ SpecialKeys::keyPath ] = canonicalKeyPath;
                }

                if ( response == printWarning )
                {
                    std::cerr << "Unidimensional array inferred for key: " << canonicalKeyPath << std::endl;
                }

                return convertedJsonObject.get< ValueType >( );
            }
        }

        return convertedObject;
    }
    catch ( const nlohmann::detail::type_error& )
    {
        // Could not convert to the requested type
        throw illegalValueError;
    }
    catch ( const UnknownEnumError& )
    {
        // Could not convert string to enum
        throw illegalValueError;
    }
    catch ( const UndefinedKeyError& error )
    {
        // Some of the keys needed during creation of object of type ValueType were missing
        throw error;
    }
    catch ( const IllegalValueError& error )
    {
        // Some of the values for the keys needed during creation of object of type ValueType were illegal
        throw error;
    }
    /*
    catch ( ... )
    {
        throw UnrecognizedValueAccessError( canonicalKeyPath, typeid( ValueType ) );
    }
    */
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
    catch ( const IllegalValueError& error )
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
    catch ( const IllegalValueError& error )
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

//! Get the value of \p jsonObject at \p keyPath, or return \p optionalValue if not defined.
/*!
 * Get the value of \p jsonObject at \p keyPath, or return \p optionalValue if not defined.
 * \param jsonObject The `json` object.
 * \param keyPath Key path from which to retrieve the value.
 * \param defaultValue The default value to be returned if the key is not defined.
 * \return Value at the requested key path, or \p optionalValue if not defined.
 * \throws UndefinedKeyError If some of the subkeys needed to create the shared pointer of `ValueType` are missing
 * (only when \p jsonObject at \p keyPath is of type object).
 * \throws IllegalValueError<ValueType> If the obtained value for the requested key path is not convertible to
 * `ValueType`.
 * \throws IllegalValueError<SubvalueType> If some of the subkeys needed to create the shared pointer of `ValueType`
 * are not convertible to `SubvalueType` (only when \p jsonObject at \p keyPath is of type object).
 */
template< typename ValueType >
ValueType getValue( const json& jsonObject, const KeyPath& keyPath, const ValueType& defaultValue )
{
    try
    {
        return getValue< ValueType >( jsonObject, keyPath );
    }
    catch ( const UndefinedKeyError& error )
    {
        error.rethrowIfNotTriggeredByMissingValueAt( keyPath );
        error.rethrowIfDefaultValuesNotAllowed(
                    getResponseToEventNamed( jsonObject, Keys::Options::defaultValueUsedForMissingKey ),
                    json( defaultValue ) );
        return defaultValue;
    }
}

//! Get the numeric value of \p jsonObject at \p keyPath, or return \p optionalValue if not defined.
/*!
 * Get the numeric value of \p jsonObject at \p keyPath, trying to parse it as a physical magnitude with units and
 * convert to SI units if it is as string, or return \p optionalValue if not defined.
 * \param jsonObject The `json` object.
 * \param keyPath Key path from which to retrieve the value.
 * \param defaultValue The default value to be returned if the key is not defined.
 * \param allowNaN Whether the returned value is allowed to be NaN. Default is `false` (cannot be NaN).
 * \return Value at the requested key path, or \p optionalValue if not defined.
 * \throws UndefinedKeyError If the requested key path is not defined AND \p optionalValue is NaN AND \p allowNaN
 * is set to `false`.
 * \throws IllegalValueError<ValueType> If the obtained value for the requested key path is not convertible to
 * `ValueType`.
 */
template< typename NumberType >
NumberType getNumeric( const json& jsonObject, const KeyPath& keyPath,
                       const NumberType& defaultValue, bool allowNaN = false )
{
    try
    {
        return getNumeric< NumberType >( jsonObject, keyPath );
    }
    catch ( const UndefinedKeyError& error )
    {
        error.rethrowIfNotTriggeredByMissingValueAt( keyPath );
        error.rethrowIfDefaultValuesNotAllowed(
                    getResponseToEventNamed( jsonObject, Keys::Options::defaultValueUsedForMissingKey ),
                    json( defaultValue ) );
        error.rethrowIfNaNNotAllowed( allowNaN, defaultValue );
        return defaultValue;
    }
}

//! Get the numeric value of \p jsonObject at \p keyPath, or return \p optionalValue if not defined.
/*!
 * Get the numeric value of \p jsonObject at \p keyPath, trying to parse it as a date and convert to seconds since
 * J2000 if it is as string, or return \p optionalValue if not defined.
 * \param jsonObject The `json` object.
 * \param keyPath Key path from which to retrieve the value.
 * \param defaultValue The default value to be returned if the key is not defined.
 * \param allowNaN Whether the returned value is allowed to be NaN. Default is `false` (cannot be NaN).
 * \return Value at the requested key path, or \p optionalValue if not defined.
 * \throws UndefinedKeyError If the requested key path is not defined AND \p optionalValue is NaN AND \p allowNaN
 * is set to `false`.
 * \throws IllegalValueError<ValueType> If the obtained value for the requested key path is not convertible to
 * `ValueType`.
 */
template< typename NumberType >
NumberType getEpoch( const json& jsonObject, const KeyPath& keyPath,
                     const NumberType& defaultValue, bool allowNaN = false )
{
    try
    {
        return getEpoch< NumberType >( jsonObject, keyPath );
    }
    catch ( const UndefinedKeyError& error )
    {
        error.rethrowIfNotTriggeredByMissingValueAt( keyPath );
        error.rethrowIfDefaultValuesNotAllowed(
                    getResponseToEventNamed( jsonObject, Keys::Options::defaultValueUsedForMissingKey ),
                    json( defaultValue ) );
        error.rethrowIfNaNNotAllowed( allowNaN, defaultValue );
        return defaultValue;
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
void updateFromJSON( T& value, const json& jsonObject, const KeyPath& keyPath = { } )
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
void updateFromJSONIfDefined( T& value, const json& jsonObject, const KeyPath& keyPath = { } )
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
