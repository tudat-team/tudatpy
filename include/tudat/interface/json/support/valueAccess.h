/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include <memory>
#include <boost/make_shared.hpp>

#include "tudat/io/basicInputOutput.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/interface/json/support/errorHandling.h"
#include "tudat/interface/json/support/keys.h"
#include "tudat/interface/json/support/utilities.h"

namespace tudat
{

namespace json_interface
{

// KEY ACCESS

//! Access/mutate a key of a `json` object or array.
/*!
 * Access/mutate a key of a `json` object or array. \p jsonObject is passed by reference and the returned value is a
 * reference, so this method can be used to mutate a \p jsonObject (e.g. `valueAt( jsonObject, myKey, true ) = newValue`).
 * This only works if the third argument is `true`. If it isn't, and the key does not exist, an error will be thrown.
 * Supports json arrays. If the field "key" does not exist, this function will try to convert it to integer and
 * access \p jsonObject at that index.
 * \param jsonObject The `json` object.
 * \param key The key to access.
 * \param mutator Whether to use mutator or accessor methods.
 * \return A reference to the value of the accessed key (will be a nullptr `json` if the key did not exist).
 */
nlohmann::json& valueAt( nlohmann::json& jsonObject, const std::string& key, const bool mutator = false );

//! Access a key path of a `json` object or array.
/*!
 * Access a key path of a `json` object or array. \p jsonObject is constant, so the returned value cannot be modified.
 * Supports json arrays. If the any key in \p keyPath does not exist, this function will try to convert it to integer
 * and access \p jsonObject at that index.
 * \param jsonObject The constant `json` object.
 * \param keyPath The key path to access.
 * \return The value at the accessed key path.
 */
nlohmann::json valueAt( nlohmann::json jsonObject, const KeyPath& keyPath );

//! Whether the key at \p keyPath is defined for \p jsonObject.
/*!
 * @copybrief isDefined
 * \param jsonObject The `json` object.
 * \param keyPath The key path specifying to key to be checked.
 * \return @copybrief isDefined
 */
bool isDefined( const nlohmann::json& jsonObject, const KeyPath& keyPath );


// SPECIAL KEYS ACCESS

//! Get the \p jsonObject at key SpecialKeys::rootObject.
/*!
 * copybrief getRootObject
 * \param jsonObject The `json` object.
 * \return JSON representation of the root object of \p jsonObject.
 * \throws UndefinedKeyError If the key SpecialKeys::rootObject is not defined for \p jsonObject.
 */
nlohmann::json getRootObject( const nlohmann::json& jsonObject );

//! Get the absolute key path from which \p jsonObject was retrieved.
/*!
 * copybrief getKeyPath
 * \param jsonObject The `json` object.
 * \return Absolute key path from which \p jsonObject was retrieved.
 * \throws UndefinedKeyError If the key SpecialKeys::keyPath is not defined for \p jsonObject.
 */
KeyPath getKeyPath( const nlohmann::json& jsonObject );

//! Get the key at which \p jsonObject was obtained.
/*!
 * copybrief getParentKey
 * \param jsonObject The `json` object.
 * \param errorMessage Error message to be printed if \p jsonObject has no key SpecialKeys::keyPath.
 * \return Key at which \p jsonObject was obtained.
 * \throws UndefinedKeyError If the key SpecialKeys::keyPath is not defined for \p jsonObject.
 */
std::string getParentKey( const nlohmann::json& jsonObject,
                          const std::string& errorMessage = "Could not determine parent key: context is missing." );

//! Get the response type to an event for a `json` object.
/*!
 * @copybrief getResponseToEvent
 * \param jsonObject The `json` object in which the type of response to the event is defined, or a `json` object with
 * a root object containing this information.
 * \param eventName The name of the event.
 * \param defaultResponse Default response if the response to the event is not defined. By default, continueSilently.
 * \return Response type to the event.
 */
ExceptionResponseType getResponseToEvent( const nlohmann::json& jsonObject, const std::string& eventName,
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
void checkUnusedKeys( const nlohmann::json& jsonObject, const ExceptionResponseType response );


// JSON ARRAY

//! Convert \p j to object if \p j is an array containing objects.
/*!
 * @copybrief convertToObjectIfContainsObjects This method does nothing if \p j is not an array, if \p j is empty, or if
 * the elements (and subelements) of \p j are not of type object.
 * \param j `json` object to be converted, or the original \p j.
 */
void convertToObjectIfContainsObjects( nlohmann::json& j );

//! Convert \p jsonObject to a json array.
/*!
 * Convert \p jsonObject to a json array. Does nothing if \p jsonObject is already an array.
 * If \p jsonObject is an object whose (non-special) keys are all convertible to unsigned int, the equivalent `json`
 * array will be returned. Otherise, the original \p jsonObject is returned.
 * \param jsonObject The `json` object to be converted.
 * \return \p jsonObject as a json array, or the original \p jsonObject if it is not convertible to array.
 */
nlohmann::json getAsArray( const nlohmann::json& jsonObject );

//! Whether \p j is convertible to a json array.
/*!
 * Whether \p j is convertible to a json array. True if \p j is an array or an object whose (non-special) keys are all
 * convertible to unsigned int.
 * \param j The `json` object to be checked.
 * \return Whether \p j is convertible to a json array.
 */
bool isConvertibleToArray( const nlohmann::json& j );


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
ValueType getValue( const nlohmann::json& jsonObject, const KeyPath& keyPath )
{
    nlohmann::json rootObject;
    nlohmann::json currentObject = jsonObject;
    KeyPath currentKeyPath = keyPath;
    KeyPath canonicalKeyPath = keyPath;

    if ( ! contains( currentKeyPath, SpecialKeys::keyPath ) )
    {
        canonicalKeyPath = currentKeyPath.canonical( getKeyPath( currentObject ) );

        if ( ! contains( currentKeyPath, SpecialKeys::rootObject ) )
        {
            if ( isDefined( currentObject, SpecialKeys::rootObject ) )
            {
                rootObject = getValue< nlohmann::json >( currentObject, SpecialKeys::rootObject );
            }

            if ( currentKeyPath.size( ) > 0 )
            {
                if ( currentKeyPath.front( ) == SpecialKeys::root || contains( currentKeyPath, SpecialKeys::up ) )
                {
                    // Path is absolute or contains up-key

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

    if ( ! containsAnyOf( currentKeyPath, SpecialKeys::all ) )
    {
        // If jsonObject is an array containing objects, convert to object
        convertToObjectIfContainsObjects( currentObject );

        if ( currentObject.is_object( ) )
        {
            // Define keys rootObject and keyPath of jsonObject to be returned
            currentObject[ SpecialKeys::rootObject ] = rootObject.is_null( ) ? jsonObject : rootObject;
            currentObject[ SpecialKeys::keyPath ] = canonicalKeyPath;
        }
    }

    try
    {
        // Try to convert to the requested type
        const ValueType convertedObject = currentObject;
        return convertedObject;
    }
    catch ( const nlohmann::detail::type_error& )
    {
        // Could not convert to the requested type
        throw IllegalValueError( canonicalKeyPath, currentObject, typeid( ValueType ) );
    }
    catch ( const UnknownEnumError& )
    {
        // Could not convert string to enum
        throw IllegalValueError( canonicalKeyPath, currentObject, typeid( ValueType ) );
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
ValueType getAs( const nlohmann::json& jsonObject )
{
    return getValue< ValueType >( jsonObject, { } );
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
ValueType getValue( const nlohmann::json& jsonObject, const KeyPath& keyPath, const ValueType& defaultValue )
{
    try
    {
        return getValue< ValueType >( jsonObject, keyPath );
    }
    catch ( const UndefinedKeyError& error )
    {
        error.rethrowIfNotTriggeredByMissingValueAt( keyPath );
        error.handleUseOfDefaultValue( nlohmann::json( defaultValue ),
                                       getResponseToEvent( jsonObject, Keys::Options::defaultValueUsedForMissingKey ) );
        return defaultValue;
    }
}

//! Get the value of \p jsonObject at the first defined key path in \p keyPaths.
/*!
 * Get the value of \p jsonObject at the first defined key path in \p keyPaths.
 * \param jsonObject The `json` object.
 * \param keyPaths Vector of key paths from which the value can be retrieved.
 * \return Value at the first defined key path in \p keyPaths.
 * \throws UndefinedKeyError If all the key paths in \p keyPaths are not defined for \p jsonObject. The printed message
 * will indicate the key path corresponding to the last of \p keyPaths.
 * \throws UndefinedKeyError If some of the subkeys needed to create the shared pointer of `ValueType` are missing
 * (only when \p jsonObject at \p keyPaths is of type object).
 * \throws IllegalValueError<ValueType> If the obtained value for the requested key path is not convertible to
 * `ValueType`.
 * \throws IllegalValueError<SubvalueType> If some of the subkeys needed to create the shared pointer of `ValueType`
 * are not convertible to `SubvalueType` (only when \p jsonObject at \p keyPaths is of type object).
 */
template< typename ValueType >
ValueType getValue( const nlohmann::json& jsonObject, const std::vector< KeyPath >& keyPaths )
{
    for ( unsigned int i = 0; i < keyPaths.size( ); ++i )
    {
        try
        {
            return getValue< nlohmann::json >( jsonObject, keyPaths.at( i ) );
        }
        catch ( const UndefinedKeyError& error )
        {
            error.rethrowIfNotTriggeredByMissingValueAt( keyPaths.at( i ) );
            if ( i == keyPaths.size( ) - 1 )
            {
                throw error;
            }
        }
    }
    throw;  // silence warning may reach end of non-void function
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
void updateFromJSON( T& value, const nlohmann::json& jsonObject, const KeyPath& keyPath = { } )
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
void updateFromJSONIfDefined( T& value, const nlohmann::json& jsonObject, const KeyPath& keyPath = { } )
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
void assignIfNotNaN( nlohmann::json& jsonObject, const std::string& key, const EquatableType& value )
{
    if ( ! isNaN( value ) )
    {
        jsonObject[ key ] = value;
    }
}

//! Assign \p object to `jsonObject[ key ]` if \p object is not nullptr.
/*!
 * @copybrief assignIfNotnullptr
 * \remark This function does nothing if \p object is `nullptr`.
 * \param jsonObject The `json` object being updated.
 * \param key The key of \p jsonObject being updated.
 * \param object Shared pointer to the object that is being to be used to update `jsonObject[ key ]`.
 */
template< typename T >
void assignIfNotnullptr( nlohmann::json& jsonObject, const std::string& key, const std::shared_ptr< T >& object )
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
void assignIfNotEmpty( nlohmann::json& jsonObject, const std::string& key, const ContainerType& container )
{
    if ( ! container.empty( ) )
    {
        jsonObject[ key ] = container;
    }
}


// Typedefs for single-body and body-to-body maps.

template < typename T >
using SingleBodyMap = std::unordered_map< std::string, std::vector< std::shared_ptr< T > > >;

template< typename T >
using BodyToBodyMap = std::unordered_map< std::string, SingleBodyMap< T > >;

template< typename T >
using NoPointerSingleBodyMap = std::unordered_map< std::string, std::vector< T > >;

template< typename T >
using NoPointerBodyToBodyMap = std::unordered_map< std::string, NoPointerSingleBodyMap< T > >;

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_VALUEACCESS_H
