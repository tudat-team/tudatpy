/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_H
#define TUDAT_JSONINTERFACE_H

#include <set>

#include "json/src/json.hpp"
using json = nlohmann::json;

#include "units.h"
#include "utilities.h"

namespace tudat
{

namespace json_interface
{

// FIXME: what about arrays?
//! Class for specifying a key tree (key.subkey.subsubkey ...) used to access data from `json` objects.
class KeyTree : public std::vector< std::string >
{
public:
    //! Inherit constructors.
    using vector< std::string >::vector;

    //! Constructor with a single string key.
    /*!
     * Constructor with a single string key.
     * \param key The key to be accessed.
     */
    KeyTree( const std::string& key ) : KeyTree( std::vector< std::string >( { key } ) ) { }

    //! Constructor with a single char key.
    /*!
     * Constructor with a single char key.
     * \param key The key to be accessed.
     */
    //! Constructor with a single char key.
    KeyTree( const char* key ) : KeyTree( std::string( key ) ) { }
};

//! String representation for `KeyTree`, as key.subkey.subsubkey ...
inline std::ostream& operator<< ( std::ostream & stringRepresentation, KeyTree const & keyTree ) {
    for ( unsigned int i = 0; i < keyTree.size(); ++i )
    {
        stringRepresentation << keyTree.at( i );
        if ( i < keyTree.size() - 1 )
        {
            stringRepresentation << ".";
        }
    }
    return stringRepresentation;
}


//! -DOC
class JSONError : public std::runtime_error
{
private:
    const KeyTree keyTree;

public:
    //! Constructor.
    JSONError( const std::string& errorMessage, const KeyTree& keyTree )
        : runtime_error( errorMessage.c_str() ), keyTree( keyTree ) { }

    //! Error message.
    virtual const char* what() const throw()
    {
        std::ostringstream stream;
        stream << runtime_error::what() << ": " << keyTree;
        std::cerr << stream.str().c_str() << std::endl;  // FIXME
        return stream.str().c_str();
    }
};

//! -DOC
class UndefinedKeyError : public JSONError
{
public:
    //! Constructor.
    UndefinedKeyError( const KeyTree& keyTree ) : JSONError( "Undefined key", keyTree ) { }
};

//! -DOC
class IllegalValueError : public JSONError
{
public:
    //! Associated (illegal) value.
    json value;

    //! Constructor.
    IllegalValueError( const KeyTree& keyTree, const json& value )
        : JSONError( "Illegal value for key", keyTree ), value( value ) { }
};


//! Get the value of a parameter defined by a `KeyTree` from a `json` object.
/*!
 * Get the value of a parameter defined by a `KeyTree` from a `json` object.
 * An error will be thrown if the requested key does not exist, or the provided value is not of the expected type `T`.
 * \param jsonObject JSON object from which to get the value.
 * \param keyTree Vector of keys defining the value to be accessed (key.subkey.subsubkey ...).
 * \return Value of the requested key.
 * \throw UndefinedKeyError If the requested key is not defined.
 * \throw IllegalValueError If the provided value for the requested key is not of type `T`.
 */
template< typename ValueType >
ValueType getValue( json jsonObject, const KeyTree& keyTree )
{
    try
    {
        // Recursively update jsonObject for every key in keyTree
        for ( unsigned int i = 0; i < keyTree.size(); ++i )
        {
            jsonObject = jsonObject.at( keyTree.at( i ) );
        }
    }
    catch ( ... )
    {
        // At least one of the keys was missing
        throw UndefinedKeyError( keyTree );
    }
    try
    {
        // Convert to the requested type
        return jsonObject.get< ValueType >( );
    }
    catch ( ... )
    {
        // Could not convert to the requested type nor parse the string as a number
        throw IllegalValueError( keyTree, jsonObject );
    }
}


//! -DOC
template< typename NumberType >
NumberType getNumber( json jsonObject, const KeyTree& keyTree )
{
    try
    {
        return getValue< NumberType >( jsonObject, keyTree );
    }
    catch ( const IllegalValueError& error )
    {
        // Could not convert to the requested type
        try
        {
            // Convert to string (with units) and then parse the number and convert to SI
            return parseString< NumberType >( error.value.get< std::string >( ) );
        }
        catch ( ... )
        {
            // Could not convert to the requested type nor parse the string as a number
            throw IllegalValueError( keyTree, error.value );
        }
    }
}


//! -DOC
template< typename ValueType >
ValueType getValue( const json& jsonObject, const KeyTree& keyTree, const ValueType& optionalValue )
{
    try
    {
        return getValue< ValueType >( jsonObject, keyTree );
    }
    catch ( const UndefinedKeyError& error )
    {
        return optionalValue;
    }
}


//! -DOC
template< typename NumberType >
NumberType getNumber( const json& jsonObject, const KeyTree& keyTree, const NumberType& optionalValue )
{
    try
    {
        return getNumber< NumberType >( jsonObject, keyTree );
    }
    catch ( const UndefinedKeyError& error )
    {
        return optionalValue;
    }
}


//! -DOC
template< typename ValueType >
boost::shared_ptr< ValueType > getValuePointer( const json& jsonObject, const KeyTree& keyTree )
{
    try
    {
        return boost::make_shared< ValueType >( getValue< ValueType >( jsonObject, keyTree ) );
    }
    catch ( const UndefinedKeyError& error )
    {
        return NULL;
    }
}


//! -DOC
template< typename NumberType >
boost::shared_ptr< NumberType > getNumberPointer( const json& jsonObject, const KeyTree& keyTree )
{
    try
    {
        return boost::make_shared< NumberType >( getNumber< NumberType >( jsonObject, keyTree ) );
    }
    catch ( const UndefinedKeyError& error )
    {
        return NULL;
    }
}


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
        std::cerr << "Unrecognized value (" + stringValue + ") for " << typeid( EnumType ).name( ) << std::endl;
        throw;
    }
}

//! -DOC
template< typename EnumType >
std::string stringFromEnum( const EnumType enumValue,
                                          const std::map< std::string, EnumType >& possibleValues )
{
    for ( auto ent : possibleValues )
    {
        if ( ent.second == enumValue )
        {
            return ent.first;
        }
    }
    std::cerr << "Unrecognized name for " << typeid( EnumType ).name( ) << " = " << enumValue << std::endl;
    throw;
}


} // namespace json_interfaces

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_VALUES_H
