/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_UTILITIES_H
#define TUDAT_JSONINTERFACE_UTILITIES_H

#include <map>
#include <type_traits>

#include <boost/lexical_cast.hpp>

#include "json/src/json.hpp"
using json = nlohmann::json;

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
    //! Constructor.
    IllegalValueError( const KeyTree& keyTree ) : JSONError( "Illegal value for key", keyTree ) { }
};

/*
//! -DOC
template< typename T >
void split( const std::string& string, char delimiter, T result ) {
    std::stringstream stream;
    stream.str( string );
    std::string item;
    while ( std::getline( stream, item, delimiter ) )
    {
        *(result++) = item;
    }
}

//! -DOC
std::vector< std::string > split( const std::string& string, char delimiter )
{
    std::vector< std::string > parts;
    split( string, delimiter, std::back_inserter( parts ) );
    return parts;
}


//! -DOC
std::map< std::string, double > SIUnits =
{
    { "d", 86400 }
};

//! -DOC
template< typename T >
T convertToSIUnits( T number, const std::string& units )
{
    if ( std::is_arithmetic< T >::value )
    {
        return number * static_cast< T >( SIUnits.at( units ) );
    }
    else
    {
        throw;
    }
}

//! -DOC
template< typename T >
T parseNumber( const std::string& text )
{
    std::vector< std::string > parts = split( text, ' ' );
    T number = boost::lexical_cast< T >( parts.at( 0 ) );
    std::string units = parts.at( 1 );
    return convertToSIUnits< T >( number, units );
}
*/


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
template< typename T >
T getValue( json jsonObject, const KeyTree& keyTree )
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
        return jsonObject.get< T >( );
    }
    catch ( ... )
    {
        // Could not convert to the requested type
        try
        {
            // Convert to string (with units) and then parse the number and convert to SI
            // return parseNumber< T >( jsonObject.get< std::string >( ) );
            throw;
        }
        catch ( ... )
        {
            // Could not convert to the requested type nor parse the string as a number
            throw IllegalValueError( keyTree );
        }
    }
}


//! -DOC
template< typename T >
T getValue( const json& jsonObject, const KeyTree& keyTree, const T& optionalValue )
{
    try
    {
        return getValue< T >( jsonObject, keyTree );
    }
    catch ( const UndefinedKeyError& error )
    {
        return optionalValue;
    }
}


//! -DOC
template< typename T >
boost::shared_ptr< T > getValuePointer( const json& jsonObject, const KeyTree& keyTree )
{
    try
    {
        return boost::make_shared< T >( getValue< T >( jsonObject, keyTree ) );
    }
    catch ( const UndefinedKeyError& error )
    {
        return NULL;
    }
}


//! -DOC
template< typename T >
std::vector< T > set2vector( std::set< T > set )
{
    return std::vector< T >( set.begin(), set.end() );
}


} // namespace json_interfaces

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_UTILITIES_H
