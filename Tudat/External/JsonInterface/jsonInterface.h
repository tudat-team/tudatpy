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

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "json/src/json.hpp"
using json = nlohmann::json;

#include <Tudat/InputOutput/basicInputOutput.h>

#include "keys.h"
#include "path.h"
#include "units.h"
#include "utilities.h"

namespace tudat
{

namespace json_interface
{

//! -DOC
class JSONError : public std::runtime_error
{
private:
    const KeyTree keyTree;

public:
    //! Constructor.
    JSONError( const std::string& errorMessage, const KeyTree& keyTree )
        : runtime_error( errorMessage.c_str( ) ), keyTree( keyTree ) { }

    //! Error message.
    virtual const char* what( ) const throw( )
    {
        std::ostringstream stream;
        stream << runtime_error::what( ) << ": " << keyTree;
        std::cerr << stream.str( ).c_str( ) << std::endl;  // FIXME
        return stream.str( ).c_str( );
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

    //! Error message.
    virtual const char* what( ) const throw( )
    {
        std::ostringstream stream;
        stream << JSONError::what( ) << " = " << value;
        // std::cerr << stream.str( ).c_str( ) << std::endl;  // FIXME
        return stream.str( ).c_str( );
    }
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
        // Could not convert to the requested type
        throw IllegalValueError( keyTree, jsonObject );
    }
}


//! -DOC
template< typename NumberType >
NumberType getNumeric( json jsonObject, const KeyTree& keyTree )
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
            return parseMagnitudeWithUnits< NumberType >( error.value.get< std::string >( ) );
        }
        catch ( ... )
        {
            // Could not convert to the requested type nor parse the string as a number
            throw IllegalValueError( keyTree, error.value );
        }
    }
}


//! -DOC
template< typename NumberType >
NumberType getEpoch( json jsonObject, const KeyTree& keyTree )
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
            // Convert to string and then parse as a formatted date
            return convertToSecondsSinceJ2000< NumberType >( error.value.get< std::string >( ) );
        }
        catch ( ... )
        {
            // Could not convert to the requested type nor parse the string as a date
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
NumberType getNumeric( const json& jsonObject, const KeyTree& keyTree,
                      const NumberType& optionalValue, bool allowNaN = false )
{
    try
    {
        return getNumeric< NumberType >( jsonObject, keyTree );
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
NumberType getEpoch( const json& jsonObject, const KeyTree& keyTree,
                      const NumberType& optionalValue, bool allowNaN = false )
{
    try
    {
        return getEpoch< NumberType >( jsonObject, keyTree );
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
boost::shared_ptr< NumberType > getNumericPointer( const json& jsonObject, const KeyTree& keyTree )
{
    try
    {
        return boost::make_shared< NumberType >( getNumeric< NumberType >( jsonObject, keyTree ) );
    }
    catch ( const UndefinedKeyError& error )
    {
        return NULL;
    }
}


//! -DOC
template< typename NumberType >
boost::shared_ptr< NumberType > getEpochPointer( const json& jsonObject, const KeyTree& keyTree )
{
    try
    {
        return boost::make_shared< NumberType >( getEpoch< NumberType >( jsonObject, keyTree ) );
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
std::string stringFromEnum( const EnumType enumValue, const std::map< std::string, EnumType >& possibleValues )
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

} // namespace json_interface

} // namespace tudat


namespace Eigen
{

//! Create a `json` object from an `Eigen::Matrix`.
//! Called automatically by `nlohmann::json` when using `jsonObject = json( matrix )`.
template< typename ScalarType, int rows, int cols >
void to_json( json& jsonObject, const Matrix< ScalarType, rows, cols >& matrix )
{
    // Convert to std::vector of std::vector's and use that to initialise json object
    jsonObject = json( tudat::json_interface::stdVectorOfVectorsFromEigenMatrix( matrix ) );
}

//! Create an Eigen matrix from a `json` object.
//! Called automatically by `nlohmann::json` when using `matrix = jsonObject.get< Eigen::Matrix >( )`.
template< typename ScalarType, int rows, int cols >
void from_json( const json& jsonObject, Matrix< ScalarType, rows, cols >& matrix )
{
    int providedRows, providedCols;
    bool transposed = false;
    // Get as std::vector of std::vector's and then convert to Eigen matrix
    try
    {
        const std::vector< std::vector< ScalarType > > vectorOfVectors =
                jsonObject.get< std::vector< std::vector< ScalarType > > >( );
        matrix = tudat::json_interface::eigenMatrixFromStdVectorOfVectors< ScalarType >( vectorOfVectors );
        providedRows = vectorOfVectors.size( );
        providedCols = vectorOfVectors.at( 0 ).size( );
    }
    catch ( ... )
    {
        const std::vector< ScalarType > vector = jsonObject.get< std::vector< ScalarType > >( );
        // Get as std::vector and then convert to Eigen column-vector
        if ( cols == 1 )
        {
            matrix.col( 0 ) = tudat::json_interface::eigenVectorFromStdVector< ScalarType >( vector );
            providedRows = vector.size( );
            providedCols = 1;
            transposed = true;
        }
        // Get as std::vector and then convert to Eigen row-vector
        else if ( rows == 1 )
        {
            matrix.row( 0 ) = tudat::json_interface::eigenRowVectorFromStdVector< ScalarType >( vector );
            providedRows = 1;
            providedCols = vector.size( );
            transposed = false;
        }
        else
        {
            std::cerr << "Could not convert JSON array (of arrays) to Eigen vector/matrix." << std::endl;
            throw;
        }
    }
    if ( ( rows >= 0 && providedRows != rows ) || ( cols >= 0 && providedCols != cols ) )
    {
        std::cerr << "Expected matrix of size " << rows << "x" << cols;
        if ( rows == 1 || cols == 1 )
        {
            std::cerr << " or " << cols << "x" << rows;
        }
        std::cerr << ", received matrix of size ";
        if ( transposed )
        {
            std::cerr << providedCols << "x" << providedRows;
        }
        else
        {
            std::cerr << providedRows << "x" << providedCols;
        }
        std::cerr << "." << std::endl;
        throw;
    }
}

}  // namespace Eigen


namespace std
{

//! Create a `json` object from a `std::map` with arbitrary key type.
//! Called automatically by `nlohmann::json` when using `jsonObject = json( map )`.
template< typename KeyType, typename ValueType >
void to_json( json& jsonObject, const std::map< KeyType, ValueType >& map )
{
    for ( auto entry : map )
    {
        jsonObject[ boost::lexical_cast< std::string >( entry.first ) ] = entry.second;
    }
}

//! Create a std::map with arbitrary key type from a `json` object.
//! Called automatically by `nlohmann::json` when using `map = jsonObject.get< std::map< KeyType, ValueType > >( )`.
template< typename KeyType, typename ValueType >
void from_json( const json& jsonObject, std::map< KeyType, ValueType >& map )
{
    json j = jsonObject;
    for ( json::iterator it = j.begin( ); it != j.end( ); ++it )
    {
        map[ boost::lexical_cast< KeyType >( it.key( ) ) ] = it.value( ).get< ValueType >( );
    }
}

}  // namespace std

#endif // TUDAT_JSONINTERFACE_VALUES_H
