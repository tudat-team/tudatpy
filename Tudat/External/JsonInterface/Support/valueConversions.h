/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_VALUECONVERSIONS_H
#define TUDAT_JSONINTERFACE_VALUECONVERSIONS_H

#include "json/src/json.hpp"
using json = nlohmann::json;

#include "path.h"
#include "utilities.h"

namespace Eigen
{

//! Create a `json` object from an `Eigen::Matrix`.
//! Called automatically by `nlohmann::json` when using `jsonObject = matrix`.
template< typename ScalarType, int rows, int cols >
void to_json( json& jsonObject, const Matrix< ScalarType, rows, cols >& matrix )
{
    // Convert to std::vector of std::vector's and use that to initialise json object
    jsonObject = tudat::json_interface::stdVectorOfVectorsFromEigenMatrix( matrix );
}

//! Create `Eigen::Matrix` from a `json` object.
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


//! Create a `json` object from an `Eigen::Quaternion`.
//! Called automatically by `nlohmann::json` when using `jsonObject = quaternion`.
template< typename ScalarType >
void to_json( json& jsonObject, const Quaternion< ScalarType >& quaternion )
{
    // Get rotation matrix from quaternion and use that to initialise json object
    jsonObject = quaternion.toRotationMatrix( );
}

//! Create `Eigen::Quaternion` from a `json` object.
//! Called automatically by `nlohmann::json` when using `quaternion = jsonObject.get< Eigen::Quaternion >( )`.
template< typename ScalarType >
void from_json( const json& jsonObject, Quaternion< ScalarType >& quaternion )
{
    // Get rotation matrix from json object and use that to initialise json object
    quaternion = jsonObject.get< Eigen::Matrix< ScalarType, 3, 3 > >( );
}

}  // namespace Eigen


namespace std
{

/// Support for maps with arbitrary key type

//! Create a `json` object from a `std::map` with arbitrary key type.
//! Called automatically by `nlohmann::json` when using `jsonObject = myMap`.
template< typename KeyType, typename ValueType >
void to_json( json& jsonObject, const map< KeyType, ValueType >& myMap )
{
    for ( auto entry : myMap )
    {
        jsonObject[ boost::lexical_cast< string >( entry.first ) ] = entry.second;
    }
}

//! Create a std::map with arbitrary key type from a `json` object.
//! Called automatically by `nlohmann::json` when using `myMap = jsonObject.get< std::map< KeyType, ValueType > >( )`.
template< typename KeyType, typename ValueType >
void from_json( const json& jsonObject, map< KeyType, ValueType >& myMap )
{
    json j = jsonObject;
    for ( json::iterator it = j.begin( ); it != j.end( ); ++it )
    {
        myMap[ boost::lexical_cast< KeyType >( it.key( ) ) ] = it.value( ).get< ValueType >( );
    }
}


/// Support for complex numbers

//! Create a `json` object from a `std::complex`.
//! Called automatically by `nlohmann::json` when using `jsonObject = complexNumber`.
template< typename T >
void to_json( json& jsonObject, const complex< T >& complexNumber )
{
    jsonObject = boost::lexical_cast< string >( complexNumber );
}

//! Create a `std::complex` from a `json` object.
//! Called automatically by `nlohmann::json` when using `complexNumber = jsonObject.get< std::complex< double > >( )`.
template< typename T >
void from_json( const json& jsonObject, complex< T >& complexNumber )
{
    complexNumber = boost::lexical_cast< complex< T > >( jsonObject.get< string >( ) );
}

}  // namespace std

#endif // TUDAT_JSONINTERFACE_VALUECONVERSIONS_H
