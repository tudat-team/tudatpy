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
#include "valueAccess.h"

namespace tudat
{

namespace json_interface
{

//! Template used by to_json methods for std::unodered_map and std::map.
//! Use of this function outside those methods is discouraged.
template< template < typename ... > class MapType, typename KeyType, typename ValueType >
json jsonFromMap( const MapType< KeyType, ValueType >& map )
{
    json jsonObject;
    for ( auto entry : map )
    {
        jsonObject[ boost::lexical_cast< std::string >( entry.first ) ] = entry.second;
    }
    return jsonObject;
}

//! Template used by from_json methods for std::unodered_map and std::map.
//! Use of this function outside those methods is discouraged.
template< template < typename ... > class MapType, typename KeyType, typename ValueType >
MapType< KeyType, ValueType > mapFromJson( const json& jsonObject )
{
    MapType< KeyType, ValueType > map;
    for ( json::const_iterator it = jsonObject.begin( ); it != jsonObject.end( ); ++it )
    {
        const std::string key = it.key( );
        if ( ! contains( SpecialKeys::all, key ) )
        {
            map[ boost::lexical_cast< KeyType >( key ) ] = getValue< ValueType >( jsonObject, key );
        }
    }
    return map;
}

}  // namespace json_interface

}  // namespace tudat


namespace std
{

// STD::UNORDERED_MAP

//! Create a `json` object from a `std::unordered_map` with arbitrary key type.
template< typename KeyType, typename ValueType >
void to_json( json& jsonObject, const unordered_map< KeyType, ValueType >& unorderedMap )
{
    jsonObject = tudat::json_interface::jsonFromMap< unordered_map, KeyType, ValueType >( unorderedMap );
}

//! Create a `std::map` with arbitrary key type from a `json` object.
template< typename KeyType, typename ValueType >
void from_json( const json& jsonObject, unordered_map< KeyType, ValueType >& unorderedMap )
{
    unorderedMap = tudat::json_interface::mapFromJson< unordered_map, KeyType, ValueType >( jsonObject );
}


// STD::MAP

//! Create a `json` object from a `std::map` with arbitrary key type.
template< typename KeyType, typename ValueType >
void to_json( json& jsonObject, const map< KeyType, ValueType >& orderedMap )
{
    jsonObject = tudat::json_interface::jsonFromMap< map, KeyType, ValueType >( orderedMap );
}

//! Create a `std::map` with arbitrary key type from a `json` object.
template< typename KeyType, typename ValueType >
void from_json( const json& jsonObject, map< KeyType, ValueType >& orderedMap )
{
    orderedMap = tudat::json_interface::mapFromJson< map, KeyType, ValueType >( jsonObject );
}


// STD::VECTOR

//! Create a `std::vector` from a `json` object.
template< typename ValueType >
void from_json( const json& jsonObject, vector< ValueType >& myVector )
{
    using namespace tudat::json_interface;

    if ( jsonObject.is_array( ) )
    {
        for ( unsigned int i = 0; i < jsonObject.size( ); ++i )
        {
            myVector.push_back( getValue< ValueType >( jsonObject, i ) );
        }
    }
    else
    {
        bool isObjectWithIntConvertibleKeys = jsonObject.is_object( );
        if ( isObjectWithIntConvertibleKeys )
        {
            for ( json::const_iterator it = jsonObject.begin( ); it != jsonObject.end( ); ++it )
            {
                const std::string key = it.key( );
                if ( ! contains( SpecialKeys::all, key ) )
                {
                    try
                    {
                        stoi( key );
                    }
                    catch ( ... )
                    {
                        isObjectWithIntConvertibleKeys = false;
                        break;
                    }
                }
            }
        }

        if ( isObjectWithIntConvertibleKeys )
        {
            // Convert to map, and use that to create the vector
            const map< string, ValueType > auxiliaryMap = getAs< map< string, ValueType > >( jsonObject );
            for ( auto entry : auxiliaryMap )
            {
                if ( ! contains( SpecialKeys::all, entry.first ) )
                {
                    myVector.push_back( entry.second );
                }
            }
        }
        else
        {
            throw nlohmann::detail::type_error::create( 0, "" );
        }
    }
}


// STD::PAIR

//! Create a `std::pair` from a `json` object.
template< typename V, typename W >
void from_json( const json& jsonObject, pair< V, W >& myPair )
{
    using namespace tudat::json_interface;
    myPair.first = getValue< V >( jsonObject, "0" );
    myPair.second = getValue< W >( jsonObject, "1" );
}


// STD::COMPLEX

//! Create a `json` object from a `std::complex`.
template< typename T >
void to_json( json& jsonObject, const complex< T >& complexNumber )
{
    jsonObject = boost::lexical_cast< string >( complexNumber );
}

//! Create a `std::complex` from a `json` object.
template< typename T >
void from_json( const json& jsonObject, complex< T >& complexNumber )
{
    complexNumber = boost::lexical_cast< complex< T > >( jsonObject.get< string >( ) );
}

}  // namespace std


namespace Eigen
{

//! Create a `json` object from an `Eigen::Matrix`.
template< typename ScalarType, int rows, int cols >
void to_json( json& jsonObject, const Matrix< ScalarType, rows, cols >& matrix )
{
    // Convert to std::vector of std::vector's and use that to initialise json object
    jsonObject = tudat::json_interface::stdVectorOfVectorsFromEigenMatrix( matrix );
}

//! Create `Eigen::Matrix` from a `json` object.
template< typename ScalarType, int rows, int cols >
void from_json( const json& jsonObject, Matrix< ScalarType, rows, cols >& matrix )
{
    bool transposed = false;
    // Get as std::vector of std::vector's and then convert to Eigen matrix
    std::vector< std::vector< ScalarType > > vectorOfVectors;
    try
    {
        vectorOfVectors = jsonObject.get< std::vector< std::vector< ScalarType > > >( );
    }
    catch ( ... )
    {
        const std::vector< ScalarType > vector = jsonObject.get< std::vector< ScalarType > >( );
        if ( cols == 1 )
        {
            for ( const ScalarType element : vector )
            {
                vectorOfVectors.push_back( { element } );
            }
            transposed = true;
        }
        else if ( rows == 1 )
        {
            vectorOfVectors = { vector };
        }
        else
        {
            std::cerr << "Could not convert JSON array (of arrays) to Eigen vector/matrix." << std::endl;
            throw nlohmann::detail::type_error::create( 0, "" );
        }
    }
    matrix = tudat::json_interface::eigenMatrixFromStdVectorOfVectors< ScalarType >( vectorOfVectors );
    const int providedRows = vectorOfVectors.size( );
    const int providedCols = vectorOfVectors.at( 0 ).size( );

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
        throw nlohmann::detail::type_error::create( 0, "" );
    }
}


//! Create a `json` object from an `Eigen::Quaternion`.
template< typename ScalarType >
void to_json( json& jsonObject, const Quaternion< ScalarType >& quaternion )
{
    // Get rotation matrix from quaternion and use that to initialise json object
    jsonObject = quaternion.toRotationMatrix( );
}

//! Create `Eigen::Quaternion` from a `json` object.
template< typename ScalarType >
void from_json( const json& jsonObject, Quaternion< ScalarType >& quaternion )
{
    // Get rotation matrix from json object and use that to initialise json object
    quaternion = jsonObject.get< Eigen::Matrix< ScalarType, 3, 3 > >( );
}

}  // namespace Eigen

#endif // TUDAT_JSONINTERFACE_VALUECONVERSIONS_H
