/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <boost/lexical_cast.hpp>

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/interface/json/support/path.h"
#include "tudat/interface/json/support/utilities.h"
#include "tudat/interface/json/support/valueAccess.h"

namespace tudat
{

namespace json_interface
{

// json <-> std::map, std::unordered_map

//! Template used by to_json methods for std::unodered_map and std::map.
//! Use of this function outside those methods is discouraged.
template< template< typename ... > class MapType, typename KeyType, typename ValueType >
nlohmann::json jsonFromMap( const MapType< KeyType, ValueType >& map )
{
    nlohmann::json jsonObject;
    for ( auto entry : map )
    {
        jsonObject[ boost::lexical_cast< std::string >( entry.first ) ] = entry.second;
    }
    return jsonObject;
}

template< template < typename ... > class VectorType, typename ValueType >
nlohmann::json jsonFromVector( const VectorType< ValueType >& vector )
{
    nlohmann::json jsonObject;

    for ( unsigned int i = 0; i < vector.size( ); i++ )
    {
        jsonObject.push_back( vector.at( i ) );
    }
    return jsonObject;
}




//! Template used by from_json methods for std::unodered_map and std::map.
//! Use of this function outside those methods is discouraged.
template< template< typename ... > class MapType, typename KeyType, typename ValueType >
MapType< KeyType, ValueType > mapFromJson( const nlohmann::json& jsonObject )
{
    MapType< KeyType, ValueType > map;
    for ( nlohmann::json::const_iterator it = jsonObject.begin( ); it != jsonObject.end( ); ++it )
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
void to_json( nlohmann::json& jsonObject, const unordered_map< KeyType, ValueType >& unorderedMap )
{
    jsonObject = tudat::json_interface::jsonFromMap< unordered_map, KeyType, ValueType >( unorderedMap );
}

//! Create a `std::map` with arbitrary key type from a `json` object.
template< typename KeyType, typename ValueType >
void from_json( const nlohmann::json& jsonObject, unordered_map< KeyType, ValueType >& unorderedMap )
{
    unorderedMap = tudat::json_interface::mapFromJson< unordered_map, KeyType, ValueType >( jsonObject );
}


// STD::MAP

//! Create a `json` object from a `std::map` with arbitrary key type.
template< typename KeyType, typename ValueType >
void to_json( nlohmann::json& jsonObject, const map< KeyType, ValueType >& orderedMap )
{
    jsonObject = tudat::json_interface::jsonFromMap< map, KeyType, ValueType >( orderedMap );
}

//! Create a `std::map` with arbitrary key type from a `json` object.
template< typename KeyType, typename ValueType >
void from_json( const nlohmann::json& jsonObject, map< KeyType, ValueType >& orderedMap )
{
    orderedMap = tudat::json_interface::mapFromJson< map, KeyType, ValueType >( jsonObject );
}


// STD::VECTOR

//! Create a `json` object from a `std::map` with arbitrary key type.
template< typename ValueType >
void to_json( nlohmann::json& jsonObject, const vector< ValueType >& vectorInput )
{
    jsonObject = tudat::json_interface::jsonFromVector< vector, ValueType >( vectorInput );
}

//! Create a `std::vector` from a `json` object.
template< typename ValueType >
void from_json( const nlohmann::json& jsonObject, vector< ValueType >& myVector )
{
    using namespace tudat::json_interface;
    const nlohmann::json jsonArray = getAsArray( jsonObject );
    if ( jsonArray.is_array( ) )
    {
        myVector.resize( jsonArray.size( ) );
        for ( unsigned int i = 0; i < jsonArray.size( ); ++i )
        {
            myVector[ i ] = getAs< ValueType >( jsonArray.at( i ) );
        }
    }
    else
    {
        std::cerr << "Expected an array." << std::endl;
        throw nlohmann::detail::type_error::create( 0, "", &jsonObject );
    }
}


// STD::PAIR

//! Create a `std::pair` from a `json` object.
template< typename V, typename W >
void from_json( const nlohmann::json& jsonObject, pair< V, W >& myPair )
{
    using namespace tudat::json_interface;
    const nlohmann::json jsonArray = getAsArray( jsonObject );
    myPair.first = jsonArray.at( 0 );
    myPair.second = jsonArray.at( 1 );
}


// STD::COMPLEX

//! Create a `json` object from a `std::complex`.
template< typename T >
void to_json( nlohmann::json& jsonObject, const complex< T >& complexNumber )
{
    jsonObject = boost::lexical_cast< string >( complexNumber );
}

//! Create a `std::complex` from a `json` object.
template< typename T >
void from_json( const nlohmann::json& jsonObject, complex< T >& complexNumber )
{
    complexNumber = boost::lexical_cast< complex< T > >( jsonObject.get< string >( ) );
}

}  // namespace std


namespace Eigen
{

//! Create a `json` object from an `Eigen::Matrix`.
template< typename ScalarType, int rows, int cols >
void to_json( nlohmann::json& jsonObject, const Matrix< ScalarType, rows, cols >& matrix )
{
    for ( int r = 0; r < matrix.rows( ); ++r )
    {
        nlohmann::json jsonArray;
        for ( int c = 0; c < matrix.cols( ); ++c )
        {
            jsonArray.push_back( matrix( r, c ) );
        }
        jsonObject.push_back( jsonArray );
    }
}

//! Create `Eigen::Matrix` from a `json` object.
template< typename ScalarType, int rows, int cols >
void from_json( const nlohmann::json& jsonObject, Matrix< ScalarType, rows, cols >& matrix )
{
    nlohmann::json jsonArray;
    if ( jsonObject.is_array( ) )
    {
        if ( jsonObject.empty( ) )
        {
            return;
        }
        jsonArray = jsonObject;
    }
    else if ( jsonObject.is_number( ) )
    {
        jsonArray.push_back( jsonObject );
    }
    else
    {
        throw nlohmann::detail::type_error::create( 0, "", &jsonObject );
    }

    nlohmann::json jsonArrayOfArrays;
    if ( jsonArray.front( ).is_array( ) )  // provided matrix
    {
        if( jsonArray.at( 0 ).size( ) == 1 )
        {
            for ( unsigned int i = 0; i < jsonArray.size( ); ++i )
            {
                jsonArrayOfArrays.push_back( { jsonArray.at( i ).at( 0 ) } );
            }
        }
        else
        {
            jsonArrayOfArrays = jsonArray;
        }
    }
    else  // provided vector
    {
        if ( rows == 1 )  // expected row vector
        {
            jsonArrayOfArrays.push_back( jsonArray );
        }
        else if ( cols == 1 )  // expected column vector
        {
            for ( unsigned int i = 0; i < jsonArray.size( ); ++i )
            {
                jsonArrayOfArrays.push_back( { jsonArray.at( i ) } );
            }
        }
        else  // expected matrix
        {
            std::cerr << "Expected a matrix, received a vector." << std::endl;
            throw nlohmann::detail::type_error::create( 0, "", &jsonObject );
        }
    }

    const unsigned int providedRows = jsonArrayOfArrays.size( );
    const unsigned int providedCols = jsonArrayOfArrays.front( ).size( );
    if ( ( rows >= 0 && int( providedRows ) != rows ) || ( cols >= 0 && int( providedCols ) != cols ) )
    {
        std::cerr << "Expected matrix of size " << rows << "x" << cols
                  << ", received matrix of size " << providedRows << "x" << providedCols << "." << std::endl;
        throw nlohmann::detail::type_error::create( 0, "", &jsonObject );
    }

    matrix.resize( providedRows, providedCols );
    for ( unsigned int r = 0; r < providedRows; ++r )
    {
        if ( jsonArrayOfArrays.at( r ).size( ) != providedCols )
        {
            std::cerr << "Unconsistent matrix size: some rows have different number of columns." << std::endl;
            throw nlohmann::detail::type_error::create( 0, "", &jsonObject );
        }
        for ( unsigned int c = 0; c < providedCols; ++c )
        {
            if( providedCols > 1 && providedRows > 1 )
            {
                matrix( r, c ) = jsonArrayOfArrays.at( r ).at( c );
            }
            else if( providedCols == 1 )
            {
                matrix( r, c ) = jsonArrayOfArrays.at( r );
            }
            else if( providedRows == 1 )
            {
                matrix( r, c ) = jsonArrayOfArrays.at( c );
            }
        }
    }
}


//! Create a `json` object from an `Eigen::Quaternion`.
inline void to_json( nlohmann::json& jsonObject, const Quaterniond& quaternion )
{
    // Get rotation matrix from quaternion and use that to initialise json object
    jsonObject = quaternion.toRotationMatrix( );
}

//! Create `Eigen::Quaternion` from a `json` object.
inline void from_json( const nlohmann::json& jsonObject, Quaterniond& quaternion )
{
    using namespace tudat::json_interface;
    using K = Keys::Body::RotationModel;

    if ( jsonObject.is_array( ) )
    {
        // Get rotation matrix from json object and use that to initialise quaternion
        quaternion = jsonObject.get< Eigen::Matrix3d >( );
    }
    else
    {
        quaternion = tudat::spice_interface::computeRotationQuaternionBetweenFrames(
                    getValue< std::string >( jsonObject, K::originalFrame ),
                    getValue< std::string >( jsonObject, K::targetFrame ),
                    getValue< double >( jsonObject, K::initialTime ) );
    }
}

}  // namespace Eigen

#endif // TUDAT_JSONINTERFACE_VALUECONVERSIONS_H
