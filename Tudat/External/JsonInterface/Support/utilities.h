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
#include <set>
#include <vector>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <string>

#include <eigen/Eigen>

namespace tudat
{

namespace json_interface
{

/*
//! Convert a `std::set` to a `std::vecotr`.
template< typename T >
std::vector< T > set2vector( std::set< T > set )
{
    return std::vector< T >( set.begin(), set.end() );
}
*/

//! Returns whether a value is NaN.
template< typename T >
bool isnan( const T& value )
{
    return value != value;
}

//! Returns whether a `vector` contains a `value`.
template< typename T >
bool contains( const std::vector< T >& vector, const T& value )
{
    return std::find( vector.begin( ), vector.end( ), value ) != vector.end( );
}

//! Returns whether a `vector` contains any of `values`.
template< typename T >
bool contains( const std::vector< T >& vector, std::vector< T > values )
{
    for ( auto value : values )
    {
        if ( contains( vector, value ) )
        {
            return true;
        }
    }
    return false;
}


//! Push back an element to a vector.
template< typename T >
void pushBackElements( std::vector< T >& vector, const T& element )
{
    vector.push_back( element );
}

//! Push back the elements of a vector to another vector.
template< typename T >
void pushBackElements( std::vector< T >& vector, const std::vector< T >& elements )
{
    for ( const T element : elements )
    {
        pushBackElements( vector, element );
    }
}



//! Get a vector containing the keys of a map / unordered_map.
template< template < typename ... > class MapType, typename KeyType, typename ValueType >
std::vector< KeyType > getMapKeys( const MapType< KeyType, ValueType >& map )
{
    std::vector< KeyType > keys;
    for ( auto entry : map )
    {
        keys.push_back( entry.first );
    }
    return keys;
}

//! Get a vector containing the values of a map / unordered_map.
template< template < typename ... > class MapType, typename KeyType, typename ValueType >
std::vector< ValueType > getMapValues( const MapType< KeyType, ValueType >& map, bool flatten = false )
{
    std::vector< ValueType > values;
    for ( auto entry : map )
    {
        if ( flatten )
        {
            pushBackElements( values, entry.second );
        }
        else
        {
            values.push_back( entry.second );
        }
    }
    return values;
}


/// EIGEN <- STD

//! Create an Eigen column-vector from a `std::vector`.
template< typename ScalarType >
Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > eigenVectorFromStdVector( const std::vector< ScalarType >& stdVector )
{
    return Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 >::Map( &stdVector[ 0 ], stdVector.size( ) );
}

//! Create an Eigen row-vector from a `std::vector`.
template< typename ScalarType >
Eigen::Matrix< ScalarType, 1, Eigen::Dynamic > eigenRowVectorFromStdVector( const std::vector< ScalarType >& stdVector )
{
    return eigenVectorFromStdVector< ScalarType >( stdVector ).transpose( );
}

//! Create an Eigen matrix from a `std::vector` of `std::vector`s.
template< typename ScalarType >
Eigen::Matrix< ScalarType, Eigen::Dynamic, Eigen::Dynamic > eigenMatrixFromStdVectorOfVectors(
        const std::vector< std::vector< ScalarType > >& vectorOfVectors )
{
    const unsigned int rows = vectorOfVectors.size( );
    const unsigned int cols = vectorOfVectors.at( 0 ).size( );
    Eigen::Matrix< ScalarType, Eigen::Dynamic, Eigen::Dynamic > matrix( rows, cols );
    for ( unsigned int r = 0; r < rows; ++r )
    {
        if ( vectorOfVectors.at( r ).size( ) != cols )
        {
            std::cerr << "Unconsistent matrix size." << std::endl;
            throw;
        }
        matrix.row( r ) = eigenRowVectorFromStdVector< ScalarType >( vectorOfVectors.at( r ) );
    }
    return matrix;
}


/// STD <- EIGEN

//! Create a `std::vector` from an Eigen matrix.
template< typename ScalarType, int rows, int cols >
std::vector< ScalarType > stdVectorFromEigenMatrix(
        const Eigen::Matrix< ScalarType, rows, cols >& matrix )
{
    return std::vector< ScalarType >( matrix.data( ), matrix.data( ) + matrix.rows( ) * matrix.cols( ) );
}

//! Create a `std::vector` of `std::vector`s from an Eigen matrix.
template< typename ScalarType, int rows, int cols >
std::vector< std::vector< ScalarType > > stdVectorOfVectorsFromEigenMatrix(
        const Eigen::Matrix< ScalarType, rows, cols >& matrix )
{
    std::vector< std::vector< ScalarType > > vectorOfVectors;
    for ( unsigned int r = 0; r < matrix.rows( ); ++r )
    {
        vectorOfVectors.push_back(
                    stdVectorFromEigenMatrix( Eigen::Matrix< ScalarType, 1, cols >( matrix.row( r ) ) ) );
    }
    return vectorOfVectors;
}


//! Encode a string for use in URL [ https://stackoverflow.com/questions/154536/encode-decode-urls-in-c ]
inline std::string url_encode( const std::string &value )
{
    std::ostringstream escaped;
    escaped.fill( '0' );
    escaped << std::hex;

    for ( std::string::const_iterator i = value.begin(), n = value.end(); i != n; ++i )
    {
        std::string::value_type c = (*i);

        // Keep alphanumeric and other accepted characters intact
        if ( isalnum(c) || c == '-' || c == '_' || c == '.' || c == '~' )
        {
            escaped << c;
            continue;
        }

        // Any other characters are percent-encoded
        escaped << std::uppercase;
        escaped << '%' << std::setw( 2 ) << int( ( unsigned char ) c );
        escaped << std::nouppercase;
    }

    return escaped.str( );
}

} // namespace json_interface

} // namespace tudat


#endif // TUDAT_JSONINTERFACE_UTILITIES_H
