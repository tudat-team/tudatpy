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

#include <boost/filesystem.hpp>

#include <Eigen/Core>

namespace tudat
{

namespace json_interface
{

//! Split \p string using \p delimiter.
/*!
 * Split \p string using \p delimiter.
 */
template< typename T >
void split( const std::string& string, char delimiter, T result )
{
    std::stringstream stream;
    stream.str( string );
    std::string item;
    while ( std::getline( stream, item, delimiter ) )
    {
        *( result++ ) = item;
    }
}

//! Get a vector containing the parts resulting from splitting \p string using \p delimiter.
/*!
 * Get a vector containing the parts resulting from splitting \p string using \p delimiter.
 * \param string The string to be splitted.
 * \param delimiter The delimiter to be used.
 * \return Vector containing the parts resulting from splitting \p string using \p delimiter.
 */
inline std::vector< std::string > split( const std::string& string, char delimiter )
{
    std::vector< std::string > parts;
    split( string, delimiter, std::back_inserter( parts ) );
    return parts;
}

//! Remove all entries of \p map except the last one.
/*!
 * @copybrief reduceToLast If \p map is empty, this function does nothing.
 */
template< typename K, typename V >
void reduceToLast( std::map< K, V > map )
{
    if ( map.size( ) > 0 )
    {
        map = { { ( --map.end( ) )->first, ( --map.end( ) )->second } };
    }
}

//! The key of \p map containing \p value.
/*!
 * \return @copybrief getKeyWithValue
 * \throws std::runtime_error If \p map does not contain \p value.
 */
template< typename K, typename V >
K getKeyWithValue( const std::map< K, V >& map, const V& value )
{
    for ( auto entry : map )
    {
        if ( entry.second == value )
        {
            return entry.first;
        }
    }
    throw std::runtime_error( "Could not find an entry in map with the requested value." );
}

//! Whether \p value is NaN.
/*!
 * \return @copybrief isNaN
 */
template< typename T >
bool isNaN( const T& value )
{
    return value != value;
}

//! Whether \p vector contains \p value.
/*!
 * \return @copybrief contains
 */
template< typename T >
bool contains( const std::vector< T >& vector, const T& value )
{
    return std::find( vector.begin( ), vector.end( ), value ) != vector.end( );
}

//! Whether \p vector contains any of \p values.
/*!
 * \return @copybrief containsAnyOf
 */
template< typename T >
bool containsAnyOf( const std::vector< T >& vector, std::vector< T > values )
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

//! Whether \p vector contains all of \p values.
/*!
 * \return @copybrief containsAllOf
 */
template< typename T >
bool containsAllOf( const std::vector< T >& vector, std::vector< T > values )
{
    for ( auto value : values )
    {
        if ( ! contains( vector, value ) )
        {
            return false;
        }
    }
    return true;
}

//! Get a vector containing the keys of an (un)ordered map.
/*!
 * @copybrief getMapKeys
 */
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

//! Get a vector containing the values of an (un)ordered map.
/*!
 * @copybrief getMapValues
 */
template< template < typename ... > class MapType, typename KeyType, typename ValueType >
std::vector< ValueType > getMapValues( const MapType< KeyType, ValueType >& map )
{
    std::vector< ValueType > values;
    for ( auto entry : map )
    {
        values.push_back( entry.second );
    }
    return values;
}

//! Get a vector containing the concatenated values of the vectors of an (un)ordered map.
/*!
 * @copybrief getFlattenedMapValues
 */
template< template < typename ... > class MapType, typename KeyType, typename ValueType >
std::vector< ValueType > getFlattenedMapValues( const MapType< KeyType, std::vector< ValueType > >& map )
{
    std::vector< ValueType > values;
    for ( auto entry : map )
    {
        for ( const ValueType value : entry.second )
        {
            values.push_back( value );
        }
    }
    return values;
}


// EIGEN <- STD

//! Create an Eigen column-vector from a `std::vector`.
/*!
 * @copybrief eigenVectorFromStdVector
 */
template< typename ScalarType >
Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > eigenVectorFromStdVector( const std::vector< ScalarType >& stdVector )
{
    return Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 >::Map( &stdVector[ 0 ], stdVector.size( ) );
}

//! Create an Eigen row-vector from a `std::vector`.
/*!
 * @copybrief eigenRowVectorFromStdVector
 */
template< typename ScalarType >
Eigen::Matrix< ScalarType, 1, Eigen::Dynamic > eigenRowVectorFromStdVector( const std::vector< ScalarType >& stdVector )
{
    return eigenVectorFromStdVector< ScalarType >( stdVector ).transpose( );
}

//! Create an Eigen matrix from a `std::vector` of `std::vector`s.
/*!
 * @copybrief eigenMatrixFromStdVectorOfVectors
 * \throws exception If all the elements of \p vectorOfVectors are not of the same size.
 */
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


// STD <- EIGEN

//! Create a `std::vector` from an Eigen matrix.
/*!
 * @copybrief stdVectorFromEigenMatrix
 */
template< typename ScalarType, int rows, int cols >
std::vector< ScalarType > stdVectorFromEigenMatrix(
        const Eigen::Matrix< ScalarType, rows, cols >& matrix )
{
    return std::vector< ScalarType >( matrix.data( ), matrix.data( ) + matrix.rows( ) * matrix.cols( ) );
}

//! Create a `std::vector` of `std::vector`s from an Eigen matrix.
/*!
 * @copybrief stdVectorOfVectorsFromEigenMatrix
 */
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


//! Encode a string for use in an URL.
/*!
 * @copybrief url_encode Reference: https://stackoverflow.com/questions/154536/encode-decode-urls-in-c
 */
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
