/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include <sstream>
#include <string>
#include <iomanip>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <Eigen/Core>

namespace tudat
{

namespace json_interface
{

//! Split \p string using \p delimiter.
template< typename T >
void split( const std::string& string, const char delimiter, const bool trimSpaces, T result )
{
    std::stringstream stream;
    stream.str( string );
    std::string item;
    while ( std::getline( stream, item, delimiter ) )
    {
        *( result++ ) = trimSpaces ? boost::trim_copy( item ) : item;
    }
}

//! Get a vector containing the parts resulting from splitting \p string using \p delimiter.
/*!
 * Get a vector containing the parts resulting from splitting \p string using \p delimiter.
 * \param string The string to be splitted.
 * \param delimiter The delimiter to be used.
 * \param trim Whether to remove leading and trailing spaces from the returned strings.
 * \return Vector containing the parts resulting from splitting \p string using \p delimiter.
 */
inline std::vector< std::string > split( const std::string& string, const char delimiter, const bool trim = true )
{
    std::vector< std::string > parts;
    split( string, delimiter, trim, std::back_inserter( parts ) );
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
template< template< typename ... > class MapType, typename KeyType, typename ValueType >
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
template< template< typename ... > class MapType, typename KeyType, typename ValueType >
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
template< template< typename ... > class MapType, typename KeyType, typename ValueType >
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
        unsigned char uc = c;
        escaped << '%' << std::setw( 2 ) << int( uc );
        escaped << std::nouppercase;
    }

    return escaped.str( );
}

} // namespace json_interface

} // namespace tudat


#endif // TUDAT_JSONINTERFACE_UTILITIES_H
