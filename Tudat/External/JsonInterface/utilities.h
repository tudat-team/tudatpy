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

#include <vector>
#include <set>

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

//! Returns whether a `vector` contains a `value`.
template< typename T >
bool contains( const std::vector< T >& vector, const T& value )
{
    return std::find( vector.begin( ), vector.end( ), value ) != vector.end( );
}

} // namespace json_interfaces

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_UTILITIES_H
