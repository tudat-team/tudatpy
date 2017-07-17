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

//! Create an `Eigen::VectorX` from an `std::vector`.
template< typename ScalarType, int rows >
Eigen::Matrix< ScalarType, rows, 1 > eigenFromStd( std::vector< ScalarType > stdVector )
{
    return Eigen::Matrix< ScalarType, rows, 1 >( stdVector.data( ) );
}

//! Create an `std::vector` from an `Eigen::VectorX`.
template< typename ScalarType >
std::vector< ScalarType > stdFromEigen( Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > eigenVector )
{
    return std::vector< ScalarType >( eigenVector.data( ),
                                      eigenVector.data( ) + eigenVector.rows( ) * eigenVector.cols( ) );
}

} // namespace json_interfaces

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_UTILITIES_H
