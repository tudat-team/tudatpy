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

/*
//! Get a vector containing the keys of a map.
template< typename KeyType, typename ValueType >
std::vector< KeyType > getKeys( const std::map< KeyType, ValueType >& map )
{
    std::vector< KeyType > keys;
    for ( auto entry : map )
    {
        keys.push_back( entry.first );
    }
    return keys;
}

//! Get a vector containing the values of a map.
template< typename KeyType, typename ValueType >
std::vector< ValueType > getValues( const std::map< KeyType, ValueType >& map )
{
    std::vector< ValueType > values;
    for ( auto entry : map )
    {
        values.push_back( entry.second );
    }
    return values;
}
*/


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

} // namespace json_interface

} // namespace tudat


#endif // TUDAT_JSONINTERFACE_UTILITIES_H
