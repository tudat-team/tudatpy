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


/// EIGEN <- STD

//! Create an Eigen column-vector from a `std::vector`.
template< typename ScalarType, int rows >
Eigen::Matrix< ScalarType, rows, 1 > eigenVectorFromStdVector( const std::vector< ScalarType >& stdVector )
{
    return Eigen::Matrix< ScalarType, rows, 1 >( stdVector.data( ) );
}

//! Create an Eigen row-vector from a `std::vector`.
template< typename ScalarType, int cols >
Eigen::Matrix< ScalarType, 1, cols > eigenRowVectorFromStdVector( const std::vector< ScalarType >& stdVector )
{
    return Eigen::Matrix< ScalarType, 1, cols >( stdVector.data( ) );
}

//! Create an Eigen matrix from a `std::vector` of `std::vector`s.
template< typename ScalarType, int rows, int cols >
Eigen::Matrix< ScalarType, rows, cols > eigenMatrixFromStdVectorOfVectors(
        const std::vector< std::vector< ScalarType > >& vectorOfVectors )
{
    Eigen::Matrix< ScalarType, rows, cols > matrix;
    for ( unsigned int r = 0; r < rows; ++r )
    {
        matrix.row( r ) = eigenRowVectorFromStdVector< ScalarType, cols >( vectorOfVectors.at( r ) );
    }
    return matrix;
}


/// STD <- EIGEN

/*
//! Create a `std::vector` from an Eigen column-vector.
template< typename ScalarType, int rows >
std::vector< ScalarType > stdVectorFromEigenColumnVector( const Eigen::Matrix< ScalarType, rows, 1 >& columnVector )
{
    return std::vector< ScalarType >( columnVector.data( ), columnVector.data( ) + columnVector.cols( ) );
}

//! Create a `std::vector` from an Eigen row-vector.
template< typename ScalarType, int cols >
std::vector< ScalarType > stdVectorFromEigenRowVector( const Eigen::Matrix< ScalarType, 1, cols >& rowVector )
{
    return std::vector< ScalarType >( rowVector.data( ), rowVector.data( ) + rowVector.rows( ) );
}
*/

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
    for ( unsigned int r = 0; r < rows; ++r )
    {
        vectorOfVectors.push_back(
                    stdVectorFromEigenMatrix( Eigen::Matrix< ScalarType, 1, cols >( matrix.row( r ) ) ) );
    }
    return vectorOfVectors;
}

} // namespace json_interface

} // namespace tudat


#endif // TUDAT_JSONINTERFACE_UTILITIES_H
