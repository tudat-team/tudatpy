/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved.
 */

#ifndef TUDAT_SERIALIZATION_EIGEN_H
#define TUDAT_SERIALIZATION_EIGEN_H

/**
 * @file serialization/eigen.h
 * @brief Archive-independent Cereal serialization support for Eigen matrices.
 */

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

#include <Eigen/Core>

#include "tudat/io/serialization/core.h"

namespace cereal
{

//! Guard against malformed archives requesting unreasonable allocations.
constexpr std::uintmax_t kMaximumSerializedEigenCoefficients = 100000000;

template< class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols >
void save( Archive& ar, const Eigen::Matrix< Scalar, Rows, Cols, Options, MaxRows, MaxCols >& matrix )
{
    Eigen::Index rows = matrix.rows( );
    Eigen::Index cols = matrix.cols( );

    ar( make_nvp( "rows", rows ) );
    ar( make_nvp( "cols", cols ) );

    for( Eigen::Index i = 0; i < rows; ++i )
    {
        for( Eigen::Index j = 0; j < cols; ++j )
        {
            ar( make_nvp( "element", matrix( i, j ) ) );
        }
    }
}

template< class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols >
void load( Archive& ar, Eigen::Matrix< Scalar, Rows, Cols, Options, MaxRows, MaxCols >& matrix )
{
    Eigen::Index rows, cols;

    ar( make_nvp( "rows", rows ) );
    ar( make_nvp( "cols", cols ) );

    if( rows < 0 || cols < 0 )
    {
        throw std::runtime_error( "Cannot deserialize Eigen matrix: archived dimensions must be non-negative." );
    }
    if( Rows != Eigen::Dynamic && rows != Rows )
    {
        throw std::runtime_error( "Cannot deserialize Eigen matrix: archived row count " + std::to_string( rows ) +
                                  " does not match fixed row count " + std::to_string( Rows ) + "." );
    }
    if( Cols != Eigen::Dynamic && cols != Cols )
    {
        throw std::runtime_error( "Cannot deserialize Eigen matrix: archived column count " + std::to_string( cols ) +
                                  " does not match fixed column count " + std::to_string( Cols ) + "." );
    }
    if( MaxRows != Eigen::Dynamic && rows > MaxRows )
    {
        throw std::runtime_error( "Cannot deserialize Eigen matrix: archived row count exceeds the matrix maximum." );
    }
    if( MaxCols != Eigen::Dynamic && cols > MaxCols )
    {
        throw std::runtime_error( "Cannot deserialize Eigen matrix: archived column count exceeds the matrix maximum." );
    }

    const std::uintmax_t unsignedRows = static_cast< std::uintmax_t >( rows );
    const std::uintmax_t unsignedCols = static_cast< std::uintmax_t >( cols );
    if( unsignedRows > kMaximumSerializedEigenCoefficients || unsignedCols > kMaximumSerializedEigenCoefficients )
    {
        throw std::runtime_error( "Cannot deserialize Eigen matrix: an archived dimension exceeds the supported maximum of " +
                                  std::to_string( kMaximumSerializedEigenCoefficients ) + "." );
    }
    if( unsignedRows != 0 && unsignedCols > kMaximumSerializedEigenCoefficients / unsignedRows )
    {
        throw std::runtime_error( "Cannot deserialize Eigen matrix: archived dimensions contain more than " +
                                  std::to_string( kMaximumSerializedEigenCoefficients ) + " coefficients." );
    }

    matrix.resize( rows, cols );

    for( Eigen::Index i = 0; i < rows; ++i )
    {
        for( Eigen::Index j = 0; j < cols; ++j )
        {
            ar( make_nvp( "element", matrix( i, j ) ) );
        }
    }
}

}  // namespace cereal

#endif  // TUDAT_SERIALIZATION_EIGEN_H
