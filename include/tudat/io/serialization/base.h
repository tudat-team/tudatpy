/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SERIALIZATION_BASE_H
#define TUDAT_SERIALIZATION_BASE_H

/**
 * @file serialization/base.h
 * @brief Core Boost serialization infrastructure: archive types, Eigen support,
 *        and serialize-to/from-string helpers.
 */

// Boost serialization headers
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

#include <Eigen/Core>

#include "tudat/basics/timeType.h"

#include <sstream>
#include <tuple>

// Provide std::tuple serialization in boost::serialization namespace
namespace boost
{
namespace serialization
{

// Helper to serialize a std::tuple element-by-element
namespace detail
{

template< std::size_t N >
struct tuple_serializer
{
    template< class Archive, typename... Args >
    static void serialize( Archive& ar, std::tuple< Args... >& t, const unsigned int version )
    {
        tuple_serializer< N - 1 >::serialize( ar, t, version );
        ar & std::get< N - 1 >( t );
    }
};

template<>
struct tuple_serializer< 0 >
{
    template< class Archive, typename... Args >
    static void serialize( Archive& ar, std::tuple< Args... >& t, const unsigned int version )
    {
        (void)ar; (void)t; (void)version;
    }
};

}  // namespace detail

template< class Archive, typename... Args >
void serialize( Archive& ar, std::tuple< Args... >& t, const unsigned int version )
{
    detail::tuple_serializer< sizeof...( Args ) >::serialize( ar, t, version );
}

}  // namespace serialization
}  // namespace boost

namespace tudat
{

namespace serialization
{

//! Serialize an Eigen matrix to a Boost archive
template< class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols >
void serialize( Archive& ar, Eigen::Matrix< Scalar, Rows, Cols, Options, MaxRows, MaxCols >& matrix, const unsigned int version )
{
    (void)version;

    Eigen::Index rows = matrix.rows( );
    Eigen::Index cols = matrix.cols( );

    ar & rows;
    ar & cols;

    if( Archive::is_loading::value )
    {
        matrix.resize( rows, cols );
    }

    for( Eigen::Index i = 0; i < rows; ++i )
    {
        for( Eigen::Index j = 0; j < cols; ++j )
        {
            ar & matrix( i, j );
        }
    }
}

//! Helper function to serialize an object to a binary string
template< typename T >
std::string serializeToBinaryString( const T& object )
{
    std::ostringstream oss( std::ios::binary );
    boost::archive::binary_oarchive oa( oss );
    oa << object;
    return oss.str( );
}

//! Helper function to deserialize an object from a binary string
template< typename T >
T deserializeFromBinaryString( const std::string& data )
{
    std::istringstream iss( data, std::ios::binary );
    boost::archive::binary_iarchive ia( iss );
    T object;
    ia >> object;
    return object;
}

//! Helper function to serialize a shared_ptr to a binary string (for polymorphic types)
template< typename T >
std::string serializeSharedPtrToBinaryString( const std::shared_ptr< T >& object )
{
    std::ostringstream oss( std::ios::binary );
    boost::archive::binary_oarchive oa( oss );
    oa << object;
    return oss.str( );
}

//! Helper function to deserialize a shared_ptr from a binary string (for polymorphic types)
template< typename T >
std::shared_ptr< T > deserializeSharedPtrFromBinaryString( const std::string& data )
{
    std::istringstream iss( data, std::ios::binary );
    boost::archive::binary_iarchive ia( iss );
    std::shared_ptr< T > object;
    ia >> object;
    return object;
}

//! Helper function to deserialize into an existing object from a binary string
template< typename T >
void deserializeFromBinaryString( const std::string& data, T& object )
{
    std::istringstream iss( data, std::ios::binary );
    boost::archive::binary_iarchive ia( iss );
    ia >> object;
}

}  // namespace serialization

}  // namespace tudat

// Provide Eigen serialization in boost::serialization namespace
namespace boost
{
namespace serialization
{

template< class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols >
void serialize( Archive& ar, Eigen::Matrix< Scalar, Rows, Cols, Options, MaxRows, MaxCols >& matrix, const unsigned int version )
{
    tudat::serialization::serialize( ar, matrix, version );
}

}  // namespace serialization
}  // namespace boost

#endif  // TUDAT_SERIALIZATION_BASE_H
