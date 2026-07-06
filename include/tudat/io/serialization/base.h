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
 * @brief Core Cereal serialization infrastructure: archive types, Eigen support,
 *        std::tuple support, and serialize-to/from-string helpers.
 */

// Cereal core
#include <cereal/cereal.hpp>
#include <cereal/access.hpp>

// Cereal archive types
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>

// Cereal type support
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>

#include <Eigen/Core>

#include "tudat/basics/timeType.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

// Provide Eigen matrix serialization for cereal
namespace cereal
{

//! Serialize an Eigen matrix (save)
template< class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols >
void save( Archive& ar, const Eigen::Matrix< Scalar, Rows, Cols, Options, MaxRows, MaxCols >& matrix )
{
    Eigen::Index rows = matrix.rows( );
    Eigen::Index cols = matrix.cols( );

    ar( rows );
    ar( cols );

    for( Eigen::Index i = 0; i < rows; ++i )
    {
        for( Eigen::Index j = 0; j < cols; ++j )
        {
            ar( matrix( i, j ) );
        }
    }
}

//! Deserialize an Eigen matrix (load)
template< class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols >
void load( Archive& ar, Eigen::Matrix< Scalar, Rows, Cols, Options, MaxRows, MaxCols >& matrix )
{
    Eigen::Index rows, cols;

    ar( rows );
    ar( cols );

    matrix.resize( rows, cols );

    for( Eigen::Index i = 0; i < rows; ++i )
    {
        for( Eigen::Index j = 0; j < cols; ++j )
        {
            ar( matrix( i, j ) );
        }
    }
}

}  // namespace cereal

namespace tudat
{

namespace serialization
{

//! Helper function to serialize an object to a binary string
template< typename T >
std::string serializeToBinaryString( const T& object )
{
    std::ostringstream oss( std::ios::binary );
    {
        cereal::BinaryOutputArchive oa( oss );
        oa( object );
    }  // archive must go out of scope to flush
    return oss.str( );
}

//! Helper function to deserialize an object from a binary string
template< typename T >
T deserializeFromBinaryString( const std::string& data )
{
    std::istringstream iss( data, std::ios::binary );
    cereal::BinaryInputArchive ia( iss );
    T object;
    ia( object );
    return object;
}

//! Helper function to serialize a shared_ptr to a binary string (for polymorphic types)
template< typename T >
std::string serializeSharedPtrToBinaryString( const std::shared_ptr< T >& object )
{
    std::ostringstream oss( std::ios::binary );
    {
        cereal::BinaryOutputArchive oa( oss );
        oa( object );
    }
    return oss.str( );
}

//! Helper function to deserialize a shared_ptr from a binary string (for polymorphic types)
template< typename T >
std::shared_ptr< T > deserializeSharedPtrFromBinaryString( const std::string& data )
{
    std::istringstream iss( data, std::ios::binary );
    cereal::BinaryInputArchive ia( iss );
    std::shared_ptr< T > object;
    ia( object );
    return object;
}

//! Helper function to deserialize into an existing object from a binary string
template< typename T >
void deserializeFromBinaryString( const std::string& data, T& object )
{
    std::istringstream iss( data, std::ios::binary );
    cereal::BinaryInputArchive ia( iss );
    ia( object );
}

//! Helper function to serialize an object to a binary file
//! Appends ".tudat" extension automatically.
template< typename T >
void saveToBinaryFile( const T& object, const std::string& path )
{
    std::string filePath = path + ".tudat";
    std::ofstream outputStream( filePath, std::ios::binary );
    if( !outputStream )
    {
        throw std::runtime_error( "Unable to open file for binary save: " + filePath );
    }

    cereal::BinaryOutputArchive archive( outputStream );
    archive( object );
}

//! Helper function to deserialize an object from a binary file
//! Appends ".tudat" extension automatically.
template< typename T >
T loadFromBinaryFile( const std::string& path )
{
    std::string filePath = path + ".tudat";
    std::ifstream inputStream( filePath, std::ios::binary );
    if( !inputStream )
    {
        throw std::runtime_error( "Unable to open file for binary load: " + filePath );
    }

    cereal::BinaryInputArchive archive( inputStream );
    T object;
    archive( object );
    return object;
}

//! Helper function to deserialize from a binary file into an existing object
//! Appends ".tudat" extension automatically.
template< typename T >
void loadFromBinaryFile( const std::string& path, T& object )
{
    std::string filePath = path + ".tudat";
    std::ifstream inputStream( filePath, std::ios::binary );
    if( !inputStream )
    {
        throw std::runtime_error( "Unable to open file for binary load: " + filePath );
    }

    cereal::BinaryInputArchive archive( inputStream );
    archive( object );
}

// =====================================================================
//  JSON helpers
// =====================================================================

//! Helper function to serialize an object to a JSON file
//! Appends ".json" extension automatically.
template< typename T >
void saveToJsonFile( const T& object, const std::string& path )
{
    std::string filePath = path + ".json";
    std::ofstream outputStream( filePath );
    if( !outputStream )
    {
        throw std::runtime_error( "Unable to open file for JSON save: " + filePath );
    }

    cereal::JSONOutputArchive archive( outputStream );
    archive( cereal::make_nvp( "root", object ) );
}

//! Helper function to deserialize an object from a JSON file
//! Appends ".json" extension automatically.
template< typename T >
T loadFromJsonFile( const std::string& path )
{
    std::string filePath = path + ".json";
    std::ifstream inputStream( filePath );
    if( !inputStream )
    {
        throw std::runtime_error( "Unable to open file for JSON load: " + filePath );
    }

    cereal::JSONInputArchive archive( inputStream );
    T object;
    archive( cereal::make_nvp( "root", object ) );
    return object;
}

//! Helper function to deserialize from a JSON file into an existing object
//! Appends ".json" extension automatically.
template< typename T >
void loadFromJsonFile( const std::string& path, T& object )
{
    std::string filePath = path + ".json";
    std::ifstream inputStream( filePath );
    if( !inputStream )
    {
        throw std::runtime_error( "Unable to open file for JSON load: " + filePath );
    }

    cereal::JSONInputArchive archive( inputStream );
    archive( cereal::make_nvp( "root", object ) );
}

//! Helper function to serialize an object to a JSON string
template< typename T >
std::string serializeToJsonString( const T& object )
{
    std::ostringstream oss;
    {
        cereal::JSONOutputArchive oa( oss );
        oa( cereal::make_nvp( "root", object ) );
    }
    return oss.str( );
}

//! Helper function to deserialize an object from a JSON string
template< typename T >
T deserializeFromJsonString( const std::string& data )
{
    std::istringstream iss( data );
    cereal::JSONInputArchive ia( iss );
    T object;
    ia( cereal::make_nvp( "root", object ) );
    return object;
}

}  // namespace serialization

}  // namespace tudat

#endif  // TUDAT_SERIALIZATION_BASE_H
