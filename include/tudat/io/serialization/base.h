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

#include <tudat/config.hpp>

#include <Eigen/Core>

#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

// Provide Eigen matrix serialization for cereal
namespace cereal
{

//! Serialize an Eigen matrix (save)
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

//! Deserialize an Eigen matrix (load)
template< class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols >
void load( Archive& ar, Eigen::Matrix< Scalar, Rows, Cols, Options, MaxRows, MaxCols >& matrix )
{
    Eigen::Index rows, cols;

    ar( make_nvp( "rows", rows ) );
    ar( make_nvp( "cols", cols ) );

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
    std::unique_ptr< T > objectPtr( cereal::access::construct< T >( ) );
    ia( *objectPtr );
    return std::move( *objectPtr );
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

// =====================================================================
//  Provenance helpers (not in anonymous namespace — template two-phase
//  lookup requires them at namespace scope).
// =====================================================================

//! Magic number for binary serialization files ("TED1").
constexpr std::uint32_t kBinaryMagic = 0x54454431;

//! Provenance string for the current build: "version@commit@date".
inline std::string currentProvenance( )
{
    static const std::string s =
            std::string( TUDAT_VERSION ) + "@" + std::string( TUDAT_BUILD_COMMIT ) + "@" + std::string( TUDAT_BUILD_TIME );
    return s;
}

//! Write the binary version header to an output stream.
inline void writeBinaryHeader( std::ostream& stream )
{
    const std::string p = currentProvenance( );
    const std::uint32_t len = static_cast< std::uint32_t >( p.size( ) );
    stream.write( reinterpret_cast< const char* >( &kBinaryMagic ), sizeof( kBinaryMagic ) );
    stream.write( reinterpret_cast< const char* >( &len ), sizeof( len ) );
    stream.write( p.data( ), len );
}

//! Read and return the provenance from a binary header, or return an
//! empty string if no header is present (legacy file).  Rewinds on
//! mismatch.
inline std::string readBinaryProvenance( std::istream& stream )
{
    const auto startPos = stream.tellg( );
    std::uint32_t magic = 0;
    stream.read( reinterpret_cast< char* >( &magic ), sizeof( magic ) );
    if( !stream || magic != kBinaryMagic )
    {
        stream.clear( );
        stream.seekg( startPos );
        return {};
    }
    std::uint32_t len = 0;
    stream.read( reinterpret_cast< char* >( &len ), sizeof( len ) );
    if( !stream )
    {
        throw std::runtime_error( "Corrupt binary header: could not read provenance length" );
    }
    std::string provenance( len, '\0' );
    stream.read( provenance.data( ), len );
    return provenance;
}

//! Extract the provenance from a raw JSON string by searching for
//! "tudat_provenance".  Returns empty string if not found.
inline std::string extractJsonProvenance( const std::string& raw )
{
    const std::string key = "\"tudat_provenance\": \"";
    auto pos = raw.find( key );
    if( pos == std::string::npos ) return {};
    pos += key.size( );
    auto end = raw.find( '\"', pos );
    if( end == std::string::npos ) return {};
    return raw.substr( pos, end - pos );
}

//! Parsed components of a provenance string.
struct ProvenanceParts {
    std::string version;
    std::string commit;
    std::string date;
};

//! Split a "version@commit@date" string into its three parts.
inline ProvenanceParts splitProvenance( const std::string& p )
{
    ProvenanceParts parts;
    auto at1 = p.find( '@' );
    if( at1 == std::string::npos ) return parts;
    parts.version = p.substr( 0, at1 );
    auto at2 = p.find( '@', at1 + 1 );
    if( at2 == std::string::npos ) return parts;
    parts.commit = p.substr( at1 + 1, at2 - at1 - 1 );
    parts.date = p.substr( at2 + 1 );
    return parts;
}

//! Check that the file has a provenance header, warn on mismatch, and
//! return the parsed parts.
inline ProvenanceParts checkProvenance( const std::string& fileProvenance )
{
    if( fileProvenance.empty( ) )
    {
        throw std::runtime_error(
                "Could not deserialize this object.\n"
                "This file is a legacy serialized file without a Tudat "
                "provenance header.\n"
                "Files created with the current version cannot be read "
                "by this build." );
    }

    static bool warnedOnce = false;
    if( !warnedOnce )
    {
        const std::string& current = currentProvenance( );
        if( fileProvenance != current )
        {
            std::cerr << "[Tudat] Warning: loading serialized file "
                         "created with:\n"
                      << "    " << fileProvenance << "\n"
                      << "  Current build is:\n"
                      << "    " << current << "\n"
                      << "  There is no guarantee of backward "
                         "compatibility beyond 6 months.\n"
                      << std::endl;
        }
        warnedOnce = true;
    }

    return splitProvenance( fileProvenance );
}

//! Helper function to serialize an object to a binary file.
//! Prepends a header with provenance "version@commit@date".
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
    writeBinaryHeader( outputStream );
    cereal::BinaryOutputArchive archive( outputStream );
    archive( object );
}

//! Helper function to deserialize an object from a binary file.
//! Warns on provenance mismatch and wraps errors with file provenance.
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

    const std::string fileProvenance = readBinaryProvenance( inputStream );
    const auto parts = checkProvenance( fileProvenance );

    try
    {
        cereal::BinaryInputArchive archive( inputStream );
        std::unique_ptr< T > objectPtr( cereal::access::construct< T >( ) );
        archive( *objectPtr );
        return std::move( *objectPtr );
    }
    catch( const std::exception& e )
    {
        std::string msg = "Could not deserialize this object.\n";
        if( !parts.commit.empty( ) )
        {
            msg += "This file was created with commit " + parts.commit + ".\n";
        }
        msg += "Error: " + std::string( e.what( ) );
        throw std::runtime_error( msg );
    }
}

//! Helper function to deserialize from a binary file into an existing
//! object.  Warns on provenance mismatch.
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

    const std::string fileProvenance = readBinaryProvenance( inputStream );
    const auto parts = checkProvenance( fileProvenance );

    try
    {
        cereal::BinaryInputArchive archive( inputStream );
        archive( object );
    }
    catch( const std::exception& e )
    {
        std::string msg = "Could not deserialize this object.\n";
        if( !parts.commit.empty( ) )
        {
            msg += "This file was created with commit " + parts.commit + ".\n";
        }
        msg += "Error: " + std::string( e.what( ) );
        throw std::runtime_error( msg );
    }
}

//! Helper function to serialize an object to a JSON file.
//! Embeds provenance "version@commit@date" as archive metadata.
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
    archive( cereal::make_nvp( "tudat_provenance", currentProvenance( ) ) );
    archive( cereal::make_nvp( "root", object ) );
}

//! Helper function to deserialize an object from a JSON file.
//! Warns on provenance mismatch and wraps errors with file provenance.
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

    // Read entire file to extract provenance (cereal JSON is DOM-based,
    // so we can't peek at the stream before handing it to the archive).
    std::string raw( ( std::istreambuf_iterator< char >( inputStream ) ), std::istreambuf_iterator< char >( ) );
    inputStream.close( );

    const std::string fileProvenance = extractJsonProvenance( raw );
    const auto parts = checkProvenance( fileProvenance );

    try
    {
        std::istringstream iss( raw );
        cereal::JSONInputArchive archive( iss );
        std::unique_ptr< T > objectPtr( cereal::access::construct< T >( ) );
        archive( cereal::make_nvp( "root", *objectPtr ) );
        return std::move( *objectPtr );
    }
    catch( const std::exception& e )
    {
        std::string msg = "Could not deserialize this object.\n";
        if( !parts.commit.empty( ) )
        {
            msg += "This file was created with commit " + parts.commit + ".\n";
        }
        msg += "Error: " + std::string( e.what( ) );
        throw std::runtime_error( msg );
    }
}

//! Helper function to deserialize from a JSON file into an existing object.
//! Warns on provenance mismatch.
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

    std::string raw( ( std::istreambuf_iterator< char >( inputStream ) ), std::istreambuf_iterator< char >( ) );
    inputStream.close( );

    const std::string fileProvenance = extractJsonProvenance( raw );
    const auto parts = checkProvenance( fileProvenance );

    try
    {
        std::istringstream iss( raw );
        cereal::JSONInputArchive archive( iss );
        archive( cereal::make_nvp( "root", object ) );
    }
    catch( const std::exception& e )
    {
        std::string msg = "Could not deserialize this object.\n";
        if( !parts.commit.empty( ) )
        {
            msg += "This file was created with commit " + parts.commit + ".\n";
        }
        msg += "Error: " + std::string( e.what( ) );
        throw std::runtime_error( msg );
    }
}

//! Helper function to serialize a shared_ptr to a JSON file (for polymorphic types)
//! Appends ".json" extension automatically.
template< typename T >
void saveSharedPtrToJsonFile( const std::shared_ptr< T >& object, const std::string& path )
{
    saveToJsonFile< std::shared_ptr< T > >( object, path );
}

//! Helper function to deserialize a shared_ptr from a JSON file (for polymorphic types)
//! Appends ".json" extension automatically.
template< typename T >
std::shared_ptr< T > loadSharedPtrFromJsonFile( const std::string& path )
{
    return loadFromJsonFile< std::shared_ptr< T > >( path );
}

//! Helper function to serialize a shared_ptr to a binary file (for polymorphic types)
//! Appends ".tudat" extension automatically.
template< typename T >
void saveSharedPtrToBinaryFile( const std::shared_ptr< T >& object, const std::string& path )
{
    saveToBinaryFile< std::shared_ptr< T > >( object, path );
}

//! Helper function to deserialize a shared_ptr from a binary file (for polymorphic types)
//! Appends ".tudat" extension automatically.
template< typename T >
std::shared_ptr< T > loadSharedPtrFromBinaryFile( const std::string& path )
{
    return loadFromBinaryFile< std::shared_ptr< T > >( path );
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
    std::unique_ptr< T > objectPtr( cereal::access::construct< T >( ) );
    ia( cereal::make_nvp( "root", *objectPtr ) );
    return std::move( *objectPtr );
}

}  // namespace serialization

}  // namespace tudat

// =====================================================================
//  Convenience macros for adding file-IO member functions to classes
//  that already have cereal save/load (private) serialization.
//
//  Usage: place the macro inside the class body (any visibility section).
//  Use __VA_ARGS__ so template types with commas (e.g. Foo<A,B>) work.
//
//  Four patterns:
//    TUDAT_DEFINE_FILE_IO(ClassName)              — binary + JSON, value
//    TUDAT_DEFINE_BINARY_IO(ClassName)            — binary only, value
//    TUDAT_DEFINE_JSON_IO(ClassName)              — JSON only, value
//    TUDAT_DEFINE_JSON_IO_POLYMORPHIC(BaseName)   — JSON only, shared_ptr<Base>
//    TUDAT_DEFINE_BINARY_IO_POLYMORPHIC(BaseName) — binary only, shared_ptr<Base>
//
//  For templated classes pass the full instantiation, e.g.:
//    TUDAT_DEFINE_BINARY_IO(SingleArcSimulationResults<StateScalarType, TimeType>)
// =====================================================================

//! Add both saveToBinary/loadFromBinary and saveToJson/loadFromJson (value types).
#define TUDAT_DEFINE_FILE_IO( ... )                                               \
    void saveToBinary( const std::string& path ) const                            \
    {                                                                             \
        ::tudat::serialization::saveToBinaryFile( *this, path );                  \
    }                                                                             \
    static __VA_ARGS__ loadFromBinary( const std::string& path )                  \
    {                                                                             \
        return ::tudat::serialization::loadFromBinaryFile< __VA_ARGS__ >( path ); \
    }                                                                             \
    void saveToJson( const std::string& path ) const                              \
    {                                                                             \
        ::tudat::serialization::saveToJsonFile( *this, path );                    \
    }                                                                             \
    static __VA_ARGS__ loadFromJson( const std::string& path )                    \
    {                                                                             \
        return ::tudat::serialization::loadFromJsonFile< __VA_ARGS__ >( path );   \
    }

//! Add saveToBinary/loadFromBinary only (value types).
#define TUDAT_DEFINE_BINARY_IO( ... )                                             \
    void saveToBinary( const std::string& path ) const                            \
    {                                                                             \
        ::tudat::serialization::saveToBinaryFile( *this, path );                  \
    }                                                                             \
    static __VA_ARGS__ loadFromBinary( const std::string& path )                  \
    {                                                                             \
        return ::tudat::serialization::loadFromBinaryFile< __VA_ARGS__ >( path ); \
    }

//! Add saveToJson/loadFromJson only (value types).
#define TUDAT_DEFINE_JSON_IO( ... )                                             \
    void saveToJson( const std::string& path ) const                            \
    {                                                                           \
        ::tudat::serialization::saveToJsonFile( *this, path );                  \
    }                                                                           \
    static __VA_ARGS__ loadFromJson( const std::string& path )                  \
    {                                                                           \
        return ::tudat::serialization::loadFromJsonFile< __VA_ARGS__ >( path ); \
    }

//! Add saveToJson/loadFromJson for polymorphic base (load returns shared_ptr<Base>).
//! Throws if deserialization produces a null pointer.
#define TUDAT_DEFINE_JSON_IO_POLYMORPHIC( ... )                                                                              \
    void saveToJson( const std::string& path ) const                                                                         \
    {                                                                                                                        \
        auto shared = std::shared_ptr< __VA_ARGS__ >( const_cast< __VA_ARGS__* >( this ), []( __VA_ARGS__* ) {} );           \
        ::tudat::serialization::saveSharedPtrToJsonFile< __VA_ARGS__ >( shared, path );                                      \
    }                                                                                                                        \
    static std::shared_ptr< __VA_ARGS__ > loadFromJson( const std::string& path )                                            \
    {                                                                                                                        \
        auto result = ::tudat::serialization::loadSharedPtrFromJsonFile< __VA_ARGS__ >( path );                              \
        if( !result )                                                                                                        \
        {                                                                                                                    \
            throw std::runtime_error( "Failed to deserialize " #__VA_ARGS__ " from JSON file '" + path +                     \
                                      "'.\n"                                                                                 \
                                      "The file was likely written through a base-class pointer of a different hierarchy.\n" \
                                      "Use the base class's loadFromJson() instead." );                                      \
        }                                                                                                                    \
        return result;                                                                                                       \
    }

//! Add saveToBinary/loadFromBinary for polymorphic base (load returns shared_ptr<Base>).
//! Throws if deserialization produces a null pointer.
#define TUDAT_DEFINE_BINARY_IO_POLYMORPHIC( ... )                                                                            \
    void saveToBinary( const std::string& path ) const                                                                       \
    {                                                                                                                        \
        auto shared = std::shared_ptr< __VA_ARGS__ >( const_cast< __VA_ARGS__* >( this ), []( __VA_ARGS__* ) {} );           \
        ::tudat::serialization::saveSharedPtrToBinaryFile< __VA_ARGS__ >( shared, path );                                    \
    }                                                                                                                        \
    static std::shared_ptr< __VA_ARGS__ > loadFromBinary( const std::string& path )                                          \
    {                                                                                                                        \
        auto result = ::tudat::serialization::loadSharedPtrFromBinaryFile< __VA_ARGS__ >( path );                            \
        if( !result )                                                                                                        \
        {                                                                                                                    \
            throw std::runtime_error( "Failed to deserialize " #__VA_ARGS__ " from binary file '" + path +                   \
                                      "'.\n"                                                                                 \
                                      "The file was likely written through a base-class pointer of a different hierarchy.\n" \
                                      "Use the base class's loadFromBinary() instead." );                                    \
        }                                                                                                                    \
        return result;                                                                                                       \
    }

//! Add saveToBinary/loadFromBinary AND saveToJson/loadFromJson for polymorphic base
//! (load returns shared_ptr<Base>).  Throws on null.
#define TUDAT_DEFINE_FILE_IO_POLYMORPHIC( ... )                                                                              \
    void saveToBinary( const std::string& path ) const                                                                       \
    {                                                                                                                        \
        auto shared = std::shared_ptr< __VA_ARGS__ >( const_cast< __VA_ARGS__* >( this ), []( __VA_ARGS__* ) {} );           \
        ::tudat::serialization::saveSharedPtrToBinaryFile< __VA_ARGS__ >( shared, path );                                    \
    }                                                                                                                        \
    static std::shared_ptr< __VA_ARGS__ > loadFromBinary( const std::string& path )                                          \
    {                                                                                                                        \
        auto result = ::tudat::serialization::loadSharedPtrFromBinaryFile< __VA_ARGS__ >( path );                            \
        if( !result )                                                                                                        \
        {                                                                                                                    \
            throw std::runtime_error( "Failed to deserialize " #__VA_ARGS__ " from binary file '" + path +                   \
                                      "'.\n"                                                                                 \
                                      "The file was likely written through a base-class pointer of a different hierarchy.\n" \
                                      "Use the base class's loadFromBinary() instead." );                                    \
        }                                                                                                                    \
        return result;                                                                                                       \
    }                                                                                                                        \
    void saveToJson( const std::string& path ) const                                                                         \
    {                                                                                                                        \
        auto shared = std::shared_ptr< __VA_ARGS__ >( const_cast< __VA_ARGS__* >( this ), []( __VA_ARGS__* ) {} );           \
        ::tudat::serialization::saveSharedPtrToJsonFile< __VA_ARGS__ >( shared, path );                                      \
    }                                                                                                                        \
    static std::shared_ptr< __VA_ARGS__ > loadFromJson( const std::string& path )                                            \
    {                                                                                                                        \
        auto result = ::tudat::serialization::loadSharedPtrFromJsonFile< __VA_ARGS__ >( path );                              \
        if( !result )                                                                                                        \
        {                                                                                                                    \
            throw std::runtime_error( "Failed to deserialize " #__VA_ARGS__ " from JSON file '" + path +                     \
                                      "'.\n"                                                                                 \
                                      "The file was likely written through a base-class pointer of a different hierarchy.\n" \
                                      "Use the base class's loadFromJson() instead." );                                      \
        }                                                                                                                    \
        return result;                                                                                                       \
    }

#endif  // TUDAT_SERIALIZATION_BASE_H
