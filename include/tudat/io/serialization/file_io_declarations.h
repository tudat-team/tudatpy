/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved.
 */

#ifndef TUDAT_SERIALIZATION_FILE_IO_DECLARATIONS_H
#define TUDAT_SERIALIZATION_FILE_IO_DECLARATIONS_H

#include <memory>
#include <string>

// Declaration-only counterparts used by non-template model classes. Definitions live in the
// owning model library so public model headers do not pull in archive and stream implementations.
#if TUDAT_BUILD_WITH_SERIALIZATION

#define TUDAT_DECLARE_FILE_IO( ... )                              \
    void saveToBinary( const std::string& path ) const;           \
    static __VA_ARGS__ loadFromBinary( const std::string& path ); \
    void saveToJson( const std::string& path ) const;             \
    static __VA_ARGS__ loadFromJson( const std::string& path );

#define TUDAT_DECLARE_BINARY_IO( ... )                  \
    void saveToBinary( const std::string& path ) const; \
    static __VA_ARGS__ loadFromBinary( const std::string& path );

#define TUDAT_DECLARE_BINARY_IO_POLYMORPHIC( ... )      \
    void saveToBinary( const std::string& path ) const; \
    static std::shared_ptr< __VA_ARGS__ > loadFromBinary( const std::string& path );

#define TUDAT_DECLARE_FILE_IO_POLYMORPHIC( ... )                                     \
    void saveToBinary( const std::string& path ) const;                              \
    static std::shared_ptr< __VA_ARGS__ > loadFromBinary( const std::string& path ); \
    void saveToJson( const std::string& path ) const;                                \
    static std::shared_ptr< __VA_ARGS__ > loadFromJson( const std::string& path );

#else

#define TUDAT_DECLARE_FILE_IO( ... )
#define TUDAT_DECLARE_BINARY_IO( ... )
#define TUDAT_DECLARE_BINARY_IO_POLYMORPHIC( ... )
#define TUDAT_DECLARE_FILE_IO_POLYMORPHIC( ... )

#endif

#endif  // TUDAT_SERIALIZATION_FILE_IO_DECLARATIONS_H
