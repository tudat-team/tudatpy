// serialization/pybind11_helpers.h
#pragma once
#include <pybind11/pybind11.h>
#include "tudat/io/serialization/base.h"

namespace py = pybind11;

namespace tudat::serialization
{

// For plain value types (non-polymorphic)
template< typename T >
auto make_pickle( )
{
    return py::pickle( []( const T& obj ) { return py::bytes( serializeToBinaryString( obj ) ); },
                       []( py::bytes data ) { return deserializeFromBinaryString< T >( data.cast< std::string >( ) ); } );
}

template< typename Base, typename Derived = Base >
auto make_pickle_polymorphic( )
{
    return py::pickle(
            []( const std::shared_ptr< Derived >& obj ) {
                // upcast to Base so cereal sees the polymorphic pointer
                return py::bytes( serializeSharedPtrToBinaryString< Base >( std::static_pointer_cast< Base >( obj ) ) );
            },
            []( py::bytes data ) {
                // cereal reconstructs the real derived type inside the shared_ptr<Base>
                // pybind11 then resolves the correct Python class from the dynamic type
                return std::static_pointer_cast< Derived >( deserializeSharedPtrFromBinaryString< Base >( data.cast< std::string >( ) ) );
            } );
}

}  // namespace tudat::serialization

// =====================================================================
//  Convenience macros for pybind11 exposure
//  Use __VA_ARGS__ so template types with commas (e.g. Foo<A,B>) work
// =====================================================================

//! Add __eq__ and __ne__ to a pybind11 class_ chain using operator==
#define TUDATPY_DEF_EQ_NE( ... )                                                                                                          \
    .def( "__eq__", &__VA_ARGS__::operator==, py::arg( "rhs" ) ).def( "__ne__", []( const __VA_ARGS__& self, const __VA_ARGS__& other ) { \
        return self != other;                                                                                                             \
    } )

//! Add pickle (value-type / standalone objects) via tudat::serialization::make_pickle
#define TUDATPY_DEF_PICKLE( ... ) .def( tudat::serialization::make_pickle< __VA_ARGS__ >( ) )

//! Add pickle (polymorphic, base-only) via tudat::serialization::make_pickle_polymorphic
#define TUDATPY_DEF_PICKLE_POLYMORPHIC( ... ) .def( tudat::serialization::make_pickle_polymorphic< __VA_ARGS__ >( ) )

//! Add pickle (polymorphic, with explicit base/derived) via tudat::serialization::make_pickle_polymorphic
#define TUDATPY_DEF_PICKLE_POLYMORPHIC_DERIVED( ... ) .def( tudat::serialization::make_pickle_polymorphic< __VA_ARGS__ >( ) )

// =====================================================================
//  File-IO exposure macros
//  These add save_to_json/load_from_json and/or save_to_binary/load_from_binary
//  to a pybind11 class_ chain.  __VA_ARGS__ handles template types.
// =====================================================================

//! Add save_to_json and load_from_json (value type — load returns T by value)
#define TUDATPY_DEF_JSON_IO( ... )                                      \
    .def( "save_to_json", &__VA_ARGS__::saveToJson, py::arg( "path" ) ) \
            .def_static(                                                \
                    "load_from_json", []( const std::string& path ) { return __VA_ARGS__::loadFromJson( path ); }, py::arg( "path" ) )

//! Add save_to_json and load_from_json (polymorphic — load returns shared_ptr<Base>)
#define TUDATPY_DEF_JSON_IO_POLYMORPHIC( ... )                          \
    .def( "save_to_json", &__VA_ARGS__::saveToJson, py::arg( "path" ) ) \
            .def_static(                                                \
                    "load_from_json", []( const std::string& path ) { return __VA_ARGS__::loadFromJson( path ); }, py::arg( "path" ) )

//! Add save_to_binary and load_from_binary (value type — load returns T by value)
#define TUDATPY_DEF_BINARY_IO( ... )                                        \
    .def( "save_to_binary", &__VA_ARGS__::saveToBinary, py::arg( "path" ) ) \
            .def_static(                                                    \
                    "load_from_binary", []( const std::string& path ) { return __VA_ARGS__::loadFromBinary( path ); }, py::arg( "path" ) )

//! Add save_to_binary and load_from_binary (polymorphic — load returns shared_ptr<Base>)
#define TUDATPY_DEF_BINARY_IO_POLYMORPHIC( ... )                            \
    .def( "save_to_binary", &__VA_ARGS__::saveToBinary, py::arg( "path" ) ) \
            .def_static(                                                    \
                    "load_from_binary", []( const std::string& path ) { return __VA_ARGS__::loadFromBinary( path ); }, py::arg( "path" ) )

//! Add both JSON and binary file IO (value type)
#define TUDATPY_DEF_FILE_IO( ... )     \
    TUDATPY_DEF_JSON_IO( __VA_ARGS__ ) \
    TUDATPY_DEF_BINARY_IO( __VA_ARGS__ )

//! Add both JSON and binary file IO (polymorphic)
#define TUDATPY_DEF_FILE_IO_POLYMORPHIC( ... )     \
    TUDATPY_DEF_JSON_IO_POLYMORPHIC( __VA_ARGS__ ) \
    TUDATPY_DEF_BINARY_IO_POLYMORPHIC( __VA_ARGS__ )
