// serialization/pybind11_helpers.h
#pragma once
#include <pybind11/pybind11.h>
#include <stdexcept>
#include <string>
#include "tudat/io/serialization/file_io.h"

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

//! Polymorphic pickle for a BASE class — deserializes through the base pointer.
//! pybind11 resolves the correct Python class automatically from the shared_ptr's
//! dynamic type (std::type_info).
template< typename Base >
auto make_pickle_polymorphic_base( )
{
    return py::pickle(
            []( const Base& obj ) {
                // Non-owning shared_ptr from const_cast (similar to TUDAT_DEFINE_FILE_IO_POLYMORPHIC macros)
                auto ptr = std::shared_ptr< Base >( const_cast< Base* >( &obj ), []( Base* ) {} );
                return py::bytes( serializeSharedPtrToBinaryString< Base >( ptr ) );
            },
            []( py::bytes data ) {
                auto ptr = deserializeSharedPtrFromBinaryString< Base >( data.cast< std::string >( ) );
                if( !ptr )
                {
                    throw std::runtime_error(
                            "Polymorphic pickle deserialization returned nullptr. "
                            "The data may be corrupted or the type not registered." );
                }
                // pybind11 inspects the dynamic type of the shared_ptr and dispatches
                // to the correct Python derived class automatically
                return ptr;
            } );
}

//! Polymorphic pickle for a DERIVED class — deserializes through the base pointer,
//! then performs a checked dynamic_pointer_cast to the derived type.
template< typename Base, typename Derived >
auto make_pickle_polymorphic_derived( )
{
    return py::pickle(
            []( const Derived& obj ) {
                // Non-owning shared_ptr via const_cast, upcast to Base
                auto ptr = std::shared_ptr< Derived >( const_cast< Derived* >( &obj ), []( Derived* ) {} );
                return py::bytes( serializeSharedPtrToBinaryString< Base >( std::static_pointer_cast< Base >( ptr ) ) );
            },
            []( py::bytes data ) {
                // cereal reconstructs the real derived type inside the shared_ptr<Base>
                auto base = deserializeSharedPtrFromBinaryString< Base >( data.cast< std::string >( ) );
                if( !base )
                {
                    throw std::runtime_error(
                            "Polymorphic pickle deserialization returned nullptr. "
                            "The data may be corrupted or the type not registered." );
                }
                auto derived = std::dynamic_pointer_cast< Derived >( base );
                if( !derived )
                {
                    throw std::runtime_error(
                            "Polymorphic pickle downcast failed: the deserialized object "
                            "is not of the expected derived type." );
                }
                return derived;
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

//! Add pickle (polymorphic, BASE class) — deserializes through base pointer with
//! py::cast for dynamic type resolution
#define TUDATPY_DEF_PICKLE_POLYMORPHIC( ... ) .def( tudat::serialization::make_pickle_polymorphic_base< __VA_ARGS__ >( ) )

//! Add pickle (polymorphic, DERIVED class) — deserializes through base pointer,
//! then dynamic_pointer_cast to derived with null check
#define TUDATPY_DEF_PICKLE_POLYMORPHIC_DERIVED( ... ) .def( tudat::serialization::make_pickle_polymorphic_derived< __VA_ARGS__ >( ) )

// =====================================================================
//  File-IO exposure macros
//  These add save_to_json/load_from_json and/or save_to_binary/load_from_binary
//  to a pybind11 class_ chain.  __VA_ARGS__ handles template types.
//
//  POLYMORPHIC variants wrap the result in py::cast() so that the true
//  dynamic type is resolved to the correct Python class.
// =====================================================================

//! Add save_to_json and load_from_json (value type — load returns T by value)
#define TUDATPY_DEF_JSON_IO( ... )                                      \
    .def( "save_to_json", &__VA_ARGS__::saveToJson, py::arg( "path" ) ) \
            .def_static(                                                \
                    "load_from_json", []( const std::string& path ) { return __VA_ARGS__::loadFromJson( path ); }, py::arg( "path" ) )

//! Add save_to_json and load_from_json (polymorphic — load returns correct dynamic Python type)
#define TUDATPY_DEF_JSON_IO_POLYMORPHIC( ... )                                                                             \
    .def( "save_to_json", &__VA_ARGS__::saveToJson, py::arg( "path" ) )                                                    \
            .def_static(                                                                                                   \
                    "load_from_json",                                                                                      \
                    []( const std::string& path ) -> py::object { return py::cast( __VA_ARGS__::loadFromJson( path ) ); }, \
                    py::arg( "path" ) )

//! Add save_to_binary and load_from_binary (value type — load returns T by value)
#define TUDATPY_DEF_BINARY_IO( ... )                                        \
    .def( "save_to_binary", &__VA_ARGS__::saveToBinary, py::arg( "path" ) ) \
            .def_static(                                                    \
                    "load_from_binary", []( const std::string& path ) { return __VA_ARGS__::loadFromBinary( path ); }, py::arg( "path" ) )

//! Add save_to_binary and load_from_binary (polymorphic — load returns correct dynamic Python type)
#define TUDATPY_DEF_BINARY_IO_POLYMORPHIC( ... )                                                                             \
    .def( "save_to_binary", &__VA_ARGS__::saveToBinary, py::arg( "path" ) )                                                  \
            .def_static(                                                                                                     \
                    "load_from_binary",                                                                                      \
                    []( const std::string& path ) -> py::object { return py::cast( __VA_ARGS__::loadFromBinary( path ) ); }, \
                    py::arg( "path" ) )

//! Add both JSON and binary file IO (value type)
#define TUDATPY_DEF_FILE_IO( ... )     \
    TUDATPY_DEF_JSON_IO( __VA_ARGS__ ) \
    TUDATPY_DEF_BINARY_IO( __VA_ARGS__ )

//! Add both JSON and binary file IO (polymorphic — load returns correct dynamic Python type)
#define TUDATPY_DEF_FILE_IO_POLYMORPHIC( ... )     \
    TUDATPY_DEF_JSON_IO_POLYMORPHIC( __VA_ARGS__ ) \
    TUDATPY_DEF_BINARY_IO_POLYMORPHIC( __VA_ARGS__ )
