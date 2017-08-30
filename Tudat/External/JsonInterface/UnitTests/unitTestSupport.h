/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_UNITTESTSUPPORT
#define TUDAT_JSONINTERFACE_UNITTESTSUPPORT

#include <boost/test/unit_test.hpp>

#include "Tudat/External/JsonInterface/Support/modular.h"
#include "Tudat/External/JsonInterface/Support/utilities.h"

namespace tudat
{

namespace json_interface
{

template< typename T = json >
T readInputFile( const std::string& filename, const std::string& extension = "json" )
{
    const path filePath = path( __FILE__ ).parent_path( ) / "inputs" / ( filename + "." + extension );
    boost::filesystem::current_path( filePath.parent_path( ) );
    return readJSON( filePath.string( ) ).get< T >( );
}


template< typename T >
void checkJsonEquivalent( const T& left, const T& right )
{
    const std::string fromFile = json( left ).dump( 2 );
    const std::string manual = json( right ).dump( 2 );
    BOOST_CHECK_EQUAL( fromFile, manual );
}

#define BOOST_CHECK_EQUAL_JSON( left, right ) tudat::json_interface::checkJsonEquivalent( left, right )


template< typename Enum >
void checkConsistentEnum( const std::string& filename,
                       const std::map< Enum, std::string >& stringValues,
                       const std::vector< Enum >& usupportedValues )
{
    using namespace json_interface;

    // Create vector of supported values
    std::vector< Enum > supportedValues;
    for ( auto entry : stringValues )
    {
        Enum value = entry.first;
        if ( ! contains( usupportedValues, value ) )
        {
            supportedValues.push_back( value );
        }
    }

    // Check that values and supportedValues are equivalent
    const std::vector< Enum > values = readInputFile< std::vector< Enum > >( filename );
    BOOST_CHECK_EQUAL_JSON( values, supportedValues );
}

#define BOOST_CHECK_EQUAL_ENUM( filename, stringValues, usupportedValues ) tudat::json_interface::checkConsistentEnum( filename, stringValues, usupportedValues )


} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_UNITTESTSUPPORT
